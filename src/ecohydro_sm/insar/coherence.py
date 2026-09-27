"""InSAR 相干性支路：把"植被结构去相干"从"含水量去相干"里剥出去。

这是本仓库相对 INSAR4SM 的改进点
--------------------------------
INSAR4SM（Karamvasis & Karathanassi 2023）用气象信息找最干时刻、
再按时序排序反演土壤水，干旱区精度很好（RMSE 0.029 m3/m3，R 0.78），
但**生长季植被结构快速变化时误差显著**——因为相干性衰减里混了结构变化。

本模块把去相干建模为两个指数项的乘积：

    γ(t) = γ_dry · exp(−k_sm · ΔSM(t)) · exp(−k_v · ΔV(t))

其中 ΔV(t) 由 Sentinel-2 的 FVC/NDVI 变化率构建。两个系数 k_sm、k_v
用区域尺度同时标定，再逐像元反演 ΔSM。

适用域：VWC < 2 kg/m2 的裸土、稀疏植被与干旱半干旱区。
农田绿洲内部去相关严重，不承诺精度，应打标降级。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

__all__ = [
    "CoherenceSMConfig",
    "structural_change_index",
    "calibrate_decay_coefficients",
    "invert_sm_from_coherence",
    "closure_phase_residual",
]


@dataclass
class CoherenceSMConfig:
    """相干性土壤水反演参数。"""

    dry_reference_quantile: float = 0.95  # 相干性最高（最干）的分位
    min_gamma: float = 0.15               # 低于此值视为完全去相关，不参与反演
    max_gamma: float = 1.0
    k_sm_init: float = 2.0                # 湿度衰减系数初值
    k_v_init: float = 1.0                 # 结构衰减系数初值
    sm_min: float = 0.0
    sm_max: float = 0.5


def structural_change_index(fvc_t: np.ndarray, baseline: float | None = None) -> np.ndarray:
    """由 Sentinel-2 的 FVC 时序构建结构变化指数 ΔV(t)。

    定义 ΔV(t) = |FVC(t) − FVC_baseline|，baseline 缺省取时序中位数
    （代表该像元的常态冠层结构）。

    参数
    ----
    fvc_t : (n_time, ...) 植被覆盖度时序
    """
    f = np.asarray(fvc_t, dtype=float)
    base = np.nanmedian(f, axis=0) if baseline is None else baseline
    dv = np.abs(f - base)
    return np.nan_to_num(dv, nan=0.0)


def calibrate_decay_coefficients(
    gamma: np.ndarray,
    delta_v: np.ndarray,
    delta_sm_ref: np.ndarray,
    cfg: CoherenceSMConfig | None = None,
) -> dict:
    """用参考土壤水时序同时标定 k_sm、k_v 与归一化偏置。

    模型：ln(γ / γ_dry_ref) = c − k_sm·ΔSM − k_v·ΔV

    其中 c 是**显式拟合的截距**，用于吸收"以时序最大值做归一化"带来的系统性
    偏置（噪声最大值会把基准抬高，使 k 估计偏低）。不做这一步会让 k 低估 20% 量级。

    参数
    ----
    gamma : (n_time, n_pixels) 相干性时序
    delta_v : (n_time, n_pixels) 结构变化指数
    delta_sm_ref : (n_time, n_pixels) 参考含水量变化（ERA5-Land 或站点）

    返回
    ----
    dict：k_sm、k_v、log_offset(c)
    """
    cfg = cfg or CoherenceSMConfig()
    g = np.clip(np.asarray(gamma, dtype=float), cfg.min_gamma, cfg.max_gamma)
    dv = np.asarray(delta_v, dtype=float)
    dsm = np.asarray(delta_sm_ref, dtype=float)

    g0 = np.nanmax(g, axis=0)  # 以最干（相干性最高）时刻为基准
    g0 = np.where(g0 <= 0, np.nan, g0)
    ratio = g / g0

    mask = np.isfinite(ratio) & np.isfinite(dv) & np.isfinite(dsm)
    r = np.log(np.clip(ratio[mask].ravel(), 1e-3, None))
    v = dv[mask].ravel()
    s = dsm[mask].ravel()

    # p = [c, k_sm, k_v]
    def resid(p):
        return r - (p[0] - p[1] * s - p[2] * v)

    sol = least_squares(resid, x0=[0.0, cfg.k_sm_init, cfg.k_v_init],
                        bounds=([-10.0, 0.0, 0.0], [10.0, 50.0, 50.0]))
    return {"log_offset": float(sol.x[0]), "k_sm": float(sol.x[1]), "k_v": float(sol.x[2])}


def invert_sm_from_coherence(
    gamma: np.ndarray,
    delta_v: np.ndarray,
    k_sm: float,
    k_v: float,
    sm_dry: float = 0.02,
    log_offset: float = 0.0,
    cfg: CoherenceSMConfig | None = None,
) -> np.ndarray:
    """由相干性反演绝对体积含水量。

    模型：ln(γ/γ_dry) = c − k_sm·ΔSM − k_v·ΔV
    反解：ΔSM(t) = [c − ln(γ/γ_dry) − k_v·ΔV(t)] / k_sm
    合成：SM(t)   = SM_dry + ΔSM(t)

    注意符号：相干性随含水量升高而**下降**，故 ΔSM 前为负号；
    合成时是 SM_dry **加** ΔSM。弄反会让结果整体偏低约 2·ΔSM。
    """
    cfg = cfg or CoherenceSMConfig()
    g = np.clip(np.asarray(gamma, dtype=float), cfg.min_gamma, cfg.max_gamma)
    dv = np.asarray(delta_v, dtype=float)

    g_dry = np.nanmax(g, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.clip(g / np.where(g_dry <= 0, np.nan, g_dry), 1e-3, 1.0)
        dsm = (log_offset - np.log(ratio) - k_v * dv) / max(k_sm, 1e-6)
    sm = sm_dry + dsm
    return np.clip(sm, cfg.sm_min, cfg.sm_max)


def closure_phase_residual(phase_triplets: np.ndarray) -> np.ndarray:
    """闭合相位残差：φ_12 + φ_23 + φ_31 理论上为 0。

    非零残差反映非对称散射体变化（介电变化 + 结构变化）。
    残差大的像元说明结构主导，相干性反演不可用 —— 这是**不依赖地面真值的
    独立质量控制量**（De Zan & Gomba, 2018）。

    参数
    ----
    phase_triplets : (n_triplets, 3, ...) 三个干涉图的相位（弧度）

    返回
    ----
    每个像元的闭合相位残差标准差（弧度）
    """
    p = np.asarray(phase_triplets, dtype=float)
    if p.ndim < 3 or p.shape[1] != 3:
        raise ValueError("phase_triplets 形状应为 (n_triplets, 3, ...)。")
    closure = p[:, 0] + p[:, 1] + p[:, 2]
    wrapped = np.angle(np.exp(1j * closure))  # 相位缠绕到 (−π, π]
    return np.nanstd(wrapped, axis=0)
