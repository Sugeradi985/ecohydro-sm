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
    "PairCoherenceSMConfig",
    "calibrate_decay_coefficients_pairwise",
    "invert_delta_sm_from_pairs",
    "assign_pair_signs",
    "integrate_pair_deltas",
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


# ---------------------------------------------------------------------------
# 干涉对模式（pairwise）
# ---------------------------------------------------------------------------
# 上面的 γ(t) 模型假设输入是"每个时相相对干参考的相干性"。
# ASF HyP3 的 INSAR_GAMMA 产品给的是**逐干涉对**的 γ_ij（每个 pair 一个
# coherence GeoTIFF），两者不是同一个量。直接把 pair 堆成 (n_time, ...)
# 会得到一个语义错误的"时序"，据此反演的 SM 没有物理意义。
#
# 对模式直接对 pair 建模：
#     ln γ_ij = c − k_sm·|ΔSM_ij| − k_v·|ΔV_ij| − k_t·Δt_ij
# Δt 项不可省：star 网络或长基线对里，时间去相干与含水量去相干无法区分。
# c 吸收"理想相干性 = 1"的基准，等价于 date_series 模式的显式截距。
#
# 反演只给出 |ΔSM|（幅度），符号必须由外部参考（站点 ΔSM / L1 / ERA5-Land）
# 提供，再用最小二乘从干锚点积分成逐日期序列。见 assign_pair_signs 与
# integrate_pair_deltas。


@dataclass
class PairCoherenceSMConfig:
    """干涉对模式参数。"""

    min_gamma: float = 0.15
    max_gamma: float = 1.0
    k_sm_init: float = 2.0
    k_v_init: float = 1.0
    k_t_init: float = 0.002      # 每天的时间去相干系数初值
    k_t_max: float = 0.2
    sm_min: float = 0.0
    sm_max: float = 0.5


def _pair_time_baseline(dt_days: np.ndarray, shape: tuple) -> np.ndarray:
    """把 (n_pairs,) 的时间基线广播成与相干性数组一致的形状（支持 1D/2D）。"""
    d = np.abs(np.asarray(dt_days, dtype=float)).reshape(-1)
    d = d.reshape((-1,) + (1,) * (len(shape) - 1))
    return np.broadcast_to(d, shape).astype(float)


def calibrate_decay_coefficients_pairwise(
    gamma_pairs: np.ndarray,
    delta_v_pairs: np.ndarray,
    delta_sm_ref_pairs: np.ndarray,
    dt_days: np.ndarray,
    cfg: PairCoherenceSMConfig | None = None,
) -> dict:
    """用干涉对标定 k_sm、k_v、k_t 与截距 c。

    模型：ln γ_ij = c − k_sm·|ΔSM_ij| − k_v·|ΔV_ij| − k_t·Δt_ij

    参数
    ----
    gamma_pairs : (n_pairs, n_pixels) 每个干涉对的相干性 0–1
    delta_v_pairs : (n_pairs, n_pixels) 该对两端的结构变化量差（S-2 FVC）
    delta_sm_ref_pairs : (n_pairs, n_pixels) 该对两端的参考含水量差
    dt_days : (n_pairs,) 时间基线（天）

    返回
    ----
    dict：log_offset(c)、k_sm、k_v、k_t、n_samples、rmse
    """
    cfg = cfg or PairCoherenceSMConfig()
    g = np.clip(np.asarray(gamma_pairs, dtype=float), cfg.min_gamma, cfg.max_gamma)
    dv = np.abs(np.asarray(delta_v_pairs, dtype=float))
    dsm = np.abs(np.asarray(delta_sm_ref_pairs, dtype=float))
    dt = _pair_time_baseline(dt_days, g.shape)

    mask = np.isfinite(g) & np.isfinite(dv) & np.isfinite(dsm) & np.isfinite(dt)
    r = np.log(g[mask].ravel())
    s = dsm[mask].ravel()
    v = dv[mask].ravel()
    t = dt[mask].ravel()
    if r.size < 10:
        raise ValueError("有效样本不足 10 个，无法标定干涉对模型。")

    # p = [c, k_sm, k_v, k_t]
    def resid(p):
        return r - (p[0] - p[1] * s - p[2] * v - p[3] * t)

    sol = least_squares(
        resid,
        x0=[0.0, cfg.k_sm_init, cfg.k_v_init, cfg.k_t_init],
        bounds=([-10.0, 0.0, 0.0, 0.0], [10.0, 50.0, 50.0, cfg.k_t_max]),
    )
    c, k_sm, k_v, k_t = (float(x) for x in sol.x)
    return {
        "log_offset": c,
        "k_sm": k_sm,
        "k_v": k_v,
        "k_t": k_t,
        "n_samples": int(r.size),
        "rmse": float(np.sqrt(np.mean(sol.fun ** 2))),
    }


def invert_delta_sm_from_pairs(
    gamma_pairs: np.ndarray,
    delta_v_pairs: np.ndarray,
    dt_days: np.ndarray,
    coef: dict,
    cfg: PairCoherenceSMConfig | None = None,
) -> np.ndarray:
    """由干涉对相干性反演 |ΔSM|（每个 pair 一个值）。

    |ΔSM_ij| = [c − ln γ_ij − k_v·|ΔV_ij| − k_t·Δt_ij] / k_sm

    注意：只给幅度，不给方向。符号必须用 assign_pair_signs 从外部参考确定。
    """
    cfg = cfg or PairCoherenceSMConfig()
    g = np.clip(np.asarray(gamma_pairs, dtype=float), cfg.min_gamma, cfg.max_gamma)
    dv = np.abs(np.asarray(delta_v_pairs, dtype=float))
    dt = _pair_time_baseline(dt_days, g.shape)

    k_sm = max(float(coef["k_sm"]), 1e-6)
    k_v = float(coef.get("k_v", 0.0))
    k_t = float(coef.get("k_t", 0.0))
    c = float(coef.get("log_offset", 0.0))

    with np.errstate(divide="ignore", invalid="ignore"):
        dsm = (c - np.log(g) - k_v * dv - k_t * dt) / k_sm
    return np.clip(np.nan_to_num(dsm, nan=0.0), 0.0, None)


def assign_pair_signs(delta_sm_ref_pairs: np.ndarray) -> np.ndarray:
    """由参考序列（站点 / L1 / ERA5-Land）确定每个 pair 的 ΔSM 符号。

    返回 ±1 数组；参考值恰好为 0 时默认 +1。
    没有可信参考时用本函数等于假设"所有 pair 都在变湿"，会得到单调上升的
    SM 序列——这不是可接受的替代方案，务必提供真实参考。
    """
    d = np.asarray(delta_sm_ref_pairs, dtype=float)
    signs = np.sign(d)
    signs[signs == 0] = 1.0
    return signs


def integrate_pair_deltas(
    delta_sm_pairs: np.ndarray,
    ref_idx: np.ndarray,
    sec_idx: np.ndarray,
    n_dates: int,
    signs: np.ndarray | None = None,
    anchor_date: int = 0,
    anchor_sm: float = 0.02,
    anchor_weight: float = 10.0,
    sm_bounds: tuple[float, float] = (0.0, 0.5),
) -> np.ndarray:
    """把逐 pair 的 |ΔSM| 积分成逐日期 SM 序列。

    对每条 pair 建方程 SM[sec] − SM[ref] = sign · |ΔSM|，
    再加一条锚点方程 SM[anchor] = anchor_sm，用最小二乘求解。
    若 pair 网络不连通，未约束的日期会由锚点方程的正则化主导，应检查
    连通性后再使用结果。

    参数
    ----
    delta_sm_pairs : (n_pairs, n_pixels) 或 (n_pairs,) 反演出的 |ΔSM|
    ref_idx, sec_idx : (n_pairs,) 参考/副影像在日期轴上的索引
    n_dates : int 日期总数
    signs : (n_pairs, n_pixels) 或 (n_pairs,) 符号，缺省全 +1（不推荐）

    返回
    ----
    (n_dates, ...) 逐日期体积含水量
    """
    d = np.asarray(delta_sm_pairs, dtype=float)
    scalar = d.ndim == 1
    if scalar:
        d = d[:, None]
    ref = np.asarray(ref_idx, dtype=int).ravel()
    sec = np.asarray(sec_idx, dtype=int).ravel()
    n_pairs = ref.size
    if not (d.shape[0] == n_pairs == sec.size):
        raise ValueError("delta_sm_pairs、ref_idx、sec_idx 长度不一致。")

    if signs is None:
        sg = np.ones_like(d)
    else:
        s = np.asarray(signs, dtype=float)
        if scalar and s.ndim == 1:
            s = s[:, None]
        sg = np.broadcast_to(s, d.shape)
    b_eq = sg * d

    flat_y = d.shape[1]
    out = np.empty((n_dates, flat_y), dtype=float)
    for j in range(flat_y):
        # A: (n_pairs + 1, n_dates)
        a = np.zeros((n_pairs + 1, n_dates), dtype=float)
        rows = np.arange(n_pairs)
        a[rows, ref] = -1.0
        a[rows, sec] = 1.0
        b = np.concatenate([b_eq[:, j], [0.0]])
        a[-1, anchor_date] = anchor_weight
        b[-1] = anchor_weight * anchor_sm
        sol, *_ = np.linalg.lstsq(a, b, rcond=None)
        out[:, j] = sol

    lo, hi = sm_bounds
    res = np.clip(out, lo, hi)
    return res[:, 0] if scalar else res
