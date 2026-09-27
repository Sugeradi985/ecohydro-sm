"""形变时序分解：把储变量信号从永久压实里分离出来。

为什么必须先分解
----------------
InSAR 时序形变 = 长期非弹性压实（不可逆，不代表储变量）
              + 季节性弹性变形（可恢复，才是储变量的载体）
              + 噪声

把总形变直接当储变量用，会把永久压实误算成可恢复储水，量级错误可达数倍。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

__all__ = [
    "DeformationConfig",
    "harmonic_trend_model",
    "fit_decomposition",
    "estimate_lag_days",
    "seasonal_amplitude_snr",
]


@dataclass
class DeformationConfig:
    """时序分解参数。"""

    periods_days: tuple[float, ...] = (365.25, 182.6)  # 年周期 + 半年周期
    trend_order: int = 1                               # 1 = 线性趋势，2 = 含加速度
    max_lag_days: int = 120                            # 时滞互相关的最大搜索延迟
    min_snr: float = 3.0                               # 季节振幅信噪比下限


def _design_matrix(t_days: np.ndarray, cfg: DeformationConfig) -> np.ndarray:
    t = np.asarray(t_days, dtype=float)
    cols = [np.ones_like(t)]
    for p in cfg.periods_days:
        cols += [np.sin(2 * np.pi * t / p), np.cos(2 * np.pi * t / p)]
    for o in range(1, cfg.trend_order + 1):
        cols.append((t / 365.25) ** o)
    return np.column_stack(cols)


def harmonic_trend_model(t_days: np.ndarray, params: np.ndarray, cfg: DeformationConfig) -> np.ndarray:
    """由参数重建"谐波 + 趋势"模型。"""
    return _design_matrix(t_days, cfg) @ params


def fit_decomposition(
    t_days: np.ndarray,
    disp_mm: np.ndarray,
    cfg: DeformationConfig | None = None,
):
    """逐像元拟合"谐波 + 趋势"，并分离季节项与趋势项。

    参数
    ----
    t_days : (n_time,) 相对天数
    disp_mm : (n_time, n_pixels) 形变序列（mm）

    返回
    ----
    dict，含 seasonal（季节项）、trend（趋势项）、residual、
    amplitude（年周期振幅）、params
    """
    cfg = cfg or DeformationConfig()
    t = np.asarray(t_days, dtype=float)
    d = np.asarray(disp_mm, dtype=float)
    n_t, n_px = d.shape
    A = _design_matrix(t, cfg)

    params = np.full((A.shape[1], n_px), np.nan)
    for i in range(n_px):
        y = d[:, i]
        m = np.isfinite(y)
        if m.sum() < A.shape[1] + 2:
            continue
        params[:, i], *_ = np.linalg.lstsq(A[m], y[m], rcond=None)

    model = A @ params
    seasonal = (A[:, 1:1 + 2 * len(cfg.periods_days)]) @ params[1:1 + 2 * len(cfg.periods_days)]
    trend = A[:, 1 + 2 * len(cfg.periods_days):] @ params[1 + 2 * len(cfg.periods_days):]
    amplitude = np.sqrt(params[1] ** 2 + params[2] ** 2)  # 年周期振幅
    return {
        "seasonal": seasonal,
        "trend": trend,
        "intercept": params[0],
        "residual": d - model,
        "amplitude_mm": amplitude,
        "params": params,
    }


def estimate_lag_days(
    seasonal_disp: np.ndarray,
    head_series: np.ndarray,
    t_days: np.ndarray,
    cfg: DeformationConfig | None = None,
) -> int:
    """用时滞互相关估计形变相对水位的滞后天数。

    弱透水层排水延迟会让季节形变滞后于水位，不校正会系统性低估相关性。

    参数
    ----
    seasonal_disp : (n_time, n_pixels) 季节形变
    head_series : (n_time,) 水位序列（区域平均）
    """
    cfg = cfg or DeformationConfig()
    t = np.asarray(t_days, dtype=float)
    h = np.asarray(head_series, dtype=float).ravel()
    sd = np.nanmean(np.asarray(seasonal_disp, dtype=float), axis=1)

    dt = float(np.median(np.diff(np.sort(t)))) if len(t) > 1 else 1.0
    max_shift = int(min(cfg.max_lag_days, len(t) // 3) / max(dt, 1e-9))

    best_lag, best_abs_corr = 0, -1.0
    for lag in range(0, max_shift + 1):
        if lag == 0:
            a, b = sd, h
        else:
            # 形变滞后于水位 L 个样本：sd[i + L] 应与 h[i] 最相关
            a, b = sd[lag:], h[:-lag]
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 10:
            continue
        c = np.corrcoef(a[m], b[m])[0, 1]
        if np.isfinite(c) and abs(c) > best_abs_corr:
            best_abs_corr, best_lag = abs(c), lag
    return int(best_lag * dt)


def seasonal_amplitude_snr(decomp: dict, cfg: DeformationConfig | None = None) -> np.ndarray:
    """季节振幅信噪比：振幅 / 残差标准差。

    **这是 P0 预检的判据**：SNR >= 3 才说明该区具备用形变约束储变量的条件。
    """
    cfg = cfg or DeformationConfig()
    resid_std = np.nanstd(decomp["residual"], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = decomp["amplitude_mm"] / np.where(resid_std <= 0, np.nan, resid_std)
    return snr
