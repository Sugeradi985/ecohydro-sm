"""SAR 斑点滤波。

注意：斑点滤波会平滑相干性估计所用的相位统计特性。
**L1（σ0 反演）支路必须滤波；L2（相干性）支路在 SLC 上单独估计，不要复用被滤波的幅度。**
"""

from __future__ import annotations

import numpy as np

__all__ = ["lee_sigma", "refined_lee", "boxcar", "cv_equivalent_number_of_looks"]


def _local_stats(img: np.ndarray, window: int):
    """用积分图计算局部均值与方差。"""
    img = np.asarray(img, dtype=float)
    pad = window // 2
    padded = np.pad(img, pad, mode="reflect")
    ker = np.ones((window, window))
    from numpy.lib.stride_tricks import sliding_window_view

    views = sliding_window_view(padded, (window, window))
    mean = np.einsum("ij,...ij->...", ker, views) / (window * window)
    var = np.einsum("ij,...ij->...", ker, (views - mean[..., None, None]) ** 2) / (window * window)
    return mean, var


def lee_sigma(img: np.ndarray, window: int = 7, sigma_range: tuple[float, float] = (0.9, 1.1)) -> np.ndarray:
    """Lee-Sigma 滤波：只在同质区做均值平滑，保留边缘与点目标。

    参数
    ----
    img : 线性域后向散射（不是 dB）
    """
    img = np.asarray(img, dtype=float)
    mean, _ = _local_stats(img, window)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = img / mean
    lo, hi = sigma_range
    keep = (ratio >= lo) & (ratio <= hi)
    out = np.where(keep, mean, img)
    return out


def refined_lee(img: np.ndarray, window: int = 7, n_looks: float = 4.0) -> np.ndarray:
    """Refined Lee 滤波：按局部变异系数划分同质/异质区，分别处理。

    参数
    ----
    img : 线性域后向散射
    n_looks : 等效视数（用于噪声方差阈值）
    """
    img = np.asarray(img, dtype=float)
    mean, var = _local_stats(img, window)
    with np.errstate(divide="ignore", invalid="ignore"):
        ci = np.sqrt(np.maximum(var, 0.0)) / mean       # 变异系数
        cu = 1.0 / np.sqrt(max(n_looks, 1.0))           # 均匀区理论值
        cmax = np.sqrt(2.0) * cu                        # 异质区阈值

    # 同质区：局部均值；异质区：反距离加权；点目标：原值
    w = np.clip((cmax - ci) / (cmax - cu), 0.0, 1.0)
    out = img * (1.0 - w) + mean * w
    out = np.where(ci <= cu, mean, out)
    out = np.where(ci >= cmax, img, out)
    return out


def boxcar(img: np.ndarray, window: int = 5) -> np.ndarray:
    """多视平均（boxcar）。最朴素但会明显降低分辨率，仅作对照。"""
    mean, _ = _local_stats(np.asarray(img, dtype=float), window)
    return mean


def cv_equivalent_number_of_looks(img: np.ndarray, window: int = 7) -> np.ndarray:
    """由局部变异系数估计等效视数 ENL = 1/CV²，用于评估滤波效果与残余斑点。"""
    mean, var = _local_stats(np.asarray(img, dtype=float), window)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv2 = var / mean**2
    return np.where(cv2 > 0, 1.0 / cv2, np.nan)
