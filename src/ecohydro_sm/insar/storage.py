"""形变 -> 含水层储变量 ΔS。

核心关系
--------
    Δd   = S_ke · Δh_w      （弹性形变与水位变化的关系）
    ΔS   = S_y  · Δh_w      （储变量与水位变化的关系）
    ⇒ ΔS = (S_y / S_ke) · Δd

两条路径
--------
路径 A（有井 + 有抽水量或 GRACE 约束）：直接由上式求绝对储变量。
路径 B（井稀疏）：**形变不作为绝对量来源**，只作为水位场的空间加密观测，
用 ΔS = S_y · Δh_w 求储变量，形变用于约束 Δh_w 的空间插值。
路径 B 更稳健、更可辩护，是默认推荐。

重要警告
--------
Δd 不等于 ΔS。S_ke 在华北平原文献中跨越 2×10⁻⁴ – 2.1×10⁻²，跨两个数量级。
**结果必须给区间并做敏感性分析，禁止单点估计。**
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "StorageConfig",
    "calibrate_ske",
    "storage_change_from_deformation",
    "storage_change_from_head",
    "sy_sensitivity_span",
]


@dataclass
class StorageConfig:
    """储变量转换参数。"""

    ske_default: float = 5.0e-3          # 弹性骨架储水系数（无量纲）
    sy_default: float = 0.10             # 给水度（无量纲）
    ske_bounds: tuple[float, float] = (2.0e-4, 2.1e-2)
    sy_bounds: tuple[float, float] = (0.02, 0.30)
    n_sensitivity_samples: int = 200
    unit_mm_per_m: float = 1000.0        # 1 m 水位 = 1000 mm
    random_state: int = 42


def calibrate_ske(
    seasonal_disp_mm: np.ndarray,
    head_change_m: np.ndarray,
    cfg: StorageConfig | None = None,
) -> float:
    """由季节形变幅度与水位幅度之比标定 S_ke。

        S_ke = Δd / Δh_w  （Δh_w 需换算为 mm）

    参数
    ----
    seasonal_disp_mm : 季节形变（mm）
    head_change_m : 同期水位变化（m，正为上升）
    """
    cfg = cfg or StorageConfig()
    d = np.asarray(seasonal_disp_mm, dtype=float)
    h = np.asarray(head_change_m, dtype=float) * cfg.unit_mm_per_m
    # 水位为单点序列而形变为栅格时，按时间轴广播
    if h.ndim == 1 and d.ndim == 2 and h.shape[0] == d.shape[0]:
        h = np.repeat(h[:, None], d.shape[1], axis=1)
    d, h = np.broadcast_arrays(d, h)
    d, h = d.ravel(), h.ravel()
    m = np.isfinite(d) & np.isfinite(h) & (np.abs(h) > 1e-6)
    if m.sum() < 3:
        raise ValueError("有效样本不足，无法标定 S_ke（至少需要 3 个）。")
    # 过原点的最小二乘：Δd = S_ke · Δh_w
    # 符号约定：水位下降（Δh_w < 0）对应沉降（Δd < 0），故 S_ke > 0。
    ske = float(np.sum(d[m] * h[m]) / np.sum(h[m] ** 2))
    lo, hi = cfg.ske_bounds
    return float(np.clip(ske, lo, hi))


def storage_change_from_deformation(
    seasonal_disp_mm: np.ndarray,
    ske: float,
    sy: float,
) -> np.ndarray:
    """路径 A：ΔS = (S_y / S_ke) · Δd。

    参数
    ----
    seasonal_disp_mm : 季节弹性形变（mm，正为抬升 = 储水增加）
    """
    if ske <= 0:
        raise ValueError("S_ke 必须为正。")
    return np.asarray(seasonal_disp_mm, dtype=float) * (sy / ske)


def storage_change_from_head(head_change_m: np.ndarray, sy: float) -> np.ndarray:
    """路径 B：ΔS = S_y · Δh_w，单位由 m 水位换算为 mm 水当量。"""
    return np.asarray(head_change_m, dtype=float) * sy * 1000.0


def sy_sensitivity_span(
    seasonal_disp_mm: np.ndarray,
    ske: float,
    cfg: StorageConfig | None = None,
):
    """对 S_y 做区间敏感性分析，返回 ΔS 的中位数与 5–95% 区间。

    参数不确定性远大于像元噪声，给区间是**不可省略的报告要求**。
    """
    cfg = cfg or StorageConfig()
    rng = np.random.default_rng(cfg.random_state)
    lo, hi = cfg.sy_bounds
    sy_samples = rng.uniform(lo, hi, cfg.n_sensitivity_samples)
    deltas = [storage_change_from_deformation(seasonal_disp_mm, ske, s) for s in sy_samples]
    stack = np.stack(deltas, axis=0)
    return {
        "median_mm": np.nanmedian(stack, axis=0),
        "p05_mm": np.nanpercentile(stack, 5, axis=0),
        "p95_mm": np.nanpercentile(stack, 95, axis=0),
        "sy_samples": sy_samples,
    }
