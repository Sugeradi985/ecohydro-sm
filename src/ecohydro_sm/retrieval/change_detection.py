"""变化检测法表层土壤水（Wagner / Brocca 路线）。

原理：在粗糙度与植被结构短时不变的假设下，σ0 的时间变化归因于含水量变化，
用干、湿参考值把 σ0 线性缩放到体积含水量。

关键前提（也是最大的隐性误差源）
--------------------------------
"粗糙度时不变"会被耕作、土壤结皮、作物残茬、物候变化系统性破坏。
因此：
1. 参考值必须按**物候窗口**分别取（用 S-2 的 FVC/NDVI 分层），不能全年统一；
2. 必须显式给出 σ0 的动态范围（dry/wet 参考之差），范围过小的像元直接标记为不可反演。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ChangeDetectionConfig", "dry_wet_references", "change_detection_sm", "dynamic_range_flag"]


@dataclass
class ChangeDetectionConfig:
    """变化检测参数。"""

    mv_dry: float = 0.02     # 干参考对应的体积含水量（需本地标定）
    mv_wet: float = 0.45     # 湿参考对应的体积含水量（需本地标定）
    min_dynamic_range_db: float = 3.0  # σ0 动态范围下限，低于此值判为不可反演
    wet_percentile: float = 95.0
    dry_percentile: float = 5.0


def dry_wet_references(
    sigma0_db: np.ndarray,
    cfg: ChangeDetectionConfig | None = None,
    axis: int = 0,
    valid: np.ndarray | None = None,
):
    """沿时间轴取干、湿参考值。

    参数
    ----
    sigma0_db : (n_time, ...) 后向散射 dB
    valid : 可选掩膜，与 sigma0_db 同形，True 参与统计

    返回
    ----
    dry, wet : 干/湿参考（dB），形状为去掉时间轴后的形状
    """
    cfg = cfg or ChangeDetectionConfig()
    arr = np.asarray(sigma0_db, dtype=float)
    if valid is not None:
        arr = np.where(np.asarray(valid, dtype=bool), arr, np.nan)
    dry = np.nanpercentile(arr, cfg.dry_percentile, axis=axis)
    wet = np.nanpercentile(arr, cfg.wet_percentile, axis=axis)
    return dry, wet


def change_detection_sm(
    sigma0_db: np.ndarray,
    dry: np.ndarray,
    wet: np.ndarray,
    cfg: ChangeDetectionConfig | None = None,
) -> np.ndarray:
    """把 σ0 线性缩放到体积含水量。

        SM = mv_dry + (σ0 − σ0_dry) / (σ0_wet − σ0_dry) · (mv_wet − mv_dry)
    """
    cfg = cfg or ChangeDetectionConfig()
    s = np.asarray(sigma0_db, dtype=float)
    d = np.asarray(dry, dtype=float)
    w = np.asarray(wet, dtype=float)
    denom = w - d
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (s - d) / np.where(np.abs(denom) < 1e-9, np.nan, denom)
    sm = cfg.mv_dry + frac * (cfg.mv_wet - cfg.mv_dry)
    return np.clip(sm, 0.0, 0.6)


def dynamic_range_flag(dry: np.ndarray, wet: np.ndarray, cfg: ChangeDetectionConfig | None = None) -> np.ndarray:
    """动态范围质检：干湿差过小的像元无法可靠反演，返回 True = 可反演。"""
    cfg = cfg or ChangeDetectionConfig()
    rng = np.asarray(wet, dtype=float) - np.asarray(dry, dtype=float)
    return np.isfinite(rng) & (rng >= cfg.min_dynamic_range_db)
