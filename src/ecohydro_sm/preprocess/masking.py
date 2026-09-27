"""质量控制掩膜：云/雪、入射角、冻融、降水、适用域分区。

掩膜是产品质量的第一道闸门。所有掩膜返回布尔数组，True = 保留。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "MaskConfig",
    "s2_cloud_mask",
    "incidence_mask",
    "freeze_thaw_mask",
    "precipitation_mask",
    "vwc_domain_mask",
    "combine_masks",
]


@dataclass
class MaskConfig:
    """掩膜阈值。默认值是保守起点，按研究区调整。"""

    theta_min: float = 25.0
    theta_max: float = 45.0
    air_temp_threshold_c: float = 1.0   # 2 m 气温低于此值判为潜在冻结
    soil_temp_threshold_c: float = 0.0  # 表层土壤温度低于此值判为冻结
    precip_threshold_mm: float = 1.0    # 当日降水超过此值需打标
    precip_lag_days: int = 1           # 降水后延迟打标天数
    vwc_limit_kg_m2: float = 2.0       # 相干性支路的适用上限
    vwc_flag_kg_m2: float = 3.0        # L1 产品精度退化需打标的上限


# Sentinel-2 Scene Classification Layer 中需要剔除的类别
SCL_INVALID = {1: "saturated/defective", 3: "cloud shadow", 8: "cloud medium prob",
               9: "cloud high prob", 10: "thin cirrus", 11: "snow/ice"}


def s2_cloud_mask(scl: np.ndarray, buffer_pixels: int = 0) -> np.ndarray:
    """由 SCL 波段生成有效观测掩膜（True = 保留）。

    参数
    ----
    scl : Sentinel-2 SCL 波段
    buffer_pixels : 云边缘缓冲（像元）。云影常被低估，建议 >= 1。
    """
    scl = np.asarray(scl)
    keep = np.isin(scl, list(SCL_INVALID), invert=True)
    keep &= np.isfinite(scl)
    if buffer_pixels > 0:
        keep = _erode(keep, buffer_pixels)
    return keep


def _erode(mask: np.ndarray, n: int) -> np.ndarray:
    """简单的 n 像元腐蚀（不依赖 scipy.ndimage，避免额外依赖）。"""
    out = mask.copy()
    for _ in range(n):
        padded = np.pad(out, 1, mode="constant", constant_values=False)
        neigh = (
            padded[:-2, 1:-1] & padded[2:, 1:-1] & padded[1:-1, :-2] & padded[1:-1, 2:]
        )
        out = out & neigh
    return out


def incidence_mask(theta_deg: np.ndarray, cfg: MaskConfig | None = None) -> np.ndarray:
    """入射角掩膜：Dubois 等模型适用域为 30°–45°，工程上放宽到 25°–45°。"""
    cfg = cfg or MaskConfig()
    theta = np.asarray(theta_deg, dtype=float)
    return (theta >= cfg.theta_min) & (theta <= cfg.theta_max) & np.isfinite(theta)


def freeze_thaw_mask(
    air_temp_c: np.ndarray | None = None,
    soil_temp_c: np.ndarray | None = None,
    cfg: MaskConfig | None = None,
) -> np.ndarray:
    """冻融期掩膜。冻结土壤的介电行为与液态水完全不同，必须剔除。

    只传其中一个亦可；两个都传时取交集。
    """
    cfg = cfg or MaskConfig()
    keep = None
    if air_temp_c is not None:
        t = np.asarray(air_temp_c, dtype=float)
        keep = np.isfinite(t) & (t > cfg.air_temp_threshold_c)
    if soil_temp_c is not None:
        t = np.asarray(soil_temp_c, dtype=float)
        m = np.isfinite(t) & (t > cfg.soil_temp_threshold_c)
        keep = m if keep is None else (keep & m)
    if keep is None:
        raise ValueError("air_temp_c 与 soil_temp_c 至少提供一个。")
    return keep


def precipitation_mask(
    precip_mm: np.ndarray,
    cfg: MaskConfig | None = None,
) -> np.ndarray:
    """降水事件掩膜。截留水与冠层水不等于土壤水，需打标。

    参数
    ----
    precip_mm : (n_dates, ...) 逐日降水，第一维为时间
    """
    cfg = cfg or MaskConfig()
    p = np.asarray(precip_mm, dtype=float)
    keep = np.ones(p.shape, dtype=bool)
    hits = np.where(p > cfg.precip_threshold_mm)
    for t_idx in np.unique(hits[0]):
        sl = slice(int(t_idx), int(t_idx) + cfg.precip_lag_days + 1)
        keep[sl] = False
    return keep


def vwc_domain_mask(vwc: np.ndarray, limit: float | None = None, cfg: MaskConfig | None = None) -> np.ndarray:
    """适用域掩膜：相干性支路只在低植被区承诺精度。"""
    cfg = cfg or MaskConfig()
    lim = cfg.vwc_limit_kg_m2 if limit is None else limit
    v = np.asarray(vwc, dtype=float)
    return np.isfinite(v) & (v < lim)


def combine_masks(masks, require_all: bool = True) -> np.ndarray:
    """合并多个掩膜。require_all=True 取交集（默认，保守）。"""
    masks = list(masks)
    if not masks:
        raise ValueError("至少提供一个掩膜。")
    stacked = np.stack([np.asarray(m, dtype=bool) for m in masks], axis=0)
    return np.all(stacked, axis=0) if require_all else np.any(stacked, axis=0)
