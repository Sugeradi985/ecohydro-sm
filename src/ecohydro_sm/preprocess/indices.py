"""Sentinel-2 光谱指数与植被含水量估算。

所有函数接受 Sentinel-2 L2A 地表反射率（0–1 或 0–10000 均可，
只要红/近红外/短波红外波段量纲一致）。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "ndvi",
    "ndmi",
    "fvc",
    "vwc_from_ndvi",
    "vwc_from_ndvi_quadratic",
]


def _safe_ratio(a, b, eps: float = 1e-6) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = b + a
    return (b - a) / np.where(np.abs(denom) < eps, np.nan, denom)


def ndvi(red, nir) -> np.ndarray:
    """NDVI = (NIR − Red) / (NIR + Red)。用 B4 / B8（或 B8A）。"""
    return _safe_ratio(red, nir)


def ndmi(nir, swir) -> np.ndarray:
    """NDMI = (NIR − SWIR) / (NIR + SWIR)，植被水分胁迫指示。用 B8A / B11。"""
    return _safe_ratio(swir, nir)


def fvc(ndvi_arr, ndvi_soil: float = 0.05, ndvi_veg: float = 0.85) -> np.ndarray:
    """植被覆盖度 FVC，由 NDVI 线性拉伸得到。

    参数
    ----
    ndvi_soil, ndvi_veg : 裸土与全植被端元，**按研究区 NDVI 直方图百分位
        （如 5% 与 95%）重新取值**，不要用默认值硬套。
    """
    x = np.asarray(ndvi_arr, dtype=float)
    f = (x - ndvi_soil) / (ndvi_veg - ndvi_soil)
    return np.clip(f, 0.0, 1.0)


def vwc_from_ndvi(ndvi_arr, slope: float = 5.0, intercept: float = -0.5) -> np.ndarray:
    """NDVI -> 植被含水量 VWC (kg/m2) 的线性经验关系。

    **必须本地标定**：slope / intercept 直接用文献值会造成系统性偏差。
    无实测生物量时至少要做敏感性分析，并在论文中声明。
    """
    x = np.asarray(ndvi_arr, dtype=float)
    return np.clip(slope * x + intercept, 0.0, None)


def vwc_from_ndvi_quadratic(ndvi_arr, a: float = 1.9, b: float = -0.3, c: float = 0.1) -> np.ndarray:
    """NDVI -> VWC 的二次形式，用于冠层密集时饱和的情形。同样需本地标定。"""
    x = np.asarray(ndvi_arr, dtype=float)
    return np.clip(a * x**2 + b * x + c, 0.0, None)
