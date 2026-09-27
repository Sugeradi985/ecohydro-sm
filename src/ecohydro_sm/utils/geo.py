"""地理工具：投影、网格重采样、点提取。

只依赖 numpy，避免给核心算法引入重依赖。
涉及栅格读写的功能放在 `ecohydro_sm.io.raster`（可选依赖 rasterio）。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "project_to_meters",
    "grid_from_extent",
    "nearest_grid_index",
    "extract_point_series",
    "block_aggregate",
]


def project_to_meters(lon: np.ndarray, lat: np.ndarray, epsg: str | None = None):
    """经纬度 -> 平面米制坐标。

    需要 pyproj（可选依赖）。**空间分块前必须投影**：
    直接按经纬度分块会让块面积随纬度变化，破坏空间 CV 的均匀性。
    """
    try:
        from pyproj import Transformer
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "投影需要 pyproj：pip install pyproj。"
            "若研究区跨度小，也可直接用局部等距近似代替。"
        ) from exc

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    if epsg is None:
        # 按中心经度自动选 UTM 带
        zone = int(np.floor((np.nanmean(lon) + 180) / 6) + 1)
        hemisphere = "north" if np.nanmean(lat) >= 0 else "south"
        epsg = f"EPSG:{32600 + zone if hemisphere == 'north' else 32700 + zone}"
    tf = Transformer.from_crs("EPSG:4326", epsg, always_xy=True)
    x, y = tf.transform(lon, lat)
    return np.asarray(x), np.asarray(y), epsg


def grid_from_extent(x_min: float, x_max: float, y_min: float, y_max: float, res: float):
    """由范围与分辨率生成规则网格中心点坐标。"""
    xs = np.arange(x_min + res / 2, x_max, res)
    ys = np.arange(y_max - res / 2, y_min, -res)
    return xs, ys


def nearest_grid_index(x: np.ndarray, y: np.ndarray, xs: np.ndarray, ys: np.ndarray):
    """点 -> 最近网格索引 (row, col)。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    col = np.argmin(np.abs(xs[None, :] - x[:, None]), axis=1)
    row = np.argmin(np.abs(ys[None, :] - y[:, None]), axis=1)
    return row, col


def extract_point_series(stack: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """从 (n_time, H, W) 数据立方体中提取若干点的时序。"""
    stack = np.asarray(stack)
    return stack[:, rows, cols]


def block_aggregate(arr: np.ndarray, factor: int, method: str = "mean") -> np.ndarray:
    """按整数倍降尺度聚合（用于把 10 m 产品聚合到 100 m / 1 km）。

    降尺度前先聚合能显著降低斑点噪声，但会牺牲空间细节；
    **产品标称分辨率必须等于实际聚合后的分辨率**，
    否则就是 Brocca et al. (2024) 批评的"名义高分辨率"。
    """
    arr = np.asarray(arr, dtype=float)
    if arr.ndim < 2:
        raise ValueError("arr 至少二维。")
    h, w = arr.shape[-2] // factor * factor, arr.shape[-1] // factor * factor
    cropped = arr[..., :h, :w]
    new_shape = cropped.shape[:-2] + (h // factor, factor, w // factor, factor)
    blocked = cropped.reshape(new_shape)
    if method == "mean":
        return np.nanmean(blocked, axis=(-3, -1))
    if method == "median":
        return np.nanmedian(blocked, axis=(-3, -1))
    raise ValueError(f"method 必须是 'mean' 或 'median'，收到 {method!r}")
