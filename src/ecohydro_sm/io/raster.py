"""栅格读写（可选依赖 rasterio）。

设计原则：**核心算法模块不 import 本模块**，保证纯数值测试与合成数据试验
不需要装 rasterio。
"""

from __future__ import annotations

import numpy as np

__all__ = ["read_stack", "write_geotiff", "stack_dates_from_names"]


def read_stack(paths, band: int = 1, masked: bool = True) -> tuple[np.ndarray, dict]:
    """按顺序读取多个 GeoTIFF 为一个 (n_time, H, W) 数据立方体。

    返回 (stack, profile)。profile 取自第一个文件。
    """
    import rasterio  # 局部导入，保持可选依赖

    arrays = []
    profile = None
    for p in paths:
        with rasterio.open(p) as src:
            arr = src.read(band, masked=masked)
            arrays.append(np.ma.filled(arr, np.nan) if masked else arr)
            if profile is None:
                profile = src.profile.copy()
    stack = np.stack(arrays, axis=0).astype(float)
    if profile is not None:
        profile.update(count=len(paths))
    return stack, (profile or {})


def write_geotiff(path, array: np.ndarray, profile: dict) -> None:
    """写单波段或多波段 GeoTIFF。array 形状为 (H, W) 或 (n, H, W)。"""
    import rasterio

    arr = np.asarray(array)
    if arr.ndim == 2:
        arr = arr[None, ...]
    prof = dict(profile)
    prof.update(count=arr.shape[0], dtype="float32")
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(arr.astype("float32"))


def stack_dates_from_names(paths, pattern: str = "%Y%m%d") -> np.ndarray:
    """从文件名解析日期。默认匹配 Sentinel 命名中的 8 位日期。"""
    from datetime import datetime
    import re

    dates = []
    for p in paths:
        m = re.search(r"\d{8}", str(p))
        if not m:
            raise ValueError(f"无法从 {p} 解析日期，请检查命名或自定义 pattern。")
        dates.append(datetime.strptime(m.group(0), pattern))
    return np.array(dates, dtype="datetime64[ns]")
