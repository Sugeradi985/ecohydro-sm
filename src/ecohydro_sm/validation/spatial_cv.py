"""空间块交叉验证。

为什么必须用它
--------------
大量已发表研究用随机划分做交叉验证，相邻像元/同一站点的样本同时出现在
训练集与测试集，造成空间自相关泄漏，R² 被系统性抬高。
改用空间块划分后，欧洲 10 m、113 个 ISMN 站点的实验里 R² 从常见的
0.7–0.8 掉到约 0.51。

**不做空间块 CV，投稿基本会被方法学审稿人拦下。**
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.model_selection import BaseCrossValidator

__all__ = [
    "SpatialBlockCVConfig",
    "assign_spatial_blocks",
    "SpatialBlockKFold",
    "leave_one_station_out",
]


@dataclass
class SpatialBlockCVConfig:
    """空间分块参数。"""

    block_size_m: float = 10_000.0   # 方块边长（米），建议 >= 10 km
    min_samples_per_fold: int = 10


def assign_spatial_blocks(x: np.ndarray, y: np.ndarray, block_size_m: float = 10_000.0) -> np.ndarray:
    """把点按规则网格分配块 ID。x, y 为**投影坐标**（米）。

    注意：必须先投影到等距/等积投影，用经纬度直接分块会导致块面积随纬度变化。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    bx = np.floor(x / block_size_m).astype(np.int64)
    by = np.floor(y / block_size_m).astype(np.int64)
    # 用字符串组合避免大坐标下的整数碰撞
    _, ids = np.unique(np.char.add(bx.astype(str), "_" + by.astype("U")), return_inverse=True)
    return ids.astype(np.int64)


class SpatialBlockKFold(BaseCrossValidator):
    """按空间块划分的 K 折交叉验证，可直接喂给 sklearn 的 cross_val_score。"""

    def __init__(self, blocks: np.ndarray, n_splits: int = 5, shuffle: bool = True, random_state: int = 42):
        self.blocks = np.asarray(blocks)
        self.n_splits = int(n_splits)
        self.shuffle = shuffle
        self.random_state = random_state

    def _unique_blocks(self):
        ub = np.unique(self.blocks)
        if self.shuffle:
            rng = np.random.default_rng(self.random_state)
            ub = rng.permutation(ub)
        return ub

    def get_n_splits(self, X=None, y=None, groups=None) -> int:
        return self.n_splits

    def split(self, X, y=None, groups=None):
        ub = self._unique_blocks()
        folds = np.array_split(ub, self.n_splits)
        idx = np.arange(len(self.blocks))
        for fold in folds:
            test_mask = np.isin(self.blocks, fold)
            if test_mask.sum() == 0:
                continue
            yield idx[~test_mask], idx[test_mask]


def leave_one_station_out(stations: np.ndarray):
    """留一站出：每个站点轮流做验证集。

    参数
    ----
    stations : 每个样本所属站点 ID
    """
    stations = np.asarray(stations)
    idx = np.arange(len(stations))
    for s in np.unique(stations):
        test = stations == s
        yield idx[~test], idx[test]
