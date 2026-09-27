"""水量平衡闭合 -> 根区土壤水（本仓库的科学核心）。

闭合方程
--------
    ΔS_vad = P − ET − R − ΔS_gw

得到包气带（含根区）储变量后，以 L1/L2 反演的表层土壤水为上边界，
沿深度分配，得到根区土壤水剖面。

这一步正是 Han et al. (2020)「绿化耗地下水」与 Tian et al. (2026)
「绿化经能量介导提升降水效率」之争所缺的闭合项——
把"地下水被消耗了多少"从推测变成有观测约束的量。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "WaterBalanceConfig",
    "vadose_storage_change",
    "closure_residual",
    "allocate_to_root_zone",
    "closure_quality_flag",
]


@dataclass
class WaterBalanceConfig:
    """水量平衡参数。所有量统一为 mm。"""

    runoff_mm: float = 0.0                # 内流干旱区可近似为 0
    irrigation_return_fraction: float = 0.0  # 灌溉回归水比例（民勤渠系渗漏需考虑）
    root_depth_m: float = 1.0
    vadose_depth_m: float = 3.0           # 包气带参与水分交换的厚度
    max_closure_residual_frac: float = 0.15  # 闭合残差 / 降水 的可接受上限
    decay_shape: float = 1.5              # 垂向分配指数（越大越集中在浅层）


def vadose_storage_change(
    precip_mm: np.ndarray,
    et_mm: np.ndarray,
    delta_s_gw_mm: np.ndarray,
    cfg: WaterBalanceConfig | None = None,
) -> np.ndarray:
    """ΔS_vad = P − ET − R − ΔS_gw。"""
    cfg = cfg or WaterBalanceConfig()
    p = np.asarray(precip_mm, dtype=float)
    e = np.asarray(et_mm, dtype=float)
    g = np.asarray(delta_s_gw_mm, dtype=float)
    return p - e - cfg.runoff_mm - g


def closure_residual(
    precip_mm: np.ndarray,
    et_mm: np.ndarray,
    delta_s_gw_mm: np.ndarray,
    delta_s_vad_mm: np.ndarray,
    cfg: WaterBalanceConfig | None = None,
) -> np.ndarray:
    """闭合残差 = P − ET − R − ΔS_gw − ΔS_vad，理论上应为 0。"""
    cfg = cfg or WaterBalanceConfig()
    return (
        np.asarray(precip_mm, dtype=float)
        - np.asarray(et_mm, dtype=float)
        - cfg.runoff_mm
        - np.asarray(delta_s_gw_mm, dtype=float)
        - np.asarray(delta_s_vad_mm, dtype=float)
    )


def allocate_to_root_zone(
    delta_s_vad_mm: np.ndarray,
    surface_sm: np.ndarray,
    cfg: WaterBalanceConfig | None = None,
) -> np.ndarray:
    """把包气带储变量沿深度分配为根区平均土壤水。

    分配权重按指数递减：w(z) ∝ exp(−z / decay_length)，
    浅层变化幅度大于深层（物理上正确：入渗波向下传播时衰减）。

    参数
    ----
    delta_s_vad_mm : 包气带储变量变化（mm 水当量）
    surface_sm : 表层体积含水量（m3/m3，作为上边界约束）

    返回
    ----
    root_zone_sm : 根区平均体积含水量（m3/m3）
    """
    cfg = cfg or WaterBalanceConfig()
    z = np.linspace(0.0, cfg.vadose_depth_m, 50)
    weights = np.exp(-z / max(cfg.decay_shape * cfg.root_depth_m, 1e-6))
    weights /= weights.sum()

    # 只取根区深度内的权重并按厚度换算为体积含水量增量
    in_root = z <= cfg.root_depth_m
    w_root = weights[in_root] / weights[in_root].sum()
    dz = np.diff(np.concatenate([z[in_root], [cfg.root_depth_m]]))
    # mm 水当量 -> 体积含水量增量：Δθ = ΔS(mm) / 厚度(mm)
    thickness_mm = cfg.root_depth_m * 1000.0
    delta_theta = np.asarray(delta_s_vad_mm, dtype=float) / thickness_mm
    sm_surface = np.asarray(surface_sm, dtype=float)

    # 表层权重最大，向根区底部递减（入渗波向下传播时幅度衰减）
    rel = w_root / w_root[0]
    profile = sm_surface[..., None] + delta_theta[..., None] * rel
    root_mean = np.average(profile, axis=-1, weights=dz)
    return np.clip(root_mean, 0.0, 0.6)


def closure_quality_flag(
    precip_mm: np.ndarray,
    residual_mm: np.ndarray,
    cfg: WaterBalanceConfig | None = None,
) -> np.ndarray:
    """闭合质量标记：|残差| <= 阈值 · 降水 才认为该时段可用。"""
    cfg = cfg or WaterBalanceConfig()
    p = np.asarray(precip_mm, dtype=float)
    r = np.abs(np.asarray(residual_mm, dtype=float))
    limit = cfg.max_closure_residual_frac * np.maximum(p, 1.0)
    return np.isfinite(r) & (r <= limit)
