#!/usr/bin/env python
"""步骤 7：L3 水量平衡闭合 -> 根区土壤水。

    ΔS_vad = P − ET − R − ΔS_gw

用法：
    python scripts/07_water_balance_closure.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.fusion.water_balance import (
    WaterBalanceConfig,
    allocate_to_root_zone,
    closure_quality_flag,
    closure_residual,
    vadose_storage_change,
)


def run_synthetic_check() -> None:
    """合成自检：闭合残差应为 0，根区土壤水应在合理范围。"""
    rng = np.random.default_rng(5)
    n_t, n_px = 60, 300

    p = rng.gamma(1.2, 8.0, (n_t, n_px))          # 干旱区降水：多数为小事件
    et = np.clip(p * 0.7 + rng.normal(5, 2, (n_t, n_px)), 0, None)
    d_gw = rng.normal(-2.0, 3.0, (n_t, n_px))     # 地下水净亏损为负
    surface_sm = np.clip(rng.normal(0.12, 0.04, (n_t, n_px)), 0.02, 0.35)

    d_vad = vadose_storage_change(p, et, d_gw)
    resid = closure_residual(p, et, d_gw, d_vad)
    assert np.nanmax(np.abs(resid)) < 1e-9, "闭合方程未成立，检查符号约定"
    print(f"[自检] 闭合残差 max={np.nanmax(np.abs(resid)):.2e} mm（应为 0）")

    flag = closure_quality_flag(p, resid)
    print(f"[自检] 闭合合格时段占比 {flag.mean() * 100:.0f}%")

    sm_rz = allocate_to_root_zone(d_vad, surface_sm)
    assert np.all((sm_rz >= 0) & (sm_rz <= 0.6)), "根区土壤水越界"
    print(f"[自检] 根区 SM 均值 {np.nanmean(sm_rz):.3f} m3/m3；"
          f"表层均值 {np.nanmean(surface_sm):.3f} m3/m3")


def main() -> int:
    ap = argparse.ArgumentParser(description="水量平衡闭合与根区土壤水")
    ap.add_argument("--config", required=True)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    w = cfg["l3_water_balance"]
    if not w["enabled"]:
        print("[跳过] 配置中 l3_water_balance.enabled = false")
        return 0

    out = Path(cfg["paths"]["output"])
    need = [out / "delta_s_gw.npy", out / "l1_surface_sm.npy"]
    missing = [str(x) for x in need if not x.exists()]
    if missing:
        print("[跳过] 缺少上游产品：" + ", ".join(missing))
        print("       请先运行 04（表层 SM）与 06（储变量）。")
        return 0

    _ = WaterBalanceConfig(**w)
    print(f"[就绪] 根区深度 {w['root_depth_m']} m；允许闭合残差 ≤ "
          f"{w['max_closure_residual_frac'] * 100:.0f}% 降水")
    print("[提醒] 民勤为内流区，但渠系渗漏与灌溉回归水需单独计入。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
