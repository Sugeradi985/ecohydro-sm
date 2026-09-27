#!/usr/bin/env python
"""步骤 5：L2 相干性支路 —— 相干性 + 光学结构变化去混反演表层土壤水。

适用域：VWC < 2 kg/m2 的裸土、稀疏植被与干旱半干旱区。
农田绿洲内部去相关严重，只打标不承诺精度。

用法：
    python scripts/05_insar_coherence_sm.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.insar.coherence import (
    CoherenceSMConfig,
    calibrate_decay_coefficients,
    closure_phase_residual,
    invert_sm_from_coherence,
    structural_change_index,
)
from ecohydro_sm.validation.metrics import full_report


def run_synthetic_check() -> None:
    """合成自检：验证"结构去混"确实能降低反演误差。"""
    rng = np.random.default_rng(3)
    n_t, n_px = 60, 400

    true_dsm = np.clip(rng.normal(0.12, 0.06, (n_t, n_px)), 0.0, 0.35)
    fvc_t = np.clip(np.outer(np.sin(np.linspace(0, np.pi, n_t)), np.ones(n_px)) * 0.6
                    + rng.normal(0, 0.02, (n_t, n_px)), 0, 1)
    dv = structural_change_index(fvc_t)

    k_sm_true, k_v_true = 2.5, 1.8
    gamma = np.exp(-k_sm_true * true_dsm - k_v_true * dv) * (0.9 + rng.normal(0, 0.03, (n_t, n_px)))
    gamma = np.clip(gamma, 0.05, 1.0)

    cal = calibrate_decay_coefficients(gamma, dv, true_dsm)
    k_sm, k_v, c0 = cal["k_sm"], cal["k_v"], cal["log_offset"]
    sm_hat = invert_sm_from_coherence(gamma, dv, k_sm, k_v, log_offset=c0)

    # 对照：忽略结构项（k_v = 0），即 INSAR4SM 式的处理
    sm_naive = invert_sm_from_coherence(gamma, np.zeros_like(dv), k_sm, 0.0, log_offset=c0)
    truth = 0.02 + true_dsm

    r1 = full_report(sm_hat, truth)
    r2 = full_report(sm_naive, truth)
    print(f"[自检] 标定 k_sm={k_sm:.2f}（真值 {k_sm_true}）、k_v={k_v:.2f}（真值 {k_v_true}）")
    print(f"[自检] 含结构去混   RMSE={r1['rmse']:.4f}、bias={r1['bias']:+.4f} m3/m3")
    print(f"[自检] 不考虑结构项 RMSE={r2['rmse']:.4f}、bias={r2['bias']:+.4f} m3/m3")
    # 关键：结构项带来的主要是**系统偏差**，ubRMSE 会把偏差剔除、从而掩盖改善。
    # 因此比较必须用总 RMSE（或同时看 bias），不能用 ubRMSE。
    assert r1["rmse"] <= r2["rmse"] + 1e-6, "结构去混未带来改善，检查模型符号"
    assert abs(r1["bias"]) < abs(r2["bias"]), "结构去混应显著降低系统偏差"

    trip = rng.normal(0, 0.2, (30, 3, n_px))
    resid = closure_phase_residual(trip)
    print(f"[自检] 闭合相位残差 std 中位数 {np.nanmedian(resid):.3f} rad（用于独立质控）")


def main() -> int:
    ap = argparse.ArgumentParser(description="InSAR 相干性土壤水反演")
    ap.add_argument("--config", required=True)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    c = cfg["l2_coherence"]
    if not c["enabled"]:
        print("[跳过] 配置中 l2_coherence.enabled = false")
        return 0

    mintpy_dir = Path(cfg["paths"]["mintpy"])
    if not mintpy_dir.exists():
        print(f"[跳过] 未找到 MintPy 工作目录：{mintpy_dir}")
        print("       流程：HyP3 干涉 -> prep_hyp3 -> smallbaselineApp.py -> coherence.h5")
        return 0

    _ = CoherenceSMConfig()
    print(f"[就绪] VWC 适用上限 {c['vwc_limit_kg_m2']} kg/m2；"
          f"闭合相位残差阈值 {c['closure_phase_max_std_rad']} rad")
    print("[提醒] 干涉对只用短时空基线（<= 12 天 / <= 150 m），几何去相干才可控。")
    print("[提醒] 相干性支路不需要大气校正；形变支路必须做。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
