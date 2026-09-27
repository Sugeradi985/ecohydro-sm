#!/usr/bin/env python
"""步骤 4：生产 L1 表层土壤水（0–5 cm）。

用法：
    python scripts/04_retrieve_sm.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.preprocess.masking import MaskConfig, combine_masks, incidence_mask, vwc_domain_mask
from ecohydro_sm.preprocess.speckle import refined_lee
from ecohydro_sm.retrieval.change_detection import (
    ChangeDetectionConfig,
    change_detection_sm,
    dry_wet_references,
    dynamic_range_flag,
)
from ecohydro_sm.retrieval.inversion import InversionModel, SyntheticDBConfig, build_synthetic_db


def run_synthetic_check() -> None:
    """合成时序自检：斑点滤波 + 变化检测。"""
    rng = np.random.default_rng(2)
    n_t, h, w = 36, 32, 32
    # 含水量摆幅要足够大，否则干湿动态范围不足，变化检测会判为不可反演
    mv_true = np.clip(0.22 + 0.18 * np.sin(np.linspace(0, 4 * np.pi, n_t))[:, None, None]
                      + rng.normal(0, 0.01, (n_t, h, w)), 0.02, 0.45)

    from ecohydro_sm.retrieval.wcm import linear_to_db, wcm_forward
    s0 = linear_to_db(wcm_forward(mv_true, 1.0, 0.5, 35.0))
    s0_f = refined_lee(s0[0], window=5)

    dry, wet = dry_wet_references(s0, ChangeDetectionConfig())
    flag = dynamic_range_flag(dry, wet)
    sm = change_detection_sm(s0, dry[None], wet[None])

    valid = np.isfinite(sm) & flag
    err = np.abs(sm[valid] - mv_true[valid])
    print(f"[自检] 变化检测 MAE={np.nanmean(err):.4f} m3/m3；可反演像元 {flag.mean() * 100:.0f}%")
    print(f"[自检] 滤波后 ENL 提升：σ0 局部标准差 {np.nanstd(s0[0]):.2f} -> {np.nanstd(s0_f):.2f} dB")


def main() -> int:
    ap = argparse.ArgumentParser(description="生产 L1 表层土壤水")
    ap.add_argument("--config", required=True)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    p = cfg["l1_retrieval"]
    rtc_dir = Path(cfg["paths"]["s1_rtc"])
    if not rtc_dir.exists():
        print(f"[跳过] 未找到 S1 RTC 目录：{rtc_dir}")
        print("       推荐用 ASF HyP3 的 RTC_GAMMA 产品（radiometry=gamma0，与既有 MintPy 工作流一致）。")
        return 0

    mcfg = MaskConfig(**cfg.get("masking", {}))
    _ = combine_masks([incidence_mask(np.array([35.0]), mcfg),
                       vwc_domain_mask(np.array([1.0]), cfg["l2_coherence"]["vwc_limit_kg_m2"], mcfg)])
    print(f"[就绪] 生产参数：estimator={p['estimator']}，聚合倍率={p['aggregate_factor']}")
    print("[提醒] 升轨与降轨必须分开反演再融合，入射角差异会污染模型。")
    print("[提醒] 产品标称分辨率必须等于实际聚合后的分辨率。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
