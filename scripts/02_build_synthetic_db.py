#!/usr/bin/env python
"""步骤 2：用 WCM 正向模型生成合成训练库。

用法：
    python scripts/02_build_synthetic_db.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.retrieval.inversion import SyntheticDBConfig, build_synthetic_db
from ecohydro_sm.retrieval.wcm import linear_to_db, wcm_forward


def run_synthetic_check() -> None:
    """校验正向模型的物理一致性与合成库规模。"""
    mv = np.linspace(0.05, 0.40, 8)
    s_dry = linear_to_db(wcm_forward(mv, 1.0, 0.0, 35.0))
    s_wetveg = linear_to_db(wcm_forward(mv, 1.0, 4.0, 35.0))
    d = s_wetveg - s_dry
    assert np.all(np.diff(s_dry) > 0), "σ0 未随含水量单调上升"
    assert np.all(d > 0), "植被未产生预期的散射贡献"
    print(f"[自检] 裸土 σ0_VV {s_dry[0]:.1f} -> {s_dry[-1]:.1f} dB（mv 0.05→0.40）")
    print(f"[自检] VWC=4 时植被贡献 {np.mean(d):.2f} dB")

    cfg = SyntheticDBConfig()
    X, y, meta = build_synthetic_db(cfg)
    assert X.shape[0] == y.shape[0] and np.all(np.isfinite(X)), "合成库含非法值"
    print(f"[自检] 合成库 {X.shape[0]} 样本 × {X.shape[1]} 特征；"
          f"mv {meta['mv_range'][0]:.2f}~{meta['mv_range'][1]:.2f}")


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 WCM 合成训练库")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    p = cfg["l1_retrieval"]
    out_dir = Path(args.output or cfg["paths"]["output"])
    out_dir.mkdir(parents=True, exist_ok=True)

    db_cfg = SyntheticDBConfig(
        mv_min=p["synthetic_db"]["mv_min"], mv_max=p["synthetic_db"]["mv_max"],
        h_rms_min=p["synthetic_db"]["h_rms_min"], h_rms_max=p["synthetic_db"]["h_rms_max"],
        vwc_max=p["synthetic_db"]["vwc_max"],
        theta_min=p["synthetic_db"]["theta_min"], theta_max=p["synthetic_db"]["theta_max"],
        freq_ghz=p["freq_ghz"], pol=p["pol"],
        a_param=p["wcm"]["a_param"], b_param=p["wcm"]["b_param"], e_param=p["wcm"]["e_param"],
        noise_db=p["synthetic_db"]["noise_db"],
    )
    X, y, meta = build_synthetic_db(db_cfg)
    np.savez_compressed(out_dir / "synthetic_db.npz", X=X, y=y)
    print(f"[完成] 合成库写入 {out_dir / 'synthetic_db.npz'}：{X.shape}")
    print("[警告] 本步仅为先验分布；A/B/E 与粗糙度范围必须用在地实测标定后再用于生产。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
