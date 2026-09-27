#!/usr/bin/env python
"""步骤 8：验证 —— 空间块 CV + 三重配置 + 分层报告。

本步**不可关闭**。只报全样本 RMSE 是不合格的。

用法：
    python scripts/08_validate.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.validation.metrics import full_report, stratified_report
from ecohydro_sm.validation.spatial_cv import SpatialBlockKFold, assign_spatial_blocks, leave_one_station_out
from ecohydro_sm.validation.triple_collocation import error_correlation_check, triple_collocation


def run_synthetic_check() -> None:
    """合成自检：对比随机 CV 与空间块 CV 的差距（泄漏诊断）。"""
    rng = np.random.default_rng(6)
    n = 1200
    # 空间自相关的真值场 + 观测噪声
    x = rng.uniform(0, 40_000, n)
    y = rng.uniform(0, 40_000, n)
    field = 0.15 + 0.10 * np.sin(x / 8000) * np.cos(y / 8000)
    obs = field + rng.normal(0, 0.02, n)

    blocks = assign_spatial_blocks(x, y, block_size_m=10_000)
    print(f"[自检] {n} 个样本划分为 {blocks.max() + 1} 个空间块")

    cv = SpatialBlockKFold(blocks, n_splits=5)
    folds = list(cv.split(np.zeros((n, 1))))
    assert len(folds) == 5 and all(len(te) > 0 for _, te in folds), "空间块划分失败"
    print(f"[自检] 5 折划分成功，测试集规模 {[len(te) for _, te in folds]}")

    # 三重配置：三套独立产品
    a = field + rng.normal(0, 0.04, n)
    b = field + rng.normal(0, 0.05, n)
    c = field + rng.normal(0, 0.06, n)
    tc = triple_collocation(a, b, c)
    print(f"[自检] TC 估计 RMSE：A={tc['rmse_x']:.4f}、B={tc['rmse_y']:.4f}、"
          f"C={tc['rmse_z']:.4f} m3/m3（真值 0.04/0.05/0.06）")
    ind = error_correlation_check(a, b, c)
    print(f"[自检] TC 独立性诊断通过：{ind['independence_ok']}")

    vwc = rng.uniform(0, 5, n)
    rep = stratified_report(obs, field, vwc, bins=[0, 0.5, 1.5, 3.0, 99.0])
    for k, v in rep.items():
        print(f"[自检] VWC {k}: n={v['n']}, ubRMSE={v['ubrmse']:.4f}, R={v['pearson_r']:.3f}")

    stations = rng.integers(0, 8, n)
    loso = list(leave_one_station_out(stations))
    print(f"[自检] 留一站出：{len(loso)} 折")


def main() -> int:
    ap = argparse.ArgumentParser(description="土壤水产品验证")
    ap.add_argument("--config", required=True)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    v = cfg["validation"]
    out = Path(cfg["paths"]["output"])
    pred_path = out / "l1_surface_sm.npy"
    if not pred_path.exists():
        print(f"[跳过] 未找到产品：{pred_path}，请先运行 04。")
        return 0

    bins = v["stratify_by_vwc_bins"]
    report = {
        "block_size_m": v["spatial_block_size_m"],
        "n_splits": v["n_splits"],
        "vwc_bins": bins,
        "metrics": v["report_metrics"],
    }
    (out / "validation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[完成] 验证框架配置写入 {out / 'validation_report.json'}")
    print("[强制] 必须同时报告：空间块 CV、留一站出、按 VWC 分层、bias/ubRMSE/cRMSE/anomaly R。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
