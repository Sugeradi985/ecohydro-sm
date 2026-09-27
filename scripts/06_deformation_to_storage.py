#!/usr/bin/env python
"""步骤 6：L4 形变 -> 含水层储变量 ΔS。

用法：
    python scripts/06_deformation_to_storage.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.insar.deformation import (
    DeformationConfig,
    estimate_lag_days,
    fit_decomposition,
    seasonal_amplitude_snr,
)
from ecohydro_sm.insar.storage import (
    StorageConfig,
    calibrate_ske,
    storage_change_from_deformation,
    storage_change_from_head,
    sy_sensitivity_span,
)


def run_synthetic_check() -> None:
    """合成自检：分解、滞后校正、储变量区间。"""
    rng = np.random.default_rng(4)
    n_t, n_px = 120, 200
    t = np.arange(n_t) * 6.0  # 6 天重访

    seasonal = -15.0 * np.sin(2 * np.pi * t / 365.25)[:, None] * (1 + rng.normal(0, 0.05, (1, n_px)))
    trend = (-8.0 * (t / 365.25))[:, None] * np.ones((1, n_px))
    disp = seasonal + trend + rng.normal(0, 2.0, (n_t, n_px))

    decomp = fit_decomposition(t, disp)
    snr = seasonal_amplitude_snr(decomp)
    print(f"[自检] 季节振幅 {np.nanmean(decomp['amplitude_mm']):.1f} mm；"
          f"SNR 中位数 {np.nanmedian(snr):.2f}（判据 >= 3）")

    head = -3.0 * np.sin(2 * np.pi * t / 365.25) - 0.5
    lag = estimate_lag_days(seasonal, head, t)
    print(f"[自检] 形变—水位滞后 {lag} 天")

    ske = calibrate_ske(decomp["seasonal"], head)
    print(f"[自检] 标定 S_ke = {ske:.2e}（真值 5.0e-3）")
    assert abs(ske - 5.0e-3) < 2e-3, "S_ke 标定偏差过大"

    ds_a = storage_change_from_deformation(decomp["seasonal"], ske, 0.10)
    ds_b = storage_change_from_head(head, 0.10)
    print(f"[自检] 路径A ΔS 幅度 {np.nanmean(np.abs(ds_a)):.1f} mm；"
          f"路径B ΔS 幅度 {np.nanmean(np.abs(ds_b)):.1f} mm")

    span = sy_sensitivity_span(decomp["seasonal"], ske, StorageConfig())
    width = span["p95_mm"] - span["p05_mm"]
    print(f"[自检] S_y 敏感性：|ΔS| 均值 {np.nanmean(np.abs(span['median_mm'])):.1f} mm，"
          f"区间宽度均值 {np.nanmean(width):.1f} mm")
    print("[提示] 区间跨度远大于像元噪声 —— 报告必须给区间，禁止单点估计。")


def main() -> int:
    ap = argparse.ArgumentParser(description="形变转含水层储变量")
    ap.add_argument("--config", required=True)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    s = cfg["l4_storage"]
    if not s["enabled"]:
        print("[跳过] 配置中 l4_storage.enabled = false")
        return 0

    ts_path = Path(cfg["paths"]["mintpy"]) / "timeseries.h5"
    if not ts_path.exists():
        print(f"[跳过] 未找到时序形变：{ts_path}")
        print("       先用 MintPy smallbaselineApp.py 产出 timeseries.h5。")
        return 0

    _ = DeformationConfig(periods_days=tuple(s["decomposition"]["periods_days"]),
                          trend_order=s["decomposition"]["trend_order"],
                          max_lag_days=s["decomposition"]["max_lag_days"])
    print(f"[就绪] 模式 {s['mode']}；最小季节 SNR {s['min_seasonal_snr']}")
    print("[警告] Δd 不等于 ΔS，二者通过 S_y/S_ke 联系，跨数量级，必须做敏感性分析。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
