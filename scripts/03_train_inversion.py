#!/usr/bin/env python
"""步骤 3：训练 NN 反演器并用实测站点做残差校正。

用法：
    python scripts/03_train_inversion.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.retrieval.inversion import InversionModel, ResidualCorrector, SyntheticDBConfig, build_synthetic_db
from ecohydro_sm.validation.metrics import full_report


def run_synthetic_check() -> None:
    """端到端：合成库 -> 训练 -> 残差校正 -> 指标。"""
    X, y, _ = build_synthetic_db(SyntheticDBConfig())
    model = InversionModel(estimator="rf").fit(X, y)
    pred = model.predict(X)

    rep = full_report(pred, y)
    print(f"[自检] 合成库自洽性 ubRMSE={rep['ubrmse']:.4f} m3/m3, R={rep['pearson_r']:.3f}")
    assert rep["ubrmse"] < 0.05, "模型在合成库上都未收敛，检查正向模型或特征"

    # 残差校正：模拟 30 个站点、含固定系统偏差的情形
    rng = np.random.default_rng(1)
    idx = rng.choice(len(y), 30, replace=False)
    resid = y[idx] - pred[idx] + 0.02  # 加入 +0.02 的系统偏差
    corr = ResidualCorrector(min_samples_for_model=20).fit(X[idx], resid)
    fixed = corr.apply(X[idx], pred[idx])
    print(f"[自检] 残差校正前 bias={np.mean(pred[idx] - y[idx]):+.4f} -> "
          f"校正后 bias={np.mean(fixed - y[idx]):+.4f} m3/m3")


def main() -> int:
    ap = argparse.ArgumentParser(description="训练土壤水反演器")
    ap.add_argument("--config", required=True)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    p = cfg["l1_retrieval"]
    db_path = Path(cfg["paths"]["output"]) / "synthetic_db.npz"
    if not db_path.exists():
        print(f"[错误] 未找到合成库 {db_path}，请先运行 02_build_synthetic_db.py")
        return 1

    z = np.load(db_path)
    model = InversionModel(estimator=p["estimator"]).fit(z["X"], z["y"])
    print(f"[完成] 反演器训练：{z['X'].shape[0]} 样本，estimator={p['estimator']}")
    print("[待接] 接实测站点后调用 ResidualCorrector.fit 完成二级校正并落盘模型。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
