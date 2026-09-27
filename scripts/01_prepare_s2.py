#!/usr/bin/env python
"""步骤 1：由 Sentinel-2 L2A 生成植被指数与质量控制掩膜。

用法：
    python scripts/01_prepare_s2.py --config configs/minqin.yaml
    python scripts/01_prepare_s2.py --config configs/minqin.yaml --selftest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config
from ecohydro_sm.preprocess.indices import fvc, ndmi, ndvi, vwc_from_ndvi
from ecohydro_sm.preprocess.masking import MaskConfig, combine_masks, s2_cloud_mask


def run_synthetic_check() -> None:
    """无数据自检：验证指数与掩膜逻辑。"""
    rng = np.random.default_rng(0)
    red = rng.uniform(0.02, 0.15, (64, 64))
    nir = rng.uniform(0.15, 0.55, (64, 64))
    swir = rng.uniform(0.10, 0.40, (64, 64))
    # 真实云是空间集聚的，不是随机噪声；随机撒点会被腐蚀算法整片吃掉
    scl = np.full((64, 64), 4, dtype=int)
    scl[10:22, 15:30] = 9    # 一块云
    scl[24:27, 17:31] = 3    # 云的阴影
    scl[45:50, 40:50] = 8

    v = ndvi(red, nir)
    m = ndmi(nir, swir)
    f = fvc(v, ndvi_soil=np.nanpercentile(v, 5), ndvi_veg=np.nanpercentile(v, 95))
    vwc = vwc_from_ndvi(v)
    keep_cloud = s2_cloud_mask(scl, buffer_pixels=1)
    keep = combine_masks([keep_cloud, np.isfinite(v)])

    assert np.all(np.abs(v) <= 1.0), "NDVI 越界"
    assert np.all((f >= 0) & (f <= 1)), "FVC 越界"
    assert np.all(vwc >= 0), "VWC 出现负值"
    print(f"[自检] NDVI {np.nanmin(v):.3f}~{np.nanmax(v):.3f} | "
          f"FVC {np.nanmean(f):.3f} | VWC {np.nanmean(vwc):.3f} kg/m2 | "
          f"有效像元 {keep.mean() * 100:.1f}%")


def main() -> int:
    ap = argparse.ArgumentParser(description="Sentinel-2 植被指数与掩膜生成")
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None, help="S2 L2A 目录（缺省用配置中的 paths.s2_l2a）")
    ap.add_argument("--output", default=None)
    ap.add_argument("--selftest", action="store_true", help="用合成数据自检，无需真实影像")
    args = ap.parse_args()

    if args.selftest:
        run_synthetic_check()
        return 0

    cfg = load_config(args.config)
    in_dir = Path(args.input or cfg["paths"]["s2_l2a"])
    out_dir = Path(args.output or cfg["paths"]["output"]) / "s2_indices"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_dir.exists():
        print(f"[跳过] 输入目录不存在：{in_dir}")
        print("       请先用 GEE 或 Copernicus Data Space 下载 S2 L2A 并放到该路径。")
        return 0

    _ = MaskConfig(**cfg.get("masking", {}))
    print(f"[就绪] S2 输入：{in_dir}")
    print(f"[就绪] 输出目录：{out_dir}")
    print(f"[参数] 云边缘缓冲 {cfg['masking']['scl_cloud_buffer_pixels']} 像元；"
          f"降水打标阈值 {cfg['masking']['precip_threshold_mm']} mm")
    print("[提示] ndvi_soil / ndvi_veg 端元必须按研究区 NDVI 直方图百分位重取，勿用默认值。")
    print("[待接] 栅格 I/O 接 ecohydro_sm.io.raster.read_stack 后按场景循环写盘。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
