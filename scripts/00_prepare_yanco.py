"""Yanco 试验区：数据清点、口径核对、时相配对（接入反演前的第 0 步）。

解决的问题
----------
1. HyP3 INSAR_GAMMA 是**逐干涉对**产品，GEE GRD / S-2 是**逐日期**产品，
   两者必须先建立一张明确的配对表，否则 L2 会把 pair 误当成时序。
2. 三类数据的 CRS / 分辨率 / 边界常常不一致（HyP3 是 40/80 m UTM，
   GEE 导出常见 10 m），必须先统一网格。
3. S-1 文件名是 UTC，OzNet 站点时间戳多为澳洲东部当地时，
   降轨清晨过境会让两者整体错开一天。

本脚本只做**只读清点**，不修改任何原始数据；产出 manifest CSV 与一份报告。

用法
----
    python scripts/00_prepare_yanco.py --root . --config configs/yanco.yaml
    python scripts/00_prepare_yanco.py --selftest

不依赖 rasterio 也能跑（此时只做文件名级清点）；装了 rasterio 会额外核对
CRS / 分辨率 / 边界并抽样检查数值范围（判断 GRD 是 dB 还是线性功率）。
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.config import load_config  # noqa: E402

DATE_RE = re.compile(r"(?<!\d)(\d{8})(?!\d)")
# 月度合成产品常按 YYYYMM / YYYY_MM 命名；只认 8 位日期会把它们的日期全丢掉
MONTH_RE = re.compile(r"(?<!\d)(\d{6})(?!\d)")
YEAR_MONTH_RE = re.compile(r"(?<!\d)(\d{4})[-_](\d{2})(?!\d)")

GLOBS = {
    "s1_grd": ["*VV*.tif", "*grd*.tif", "S1*.tif"],
    "s2_l2a": ["*S2*.tif", "*ndvi*.tif", "*s2*.tif"],
    "insar_pairs": ["*corr*.tif", "*_corr.tif", "S1??_*_corr*.tif"],
}


@dataclass
class Check:
    level: str          # PASS | WARN | FAIL | NOTE
    item: str
    message: str


@dataclass
class Inventory:
    root: Path
    grd: list[dict] = field(default_factory=list)
    s2: list[dict] = field(default_factory=list)
    pairs: list[dict] = field(default_factory=list)
    checks: list[Check] = field(default_factory=list)
    match_rows: list[dict] = field(default_factory=list)


def parse_dates(name: str, limit: int = 1) -> list[datetime]:
    """从文件名解析 8 位日期；HyP3 对产品名含两个时间戳（参考 / 副影像）。"""
    out = []
    for token in DATE_RE.findall(name)[:limit]:
        try:
            out.append(datetime.strptime(token, "%Y%m%d"))
        except ValueError:
            continue
    return out


def parse_scene_date(name: str) -> tuple[datetime | None, str]:
    """解析 GRD / S-2 产品日期，返回 (日期, 时间粒度)。

    ``scene`` = 逐景（YYYYMMDD）；``monthly`` = 月度合成（YYYYMM / YYYY_MM）。
    月度产品的日期统一记为该月 1 日，**只是占位**；真实代表日期在
    ``00b_align_temporal.py`` 里按配置取月中，避免把月首当观测日。
    文件名里的 6 位时间戳（如 ``T200000``）月份为 00，会被校验排除。
    """
    m = DATE_RE.search(name)
    if m:
        try:
            return datetime.strptime(m.group(1), "%Y%m%d"), "scene"
        except ValueError:
            pass
    for y, mm in YEAR_MONTH_RE.findall(name):
        if 1990 <= int(y) <= 2100 and 1 <= int(mm) <= 12:
            return datetime(int(y), int(mm), 1), "monthly"
    for tok in MONTH_RE.findall(name):
        y, mm = int(tok[:4]), int(tok[4:])
        if 1990 <= y <= 2100 and 1 <= mm <= 12:
            return datetime(y, mm, 1), "monthly"
    return None, ""


def scan(root: Path, patterns: list[str]) -> list[Path]:
    hits: set[Path] = set()
    for pat in patterns:
        hits.update(root.rglob(pat))
    return sorted(p for p in hits if p.is_file())


def inspect_raster(path: Path) -> dict:
    """可选：读取 CRS/分辨率/范围/数值范围。无 rasterio 时返回空 dict。"""
    try:
        import rasterio
        import numpy as np
    except ImportError:
        return {}
    try:
        with rasterio.open(path) as src:
            arr = src.read(1, masked=True)
            valid = np.asarray(arr.compressed(), dtype=float)
            sample = valid[:: max(1, valid.size // 5000)]
            return {
                "crs": str(src.crs) if src.crs else "NONE",
                "res_m": round(float(src.res[0]), 3),
                "bounds": ",".join(f"{v:.2f}" for v in src.bounds),
                "dtype": str(src.dtypes[0]),
                "vmin": round(float(sample.min()), 4) if sample.size else None,
                "vmax": round(float(sample.max()), 4) if sample.size else None,
                "vmean": round(float(sample.mean()), 4) if sample.size else None,
            }
    except Exception as exc:  # 损坏文件也要让清点继续
        return {"error": f"{type(exc).__name__}: {exc}"}


def build_manifests(inv: Inventory, cfg: dict) -> None:
    paths = cfg.get("paths", {})
    for key in ("s1_grd", "s2_l2a"):
        base = inv.root / paths.get(key, f"data_raw/yanco/{key}")
        for p in scan(base, GLOBS[key]):
            d, gran = parse_scene_date(p.name)
            row = {"path": str(p),
                   "date": d.strftime("%Y-%m-%d") if d else "",
                   "temporal_granularity": gran,
                   "year_month": d.strftime("%Y-%m") if d else "",
                   "file": p.name}
            row.update(inspect_raster(p))
            (inv.grd if key == "s1_grd" else inv.s2).append(row)

    base = inv.root / paths.get("insar_pairs", "data_raw/yanco/hyp3_insar")
    for p in scan(base, GLOBS["insar_pairs"]):
        d = parse_dates(p.name, limit=2)
        row = {"path": str(p), "file": p.name}
        if len(d) >= 2:
            row["ref_date"] = d[0].strftime("%Y-%m-%d")
            row["sec_date"] = d[1].strftime("%Y-%m-%d")
            row["dt_days"] = (d[1] - d[0]).days
        else:
            row["ref_date"] = row["sec_date"] = ""
            row["dt_days"] = ""
        row.update(inspect_raster(p))
        inv.pairs.append(row)


def build_match_table(inv: Inventory, max_gap_days: int = 10) -> None:
    s2_dates = sorted({datetime.strptime(r["date"], "%Y-%m-%d")
                       for r in inv.s2 if r.get("date")})
    pair_dates = sorted({(datetime.strptime(r["ref_date"], "%Y-%m-%d"),
                          datetime.strptime(r["sec_date"], "%Y-%m-%d"))
                         for r in inv.pairs if r.get("ref_date") and r.get("sec_date")})
    for g in sorted(inv.grd, key=lambda x: x.get("date", "")):
        if not g.get("date"):
            continue
        d = datetime.strptime(g["date"], "%Y-%m-%d")
        nearest_s2 = min(s2_dates, key=lambda x: abs((x - d).days)) if s2_dates else None
        gap = abs((nearest_s2 - d).days) if nearest_s2 else None
        linked = [(a, b) for a, b in pair_dates if a <= d <= b]
        inv.match_rows.append({
            "s1_date": g["date"],
            "s2_date": nearest_s2.strftime("%Y-%m-%d") if nearest_s2 else "",
            "gap_days": gap if gap is not None else "",
            "s2_usable": "YES" if gap is not None and gap <= max_gap_days else "NO",
            "n_insar_pairs": len(linked),
            "pair_spans": ";".join(f"{a:%Y%m%d}-{b:%Y%m%d}" for a, b in linked[:3]),
        })


def run_checks(inv: Inventory, cfg: dict) -> None:
    c = inv.checks.append
    l2 = cfg.get("l2_coherence", {})
    inp = cfg.get("input_spec", {})
    timing = cfg.get("timing", {})

    c(Check("NOTE", "时区",
            f"S-1 文件名是 {timing.get('sar_time_basis', 'UTC')}；站点为 "
            f"{timing.get('site_timezone', 'Australia/Sydney')}。降轨约当地 "
            f"{timing.get('descending_local_hour', 6.5)} 点 = 前一日 19-21 UTC，"
            "按文件名日期直接匹配站点当地日期会整体错开一天。"))

    n_grd = len(inv.grd)
    c(Check("PASS" if n_grd >= 30 else "WARN", "GRD 景数",
            f"{n_grd} 景；建议 ≥ 30（约 2 年 × 12 天重访的 60% 以上）。"))

    dated = [r for r in inv.grd if r.get("date")]
    if len(dated) >= 2:
        span = (datetime.strptime(dated[-1]["date"], "%Y-%m-%d")
                - datetime.strptime(dated[0]["date"], "%Y-%m-%d")).days / 365.25
        c(Check("PASS" if span >= 2 else "WARN", "时间跨度",
                f"{span:.1f} 年；< 2 年时有效独立样本不足，ubRMSE 差异不可辨。"))

    if inv.grd and inv.grd[0].get("crs"):
        groups: dict[str, int] = {}
        for r in inv.grd + inv.s2 + inv.pairs:
            if r.get("crs") and not r.get("error"):
                k = f"{r['crs']}|{r.get('res_m')}|{r.get('bounds')}"
                groups[k] = groups.get(k, 0) + 1
        if len(groups) <= 1:
            c(Check("PASS", "网格一致性", "所有栅格 CRS / 分辨率 / 范围一致。"))
        else:
            c(Check("FAIL", "网格一致性",
                    f"发现 {len(groups)} 组不同网格：{list(groups)[:3]}。"
                    "必须重投影到单一参考网格（同名点残差 < 0.5 像元）后再往下走。"))
    else:
        c(Check("NOTE", "网格一致性",
                "未读到栅格元数据（rasterio 不可用或目录为空）；"
                "装 rasterio 后重跑以核对 CRS / 分辨率 / 边界。"))

    vals = [r.get("vmean") for r in inv.grd if r.get("vmean") is not None]
    if vals:
        m = sum(vals) / len(vals)
        looks_db = -30 <= m <= 5
        c(Check("NOTE", "GRD 辐射口径",
                f"抽样均值 {m:.2f}，看起来是{'dB' if looks_db else '线性功率'}。"
                + ("需 10**(dB/10) 转线性后再进 Dubois/WCM 合成库。"
                   if looks_db else "若应为 dB，请检查 GEE 导出时是否已做转换。")))
    c(Check("NOTE", "辐射口径",
            f"配置声明 grd_radiometry={inp.get('grd_radiometry')}；"
            "仓库 dubois_sigma0 预测的是 σ0 线性功率，GEE GRD 的 σ0 与其一致；"
            "若混用 HyP3 RTC γ0，需 γ0≈σ0/cosθ 转换，两者不可直接拼进同一模型。"))

    n_pairs = len(inv.pairs)
    c(Check("PASS" if n_pairs >= l2.get("min_pairs", 30) else "WARN", "干涉对数量",
            f"{n_pairs} 对；L2 建议 ≥ {l2.get('min_pairs', 30)}。"))
    dts = [r["dt_days"] for r in inv.pairs if isinstance(r.get("dt_days"), int)]
    if dts:
        bad = [d for d in dts if d > l2.get("max_pair_baseline_days", 24) or d <= 0]
        c(Check("PASS" if not bad else "WARN", "时间基线",
                f"范围 {min(dts)}–{max(dts)} 天；{len(bad)} 对超出 "
                f"{l2.get('max_pair_baseline_days', 24)} 天或非正，建议剔除。"))
    c(Check("NOTE" if l2.get("mode") == "pairwise" else "WARN", "L2 模式",
            "HyP3 给的是逐干涉对 γ_ij，不是逐时相 γ(t)。"
            "pairwise 模式请用 calibrate_decay_coefficients_pairwise 等新函数；"
            "直接把 pair 堆成 (n_time,...) 喂给 invert_sm_from_coherence 会得到无物理意义的结果。"))

    # ---- 时间支撑对齐（逐对 vs 月度） ----
    grans = [r.get("temporal_granularity") for r in inv.grd + inv.s2]
    n_monthly = sum(1 for g in grans if g == "monthly")
    n_dated = sum(1 for g in grans if g)
    if n_monthly:
        c(Check("WARN", "时间支撑错配",
                f"{n_monthly}/{n_dated} 个光学/GRD 时相是**月度合成**（30 天状态量），"
                "而相干性是 12 天**差分**观测量。月值既不是端点值也不是 ΔV，"
                "且月内多个干涉对会共享同一月值。必须先跑 00b_align_temporal.py，"
                "按路线 B（重导逐景）/ A（γ 降采样到月 + 月际差分）对齐。"))
    else:
        c(Check("NOTE", "时间支撑",
                "未检出月度产品；仍须确认 S-2 / GRD 是逐景而非按月命名但被误判。"
                "逐景时用 pair_endpoint_difference() 插值到干涉对端点。"))

    if dts:
        if len(set(dts)) == 1:
            c(Check("WARN", "k_t 可识别性",
                    f"全部 {len(dts)} 对 Δt 都是 {dts[0]} 天：k_t·Δt 与截距 c 完全共线，"
                    "**必须固定 k_t=0**（或删去该项），否则标定出的系数无物理意义。"
                    "HyP3 配对策略由用户自选，若想要可变 Δt 需改用 nearest-neighbor depth=2 或 SBAS 混合基线。"))
        else:
            c(Check("PASS", "k_t 可识别性",
                    f"Δt 有 {len(set(dts))} 种取值（{min(dts)}–{max(dts)} 天），"
                    "时间与含水量去相干可分离，保留 k_t 项。"))

    cross = 0
    n_pair_dated = 0
    for p in inv.pairs:
        if p.get("ref_date") and p.get("sec_date"):
            n_pair_dated += 1
            if p["ref_date"][:7] != p["sec_date"][:7]:
                cross += 1
    if n_pair_dated:
        c(Check("NOTE", "跨月干涉对",
                f"{cross}/{n_pair_dated} 对跨月（{100*cross/n_pair_dated:.0f}%）。"
                "月度合成下月内对的 ΔV 恒为 0，是纯噪声样本；路线 C 只用跨月对。"))
    c(Check("NOTE", "闭合相位质控",
            "单对 INSAR_GAMMA 产品不含闭合相位（需三景三元组），本轮跳过该质控，"
            "论文中须明确声明。"))

    gaps = [r["gap_days"] for r in inv.match_rows if isinstance(r["gap_days"], int)]
    if gaps:
        med = sorted(gaps)[len(gaps) // 2]
        c(Check("PASS" if med <= 5 else "WARN", "S-1 / S-2 时相配对",
                f"中位间隔 {med} 天；> 6 天说明光学有效观测不足，"
                "须把「距最近有效光学观测天数」作为显式特征，且不要做月度最大值合成。"))

    c(Check("NOTE", "站点分层",
            "OzNet YA（约 13 站）为灌溉作物、YB（约 11 站）为旱草地，必须分层报告；"
            "无灌溉记录时事件响应检验会把灌溉误判为降水响应。"))
    c(Check("NOTE", "真值类型",
            "确认「公开土壤湿度观测产品」是 OzNet 站点观测还是 SMAP/SMOS 栅格产品。"
            "前者可作为标签与验证真值；后者只能做三重配置 / 相对比较，不能当标签训练。"))
    c(Check("NOTE", "L3 / L4",
            "本轮关闭：单对干涉产品不构成时序形变，Yanco 也非沉降区、无水位井。"
            "聚焦 L1 + L2 与验证纪律。"))


def write_outputs(inv: Inventory, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, rows in (("manifest_grd.csv", inv.grd),
                       ("manifest_s2.csv", inv.s2),
                       ("manifest_insar_pairs.csv", inv.pairs),
                       ("match_table.csv", inv.match_rows)):
        if not rows:
            continue
        fields: list[str] = []
        for r in rows:
            for k in r:
                if k not in fields:
                    fields.append(k)
        with (out_dir / name).open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    icon = {"PASS": "OK", "WARN": "!!", "FAIL": "XX", "NOTE": "--"}
    lines = ["# Yanco 数据清点报告", "",
             f"生成时间：{datetime.now():%Y-%m-%d %H:%M}", "",
             f"- GRD：{len(inv.grd)} 景", f"- S-2：{len(inv.s2)} 景",
             f"- 干涉对：{len(inv.pairs)} 对", "", "## 检查结果", ""]
    for c in inv.checks:
        lines.append(f"- `[{icon.get(c.level, '--')}] {c.level} · {c.item}` {c.message}")
    lines += ["", "## 下一步", "",
              "1. 先消除所有 FAIL 项（尤其网格一致性与时间基线）。",
              "2. **运行 scripts/00b_align_temporal.py**：把逐干涉对与逐景/月度光学、GRD "
              "落到同一张对齐表，产出 pair_alignment.csv 与 alignment_report.md。",
              "3. 光学/GRD 为月度合成时，按报告决定路线 B（重导逐景）或 A（γ 降采样到月 + 月际差分）。",
              "4. 再跑 scripts/01_prepare_s2.py → 02 → 03 → 04（L1）。",
              "5. L2 走 pairwise 模式，见 docs/yanco_workflow.md 与 docs/temporal_alignment.md。", ""]
    (out_dir / "inventory_report.md").write_text("\n".join(lines), encoding="utf-8")


def selftest() -> int:
    """构造合成目录与文件名，验证清点/配对链路（无需真实数据）。"""
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        for sub in ("s1_grd", "s2_l2a", "hyp3_insar"):
            (root / f"data_raw/yanco/{sub}").mkdir(parents=True)
        d0 = datetime(2021, 1, 1)
        grd_dir = root / "data_raw/yanco/s1_grd"
        s2_dir = root / "data_raw/yanco/s2_l2a"
        pair_dir = root / "data_raw/yanco/hyp3_insar"
        for i in range(0, 40):
            d = d0 + timedelta(days=12 * i)
            (grd_dir / f"S1A_IW_GRD_{d:%Y%m%d}_VV.tif").touch()
            (s2_dir / f"S2_L2A_NDVI_{d + timedelta(days=2):%Y%m%d}.tif").touch()
            (pair_dir / f"S1AA_{d:%Y%m%d}T200000_{(d + timedelta(days=12)):%Y%m%d}T200000_VV_012d_INT80_corr.tif").touch()
        cfg = load_config(Path(__file__).resolve().parents[1] / "configs" / "yanco.yaml")
        cfg["paths"] = {**cfg.get("paths", {}),
                        "s1_grd": "data_raw/yanco/s1_grd",
                        "s2_l2a": "data_raw/yanco/s2_l2a",
                        "insar_pairs": "data_raw/yanco/hyp3_insar"}
        inv = Inventory(root=root)
        build_manifests(inv, cfg)
        build_match_table(inv)
        run_checks(inv, cfg)
        write_outputs(inv, root / "outputs/yanco")

        assert len(inv.grd) == 40, f"GRD 清点错误：{len(inv.grd)}"
        assert len(inv.s2) == 40, f"S-2 清点错误：{len(inv.s2)}"
        assert len(inv.pairs) == 40, f"干涉对清点错误：{len(inv.pairs)}"
        assert all(p["dt_days"] == 12 for p in inv.pairs), "时间基线解析错误"
        assert all(r["gap_days"] == 2 for r in inv.match_rows), "时相配对错误"
        assert (root / "outputs/yanco/inventory_report.md").exists()
        assert (root / "outputs/yanco/manifest_insar_pairs.csv").exists()
        # 月度合成文件名（YYYYMM / YYYY_MM）也要能解析，且 6 位时间戳不能被误判
        d1, g1 = parse_scene_date("S2_NDVI_202103.tif")
        d2, g2 = parse_scene_date("S2_NDVI_2021_03.tif")
        d3, g3 = parse_scene_date("S2_NDVI_20210315.tif")
        d4, g4 = parse_scene_date("S1_T200000_corr.tif")
        assert (d1, g1) == (datetime(2021, 3, 1), "monthly")
        assert (d2, g2) == (datetime(2021, 3, 1), "monthly")
        assert g3 == "scene" and d3 == datetime(2021, 3, 15)
        assert (d4, g4) == (None, "")
        print("selftest OK: 40 GRD / 40 S-2 / 40 干涉对，dt=12d，光学间隔=2d；YYYYMM 解析正常")
        for c in inv.checks[:4]:
            print(f"  [{c.level}] {c.item}: {c.message[:70]}")
        return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Yanco 数据清点与口径核对")
    ap.add_argument("--root", default=".", help="项目根目录（配置文件里的相对路径以此为基准）")
    ap.add_argument("--config", default="configs/yanco.yaml")
    ap.add_argument("--out", default=None, help="输出目录，默认取配置 paths.output")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    root = Path(args.root).resolve()
    cfg = load_config(root / args.config if not Path(args.config).is_absolute()
                      else Path(args.config))
    inv = Inventory(root=root)
    build_manifests(inv, cfg)
    build_match_table(inv)
    run_checks(inv, cfg)
    out = Path(args.out) if args.out else root / cfg.get("paths", {}).get("output", "outputs/yanco")
    write_outputs(inv, out)
    print(f"清点完成：GRD {len(inv.grd)} / S-2 {len(inv.s2)} / 干涉对 {len(inv.pairs)}")
    for c in inv.checks:
        if c.level in ("FAIL", "WARN"):
            print(f"  [{c.level}] {c.item}: {c.message}")
    print(f"报告：{out / 'inventory_report.md'}")
    return 1 if any(c.level == "FAIL" for c in inv.checks) else 0


if __name__ == "__main__":
    raise SystemExit(main())
