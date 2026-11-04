"""Yanco 第 0b 步：把「逐干涉对」与「月度光学 / GRD」落到同一张对齐表。

为什么需要这一步
----------------
HyP3 相干性是逐干涉对（12 天差分），GEE 的 S-2 / GRD 是月度合成（30 天状态）。
本脚本不改动任何原始数据，只产出：

    pair_alignment.csv   每个干涉对一张时间对齐户口本
    alignment_report.md  时间支撑诊断、k_t 可识别性、路线建议

用法
----
    python scripts/00b_align_temporal.py --root .
    python scripts/00b_align_temporal.py --station-csv data_raw/yanco/oznet/YA3.csv \
        --time-col timestamp --value-col sm_0_5cm
    python scripts/00b_align_temporal.py --selftest

依赖：numpy；不依赖 rasterio / GEE。输入是 00_prepare_yanco.py 的 manifest CSV。
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ecohydro_sm.preprocess.temporal_alignment import (  # noqa: E402
    make_pair_dates,
    temporal_support_diagnostic,
)

DATE_FMT = "%Y-%m-%d"


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def _parse(s: str):
    if not s:
        return None
    for fmt in (DATE_FMT, "%Y%m%d", "%Y-%m-%d %H:%M:%S", "%Y/%m/%d"):
        try:
            return datetime.strptime(s.strip()[:19], fmt)
        except ValueError:
            continue
    return None


def _representative_day(d: datetime, monthly: bool, day: int = 15) -> datetime:
    """月度合成产品的代表日期取月中（默认 15 日）。

    GEE 月度合成常把 system:time_start 设为月首或首个输入景时间；若把月首
    当作"观测日"去插值，会凭空引入最多 30 天偏差。统一取月中并显式记录，
    是三种 evils 里最小、也最容易在论文里说清楚的一种。
    """
    if not monthly:
        return d
    return datetime(d.year, d.month, min(day, 28))


def _detect_monthly(rows: list[dict]) -> bool:
    """优先用 00 脚本写下的 ``temporal_granularity``；缺省时按 YYYY-MM 去重回退。

    不能只靠「每月恰好一个文件」判断：云多的月份逐景 S-2 也可能每月只剩一景，
    那样会被误判成月度合成。
    """
    grans = [r.get("temporal_granularity") for r in rows
             if r.get("temporal_granularity")]
    if grans:
        return all(g == "monthly" for g in grans)
    ms = [r.get("year_month") for r in rows if r.get("year_month")]
    return bool(ms) and len(set(ms)) == len(ms)


def build_alignment(pair_rows: list[dict], s2_rows: list[dict],
                    grd_rows: list[dict], cfg: dict) -> dict:
    ta = cfg.get("temporal_alignment", {})
    max_gap = float(ta.get("max_optical_gap_days", 10))
    rep_day = int(ta.get("month_representative_day", 15))
    assign_by = ta.get("monthly_assign_by", "sec")

    s2_monthly = _detect_monthly(s2_rows)
    grd_monthly = _detect_monthly(grd_rows)

    def prep(rows, monthly):
        out = []
        for r in rows:
            raw = (r.get("date") or "").strip()
            if not raw and r.get("year_month"):
                raw = f"{r['year_month']}-01"
            d = _parse(raw)
            if d is None:
                continue
            out.append((_representative_day(d, monthly, rep_day), r))
        return sorted(out, key=lambda x: x[0])

    s2 = prep(s2_rows, s2_monthly)
    grd = prep(grd_rows, grd_monthly)

    valid_pairs = []
    for r in pair_rows:
        a = _parse(r.get("ref_date", ""))
        b = _parse(r.get("sec_date", ""))
        if a and b and b > a:
            valid_pairs.append((r, a, b))
    if not valid_pairs:
        raise ValueError("manifest_insar_pairs.csv 里没有可用的 (ref_date, sec_date)。")

    refs = [x[1] for x in valid_pairs]
    secs = [x[2] for x in valid_pairs]
    dt = [float(r.get("dt_days") or (b - a).days) for r, a, b in valid_pairs]
    pd = make_pair_dates(refs, secs, dt_days=dt)

    def nearest(target: datetime, series):
        if not series:
            return "", None
        d, row = min(series, key=lambda x: abs((x[0] - target).days))
        return d.strftime(DATE_FMT), abs((d - target).days)

    key_all = {"ref": pd.ref_month, "sec": pd.sec_month,
               "mid": pd.mid_month}[assign_by]
    month_index = {m: i for i, m in enumerate(sorted(set(key_all)))}
    rows_out = []
    for i, (raw, a, b) in enumerate(valid_pairs):
        s2r, g_ref = nearest(a, s2)
        s2s, g_sec = nearest(b, s2)
        gr, gg_ref = nearest(a, grd)
        gs, gg_sec = nearest(b, grd)
        key_month = key_all[i]
        gap_max_s2 = max(v for v in (g_ref, g_sec) if v is not None) if s2 else None
        rows_out.append({
            "pair_id": i,
            "ref_date": a.strftime(DATE_FMT),
            "sec_date": b.strftime(DATE_FMT),
            "dt_days": int(pd.dt_days[i]),
            "ref_month": str(pd.ref_month[i]),
            "sec_month": str(pd.sec_month[i]),
            "cross_month": int(bool(pd.cross_month_mask[i])),
            "assign_month": str(key_month),
            "monthly_index": month_index.get(key_month, ""),
            "s2_ref_date": s2r, "s2_sec_date": s2s,
            "s2_gap_ref_days": g_ref if g_ref is not None else "",
            "s2_gap_sec_days": g_sec if g_sec is not None else "",
            "s2_usable": int(gap_max_s2 is not None and gap_max_s2 <= max_gap),
            "grd_ref_date": gr, "grd_sec_date": gs,
            "grd_gap_ref_days": gg_ref if gg_ref is not None else "",
            "grd_gap_sec_days": gg_sec if gg_sec is not None else "",
            "usable_for_A": int(True),                    # 月度聚合总有归属
            "usable_for_C": int(bool(pd.cross_month_mask[i])),
            # 月度产品的"代表日"是人为约定，不是观测日：即使离端点很近，
            # 也绝不能当端点值插值。B 路线在月度输入下恒不可用。
            "usable_for_B": int((not s2_monthly) and gap_max_s2 is not None
                                and gap_max_s2 <= max_gap),
            "file": raw.get("file", ""),
        })

    return {
        "rows": rows_out,
        "pair_dates": pd,
        "s2_monthly": s2_monthly,
        "grd_monthly": grd_monthly,
        "max_gap": max_gap,
        "rep_day": rep_day,
        "n_s2": len(s2), "n_grd": len(grd),
    }


def diagnose(al: dict, station_csv: Path | None, time_col: str,
             value_col: str, r2_threshold: float) -> dict | None:
    if station_csv is None or not station_csv.exists():
        return None
    rows = _read_csv(station_csv)
    ts, vs = [], []
    for r in rows:
        d = _parse(r.get(time_col, ""))
        try:
            v = float(r.get(value_col, ""))
        except (TypeError, ValueError):
            continue
        if d is not None and np.isfinite(v):
            ts.append(d)
            vs.append(v)
    if len(ts) < 30:
        return {"recommendation": f"站点有效记录仅 {len(ts)} 条，不足以做时间支撑诊断。"}
    return temporal_support_diagnostic(
        np.asarray(vs), ts,
        [datetime.strptime(r["ref_date"], DATE_FMT) for r in al["rows"]],
        [datetime.strptime(r["sec_date"], DATE_FMT) for r in al["rows"]],
        r2_threshold=r2_threshold,
    )


def write_outputs(al: dict, diag: dict | None, cfg: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = al["rows"]
    if rows:
        with (out_dir / "pair_alignment.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    pd = al["pair_dates"]
    n = len(rows)
    n_cross = sum(r["cross_month"] for r in rows)
    n_b = sum(r["usable_for_B"] for r in rows)
    gaps = [r["s2_gap_sec_days"] for r in rows
            if isinstance(r["s2_gap_sec_days"], (int, float))]
    med_gap = float(np.median(gaps)) if gaps else float("nan")

    L = ["# 时间支撑对齐报告", "",
         f"生成时间：{datetime.now():%Y-%m-%d %H:%M}", "",
         "## 输入时间粒度", "",
         f"- 干涉对：{n} 对，Δt 范围 {int(pd.dt_days.min())}–{int(pd.dt_days.max())} 天",
         f"- S-2：{al['n_s2']} 个时相，判定为**{'月度合成' if al['s2_monthly'] else '逐景'}**",
         f"- GRD：{al['n_grd']} 个时相，判定为**{'月度合成' if al['grd_monthly'] else '逐景'}**",
         f"- 月度产品代表日期：每月 {al['rep_day']} 日（显式约定，不是月首）", "",
         "## 关键判据", ""]

    if al["s2_monthly"] or al["grd_monthly"]:
        L += ["- **[!!] 时间支撑错配**：相干性是 12 天**差分**，光学/GRD 是 30 天**状态**。"
              "月值不是端点值，也不是 ΔV；直接按「干涉对所在月」取值会造成"
              "量纲错 + 月内对 ΔV=0 的重复样本 + 事件被平滑。", ""]

    if pd.dt_is_constant:
        L += [f"- **[!!] k_t 不可识别**：全部 {n} 对 Δt 均为 "
              f"{int(pd.dt_days[0])} 天，k_t·Δt 与截距 c 完全共线。"
              "**必须固定 k_t=0**（或删掉该项），否则 least_squares 会把任意值"
              "分配给 k_t 而 c 反向补偿，系数没有物理意义。"]
    else:
        L += ["- **[OK] k_t 可识别**：Δt 有变化，时间与含水量去相干可分离，保留 k_t 项。"]
    n_months = len({r["assign_month"] for r in rows})
    L += [f"- 跨月干涉对 {n_cross}/{n}（{100*n_cross/max(n,1):.0f}%）："
          "月内对在月度合成下 ΔV 恒为 0，路线 C 必须剔除。",
          f"- 副影像到最近光学代表日的中位间隔 {med_gap:.0f} 天；"
          f"距代表日在 ±{al['max_gap']:.0f} 天内的对 {n_b}/{n}。"
          "**该数字仅用于量化错配程度，不代表月值可作端点**（月中 15 日是人为约定，"
          "不是观测日）。", ""]

    L += ["## 路线可用性", "",
          "| 路线 | 可用样本 | 说明 |", "|---|---|---|",
          f"| A 月度降采样 | {n_months} 个月 | "
          "γ_ij 聚合成 γ_M，解释变量改**月际差分**，保留 Δt̄ 协变量；须重新标定 k |",
          f"| B 逐景插值到端点 | {sum(r['usable_for_B'] for r in rows)} 对 | "
          + ("**当前为月度输入，恒为 0**：必须先回 GEE 重导逐景产品" if al["s2_monthly"]
             else "逐景产品已就位，按 gap 窗口筛选") + " |",
          f"| C 仅跨月对 | {n_cross} 对 | 时间支撑仍错配，只能作敏感性分析 |", ""]

    if diag:
        L += ["## 站点时间支撑诊断（决定 A/C 是否可用）", ""]
        for k in ("r2", "rmse", "scale_ratio", "n_pairs",
                  "median_abs_pair_delta", "median_abs_monthly_delta"):
            if k in diag:
                L.append(f"- `{k}` = {diag[k]}")
        L += ["", f"> {diag.get('recommendation', '')}", ""]
    else:
        L += ["## 站点时间支撑诊断", "",
              "未执行（缺少 `--station-csv`）。**强烈建议用 OzNet 20 min 站点跑一次**：",
              "它能定量回答「月度合成丢掉了多少逐对信号」，这个数字可以直接写进论文。", ""]

    L += ["## 下一步", "",
          "1. 若诊断为 B（或诊断未跑且光学为月值）→ 回 GEE 重导**逐景** S-2 指数与 GRD，"
          "保留 `system:time_start`，不要 `calendarRange + mean/max`。",
          "2. 逐景产品到位后，用 `pair_endpoint_difference()` 求 |ΔV_ij|，"
          "喂给 `calibrate_decay_coefficients_pairwise()`。",
          "3. 站点取过境时刻值：`sar_acquisition_utc()` 处理 AEST/AEDT 夏令时，"
          "不要按文件名日期直接配对。",
          "4. A 路线作为稳健性检验与 B 并列报告，两者差异本身就是时间支撑的敏感度分析。", ""]
    (out_dir / "alignment_report.md").write_text("\n".join(L), encoding="utf-8")


def selftest() -> int:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        md = root / "outputs/yanco"
        md.mkdir(parents=True)
        d0 = datetime(2021, 1, 1)
        with (md / "manifest_insar_pairs.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["file", "ref_date", "sec_date", "dt_days"])
            w.writeheader()
            for i in range(12):
                a = d0 + timedelta(days=12 * i)
                w.writerow({"file": f"p{i}.tif", "ref_date": a.strftime(DATE_FMT),
                            "sec_date": (a + timedelta(days=12)).strftime(DATE_FMT),
                            "dt_days": 12})
        for name, prefix in (("manifest_s2.csv", "S2"), ("manifest_grd.csv", "GRD")):
            with (md / name).open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["file", "date", "year_month",
                                                  "temporal_granularity"])
                w.writeheader()
                for m in range(1, 5):
                    ym = f"2021-{m:02d}"
                    w.writerow({"file": f"{prefix}_{ym}.tif", "date": f"{ym}-01",
                                "year_month": ym,
                                "temporal_granularity": "monthly"})
        cfg = {"temporal_alignment": {"max_optical_gap_days": 10,
                                      "month_representative_day": 15,
                                      "monthly_assign_by": "sec"}}
        al = build_alignment(_read_csv(md / "manifest_insar_pairs.csv"),
                             _read_csv(md / "manifest_s2.csv"),
                             _read_csv(md / "manifest_grd.csv"), cfg)
        diag = None
        write_outputs(al, diag, cfg, root / "outputs/yanco")
        assert len(al["rows"]) == 12
        assert al["s2_monthly"] and al["grd_monthly"]
        assert al["pair_dates"].dt_is_constant is True
        assert 0 < sum(r["cross_month"] for r in al["rows"]) < 12
        # 月度输入下 B 路线必须恒不可用（代表日不是观测日，不能当端点）
        assert sum(r["usable_for_B"] for r in al["rows"]) == 0
        assert (root / "outputs/yanco/pair_alignment.csv").exists()
        txt = (root / "outputs/yanco/alignment_report.md").read_text(encoding="utf-8")
        assert "k_t 不可识别" in txt
        assert "恒为 0" in txt
        print("selftest OK: 12 干涉对 / 4 个月度时相；月值识别、Δt 共线判定、报告生成均通过")
        return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="干涉对 ↔ 月度遥感时间支撑对齐")
    ap.add_argument("--root", default=".")
    ap.add_argument("--manifest-dir", default=None, help="manifest 所在目录，默认 <root>/outputs/yanco")
    ap.add_argument("--config", default="configs/yanco.yaml")
    ap.add_argument("--out", default=None)
    ap.add_argument("--station-csv", default=None, help="高频站点 CSV（可选，用于时间支撑诊断）")
    ap.add_argument("--time-col", default="timestamp")
    ap.add_argument("--value-col", default="sm_0_5cm")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    root = Path(args.root).resolve()
    from ecohydro_sm.config import load_config
    cfg = load_config(root / args.config)
    md = Path(args.manifest_dir) if args.manifest_dir else root / cfg.get("paths", {}).get("output", "outputs/yanco")
    out = Path(args.out) if args.out else md
    try:
        al = build_alignment(_read_csv(md / "manifest_insar_pairs.csv"),
                             _read_csv(md / "manifest_s2.csv"),
                             _read_csv(md / "manifest_grd.csv"), cfg)
    except ValueError as exc:
        print(f"错误：{exc}\n请先运行 00_prepare_yanco.py 生成 manifest。")
        return 1
    diag = diagnose(al, Path(args.station_csv) if args.station_csv else None,
                    args.time_col, args.value_col,
                    float(cfg.get("temporal_alignment", {}).get("diagnostic_r2_threshold", 0.5)))
    write_outputs(al, diag, cfg, out)
    print(f"对齐完成：{len(al['rows'])} 对 → {out / 'pair_alignment.csv'}")
    print(f"报告：{out / 'alignment_report.md'}")
    if diag:
        print(f"诊断：{diag.get('recommendation', '')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
