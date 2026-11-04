"""时间支撑对齐：逐干涉对 γ_ij ↔ 逐景 / 月度光学与 SAR GRD。

问题
----
HyP3 ``INSAR_GAMMA`` 的相干性是**逐干涉对**的差分观测量 γ_ij，时间支撑是
[t_i, t_j]（通常 12 天）；而 GEE 导出的 S-2 / S-1 GRD 常被做成**月度合成**，
时间支撑是约 30 天的状态量。二者不在同一时间轴上。

直接按「干涉对所在月份」取月值会造成三类错误：

1. **量纲错**：γ_ij 响应的是 ΔV = V(t_j) − V(t_i)（两端点的差分），
   月值给的是 ⟨V⟩_M（月内平均 / 最大状态），既不是差分也不是端点值；
2. **解释变量重复**：一个月约 2–3 个干涉对。月内对共享同一个月值，
   ΔV 被算成 0，出现「同一 x、不同 y」，回归稀释使 k_sm / k_v 系统性偏低；
3. **事件被平滑掉**：月度合成（尤其最大值合成 MVC）抹去灌溉 / 降水脉冲，
   而 12 天相干性对这类脉冲最敏感 —— 恰好丢掉 L2 相对 L1 的增量信息。

三条对齐路线（按信息保真度排序）
--------------------------------
B（主结果）**逐景插值到端点**：在 GEE 重导逐景产品，插值出 V(t_i)、V(t_j)，
           再做端点差分。物理正确，保住事件响应。
A（稳健性）**月度降采样**：把 γ_ij 聚合成月度 γ_M，解释变量改用**月际差分**，
           并保留 Δt̄ 协变量。统计诚实，但月内事件丢失，不能报逐日精度。
C（敏感性）**只用跨月对**：以月际差分近似 ΔV，样本量约减半，时间支撑仍错配，
           只能作补充，不能作主结果。

选哪条不是拍脑袋：``temporal_support_diagnostic()`` 用高频站点观测量化
「月度代理能解释多少逐对 ΔSM」，R² 低就必须走 B。

另见 ``docs/temporal_alignment.md``。
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone

import numpy as np

__all__ = [
    "PairDates",
    "make_pair_dates",
    "interpolate_series_to_dates",
    "pair_endpoint_difference",
    "select_cross_month_pairs",
    "aggregate_pairs_to_monthly",
    "monthly_state_to_difference",
    "temporal_support_diagnostic",
    "sar_acquisition_utc",
]


# ---------------------------------------------------------------------------
# 日期工具（全部内部转成"自 1970-01-01 的天数"或"自 1970-01 的月数"）
# ---------------------------------------------------------------------------
def _date_string(d) -> str:
    if isinstance(d, (datetime, date)):
        return d.strftime("%Y-%m-%d")
    if isinstance(d, np.datetime64):
        return str(np.datetime64(d, "D"))
    return str(d)[:10]


def _as_day_number(dates) -> np.ndarray:
    out = [float(np.datetime64(_date_string(d)[:10], "D").astype("int64"))
           for d in dates]
    return np.asarray(out, dtype=float)


def _as_month_number(dates) -> np.ndarray:
    """月份序号（自 1970-01 起的月数）。

    整数数组视为**已经是**月份序号（``aggregate_pairs_to_monthly`` 的
    ``months`` 返回值就是这种）；其余（字符串 / datetime / np.datetime64）
    按日期解析。
    """
    vals = np.asarray(dates).ravel()
    if np.issubdtype(vals.dtype, np.integer):
        return vals.astype(np.int64)
    out = [int(np.datetime64(_date_string(d)[:7], "M").astype("int64"))
           for d in vals]
    return np.asarray(out, dtype=np.int64)


@dataclass
class PairDates:
    """一张干涉对时间表：所有下游对齐都从它出发。"""

    ref_days: np.ndarray       # (n_pairs,) 参考影像日期（天数）
    sec_days: np.ndarray       # (n_pairs,) 副影像日期
    dt_days: np.ndarray        # (n_pairs,) 时间基线，天
    mid_days: np.ndarray       # (n_pairs,) 中点日期
    ref_month: np.ndarray      # (n_pairs,) 参考影像月份序号
    sec_month: np.ndarray      # (n_pairs,) 副影像月份序号
    mid_month: np.ndarray      # (n_pairs,) 中点月份序号

    @property
    def n_pairs(self) -> int:
        return int(self.ref_days.size)

    @property
    def dt_is_constant(self) -> bool:
        """Δt 全一致时 k_t 与截距 c 完全共线，**不可识别**，必须固定 k_t=0。"""
        return bool(np.allclose(self.dt_days, self.dt_days[0])) if self.n_pairs else True

    @property
    def cross_month_mask(self) -> np.ndarray:
        return self.ref_month != self.sec_month


def make_pair_dates(ref_dates, sec_dates, dt_days=None) -> PairDates:
    """由参考 / 副影像日期构建干涉对时间表。

    ``dt_days`` 缺省由日期差计算；若提供了 HyP3 文件名里的 ``nnn`` 字段
    （时间 separation in days），以它为准更可靠。
    """
    r = _as_day_number(ref_dates)
    s = _as_day_number(sec_dates)
    if r.size != s.size:
        raise ValueError("ref_dates 与 sec_dates 长度不一致。")
    dt = (s - r) if dt_days is None else np.asarray(dt_days, dtype=float).ravel()
    if dt.size != r.size:
        raise ValueError("dt_days 长度与干涉对数量不一致。")
    if np.any(dt <= 0):
        raise ValueError("存在非正的时间基线，请检查参考/副影像顺序。")
    mid = 0.5 * (r + s)
    return PairDates(
        ref_days=r, sec_days=s, dt_days=dt, mid_days=mid,
        ref_month=_as_month_number(ref_dates),
        sec_month=_as_month_number(sec_dates),
        # 中点月份由天数反推，避免对浮点天数直接取月
        mid_month=np.asarray(
            [int(np.datetime64(date(1970, 1, 1) + timedelta(days=float(m)), "M")
                 .astype("int64")) for m in mid], dtype=np.int64),
    )


# ---------------------------------------------------------------------------
# 路线 B：逐景插值到干涉对端点
# ---------------------------------------------------------------------------
def interpolate_series_to_dates(
    scene_values: np.ndarray,
    scene_dates,
    target_dates,
    max_gap_days: float = 10.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """把逐景时序线性插值到任意目标日期（如干涉对的两个端点）。

    参数
    ----
    scene_values : (n_scenes, ...) 逐景指数 / 幅度，第一维是时间
    scene_dates : (n_scenes,) 每景日期（**必须是该景的实际观测日**，
                  月度合成不能冒充逐景日期）
    target_dates : (n_targets,) 目标日期
    max_gap_days : 目标日期到最近有效景的最大允许间隔；超出则该目标置 NaN

    返回
    ----
    values : (n_targets, ...) 插值结果
    gap_days : (n_targets,) 到最近景的间隔（不区分像元是否有效）
    valid : (n_targets,) 时间窗口是否合格

    说明
    ----
    云洞造成的像元级 NaN：插值时若一侧无效则退化为另一侧最近值，两侧都
    无效则 NaN。``gap_days`` 只反映时间上的最近景距离，**不反映该景在某一
    像元上是否有效**；像元级有效性请在插值前用云掩膜处理，或在下游用
    ``np.isfinite`` 逐像元屏蔽（pairwise 标定函数已内置该屏蔽）。
    """
    sv = np.asarray(scene_values, dtype=float)
    if sv.ndim == 1:
        sv = sv[:, None]
    sd = _as_day_number(scene_dates)
    td = _as_day_number(target_dates)
    if sv.shape[0] != sd.size:
        raise ValueError("scene_values 第一维长度与 scene_dates 不一致。")

    order = np.argsort(sd)
    sd, sv = sd[order], sv[order]
    n = sd.size

    pos = np.searchsorted(sd, td, side="left")
    li = np.clip(pos - 1, 0, n - 1)
    ri = np.clip(pos, 0, n - 1)
    dl = np.abs(td - sd[li])
    dr = np.abs(sd[ri] - td)
    gap = np.minimum(dl, dr)

    vl, vr = sv[li], sv[ri]
    span = sd[ri] - sd[li]
    w = np.where(span > 0, (td - sd[li]) / np.where(span > 0, span, 1.0), 0.0)
    w = np.clip(w, 0.0, 1.0).reshape((-1,) + (1,) * (sv.ndim - 1))

    interp = (1.0 - w) * vl + w * vr
    both = np.isfinite(vl) & np.isfinite(vr)
    out = np.where(both, interp,
                   np.where(np.isfinite(vl), vl,
                            np.where(np.isfinite(vr), vr, np.nan)))
    valid = gap <= max_gap_days
    out[~valid] = np.nan
    return out, gap, valid


def pair_endpoint_difference(
    scene_values: np.ndarray,
    scene_dates,
    ref_dates,
    sec_dates,
    max_gap_days: float = 10.0,
) -> dict:
    """路线 B：由逐景光学 / GRD 求每个干涉对的端点差分。

    返回 ``delta``（带符号，可用来给 ΔSM 赋符号）、``abs_delta``（喂给
    pairwise 模型，模型用的是 |ΔV|）、两端点值、两端的光学间隔与按对的有效性。
    """
    vr, gr, qr = interpolate_series_to_dates(scene_values, scene_dates, ref_dates,
                                             max_gap_days=max_gap_days)
    vs, gs, qs = interpolate_series_to_dates(scene_values, scene_dates, sec_dates,
                                             max_gap_days=max_gap_days)
    delta = vs - vr
    return {
        "delta": delta,
        "abs_delta": np.abs(delta),
        "value_ref": vr,
        "value_sec": vs,
        "gap_ref_days": gr,
        "gap_sec_days": gs,
        "valid": qr & qs,
    }


# ---------------------------------------------------------------------------
# 路线 C：只用跨月干涉对
# ---------------------------------------------------------------------------
def select_cross_month_pairs(pair_dates: PairDates) -> np.ndarray:
    """返回跨月干涉对的布尔掩码。

    月内对在月度合成下 ΔV 恒为 0（同一个月值减自己），是纯噪声样本，
    必须剔除。跨月对虽仍有「30 天状态 vs 12 天差分」的时间支撑错配，
    但至少月际差分是一个真实的差分量。
    """
    return pair_dates.cross_month_mask


# ---------------------------------------------------------------------------
# 路线 A：把相干性降采样到月度
# ---------------------------------------------------------------------------
def aggregate_pairs_to_monthly(
    gamma_pairs: np.ndarray,
    pair_dates: PairDates,
    assign_by: str = "sec",
    reducer: str = "median",
    min_pairs: int = 2,
) -> dict:
    """路线 A：把逐干涉对 γ_ij 聚合成月度 γ_M。

    归属规则缺省按**副影像**月份（γ 反映「到 t_j 为止」的地表变化）。
    用中位数而非均值：γ 分布有偏，且长基线 / 雨日异常对会拉偏均值。

    返回
    ----
    months : (n_months,) 月份序号（可直接传给 ``monthly_state_to_difference``）
    gamma_monthly : (n_months, ...) 月度相干性
    gamma_std, n_pairs, dt_mean, valid
        ``gamma_std`` 是月内 γ 的离散度 —— 它本身就是诊断量：
        月内离散度大说明月度聚合把显著的月内变化平滑掉了。
        ``dt_mean`` 必须作为协变量保留：各月干涉对的平均时间基线不同时，
        不建模会把时间去相干误判成含水量变化。
    """
    if assign_by not in ("ref", "sec", "mid"):
        raise ValueError("assign_by 必须是 ref / sec / mid。")
    g = np.asarray(gamma_pairs, dtype=float)
    if g.shape[0] != pair_dates.n_pairs:
        raise ValueError("gamma_pairs 第一维长度与干涉对数量不一致。")

    key = {"ref": pair_dates.ref_month,
           "sec": pair_dates.sec_month,
           "mid": pair_dates.mid_month}[assign_by]
    months = np.unique(key)
    shape = (months.size,) + g.shape[1:]
    out = np.full(shape, np.nan)
    gamma_std = np.full(shape, np.nan)
    n_pairs = np.zeros(months.size, dtype=int)
    dt_mean = np.full(months.size, np.nan)
    std_mean = np.zeros(months.size, dtype=float)

    for i, m in enumerate(months):
        block = g[key == m]
        n_pairs[i] = block.shape[0]
        with np.errstate(all="ignore"):
            out[i] = (np.nanmedian(block, axis=0) if reducer == "median"
                      else np.nanmean(block, axis=0))
            gamma_std[i] = np.nanstd(block, axis=0)
        dt_mean[i] = float(np.nanmean(pair_dates.dt_days[key == m]))
        std_mean[i] = float(np.nanmean(np.abs(gamma_std[i])))

    return {
        "months": months,
        "gamma_monthly": out,
        "gamma_std": gamma_std,
        "within_month_std_mean": std_mean,
        "n_pairs": n_pairs,
        "dt_mean": dt_mean,
        "valid": n_pairs >= min_pairs,
        "contiguous": bool(np.all(np.diff(months) == 1)),
    }


def monthly_state_to_difference(
    monthly_values: np.ndarray,
    months=None,
    absolute: bool = False,
) -> np.ndarray:
    """路线 A/C 的必需转换：月度**状态量** → 月际**变化量**。

    月度合成给的是状态（本月平均 / 最大 NDVI），而 pairwise 模型要的是
    |ΔV|。直接把状态量喂进去是量纲错误；必须先做月际差分。

    ``months`` 提供时，缺月处的差分会被置 NaN（跨 2 个月的差分时间支撑
    与相邻月不同，不能混用）。首月没有前驱，恒为 NaN。
    """
    v = np.asarray(monthly_values, dtype=float)
    out = np.full_like(v, np.nan, dtype=float)
    if v.shape[0] < 2:
        return np.abs(out) if absolute else out
    out[1:] = np.diff(v, axis=0)
    if months is not None:
        m = _as_month_number(months)
        if m.size != v.shape[0]:
            raise ValueError("months 长度与 monthly_values 第一维不一致。")
        bad = np.nonzero(np.diff(m) != 1)[0] + 1
        out[bad] = np.nan
    return np.abs(out) if absolute else out


# ---------------------------------------------------------------------------
# 诊断：用高频站点决定走 A/C 还是必须重导逐景（B）
# ---------------------------------------------------------------------------
def temporal_support_diagnostic(
    sm_high_freq: np.ndarray,
    times_high_freq,
    ref_dates,
    sec_dates,
    r2_threshold: float = 0.5,
    assign_by: str = "sec",
) -> dict:
    """量化「月度代理」能解释多少逐干涉对的真实 ΔSM。

    做法：用 20 min / 逐时站点观测算出每个 pair 端点的真实 ΔSM_ij，
    再把同一序列聚成月均值、取月际差分作为「月度代理」，比较二者。

    判据用 **相关系数平方**（尺度无关的解释方差），不是「直接拿月度代理
    当逐对 ΔSM」的 R²：月度差分的时间支撑是 30 天、逐对差分是 12 天，
    二者幅度本就不同，固定尺度比较会把尺度失配误判成信息丢失。
    ``scale_ratio`` 单独报告这个失配，提醒 A 路线**必须重新标定 k_sm/k_v**，
    不能复用逐对系数。

    R² 低说明月度合成丢掉了相干性所响应的短周期变化 —— 此时 A / C
    路线的解释变量本质上没有信息，必须回 GEE 重导逐景产品（路线 B）。
    这个诊断本身可以写进论文，作为「为什么必须逐景对齐」的定量依据。

    返回
    ----
    dict：r2、rmse、scale_ratio、n_pairs、中位绝对幅度（真值 / 代理）、
          amplitude_attenuation、recommendation、strategy
    """
    s = np.asarray(sm_high_freq, dtype=float).ravel()
    t = _as_day_number(times_high_freq)
    if s.size != t.size:
        raise ValueError("sm_high_freq 与 times_high_freq 长度不一致。")
    order = np.argsort(t)
    s, t = s[order], t[order]

    pd = make_pair_dates(ref_dates, sec_dates)
    # 高频站点：端点当天必有观测，窗口设 1 天
    vr, _, qr = interpolate_series_to_dates(s[:, None], times_high_freq,
                                            ref_dates, max_gap_days=1.0)
    vs, _, qs = interpolate_series_to_dates(s[:, None], times_high_freq,
                                            sec_dates, max_gap_days=1.0)
    true_delta = (vs - vr).ravel()
    true_delta[~(qr & qs)] = np.nan

    # 月度代理
    mnum = np.asarray(
        [int(np.datetime64(date(1970, 1, 1) + timedelta(days=float(d)), "M")
             .astype("int64")) for d in t], dtype=np.int64)
    months = np.unique(mnum)
    monthly = np.asarray([np.nanmean(s[mnum == m]) for m in months])
    monthly_diff = monthly_state_to_difference(monthly[:, None], months,
                                               absolute=False).ravel()
    key = {"ref": pd.ref_month, "sec": pd.sec_month,
           "mid": pd.mid_month}[assign_by]
    idx = np.clip(np.searchsorted(months, key), 0, months.size - 1)
    proxy = np.where(months[idx] == key, monthly_diff[idx], np.nan)

    mask = np.isfinite(true_delta) & np.isfinite(proxy)
    if mask.sum() < 5:
        return {"r2": np.nan, "rmse": np.nan, "n_pairs": int(mask.sum()),
                "recommendation": "有效对不足 5 个，无法诊断；请检查站点时间覆盖。"}
    y, yhat = true_delta[mask], proxy[mask]
    with np.errstate(invalid="ignore"):
        corr = float(np.corrcoef(y, yhat)[0, 1])
    r2 = corr ** 2 if np.isfinite(corr) else np.nan
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    med_true = float(np.median(np.abs(y)))
    med_proxy = float(np.median(np.abs(yhat)))
    scale_ratio = med_proxy / med_true if med_true > 0 else np.nan
    atten = 1.0 - scale_ratio if np.isfinite(scale_ratio) else np.nan
    ok = bool(np.isfinite(r2) and r2 >= r2_threshold)
    r2_txt = f"{r2:.2f}" if np.isfinite(r2) else "nan"
    return {
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "rmse": rmse,
        "scale_ratio": float(scale_ratio) if np.isfinite(scale_ratio) else np.nan,
        "n_pairs": int(mask.sum()),
        "median_abs_pair_delta": med_true,
        "median_abs_monthly_delta": med_proxy,
        "amplitude_attenuation": float(atten) if np.isfinite(atten) else np.nan,
        "recommendation": (
            f"月度代理解释方差 R²={r2_txt} ≥ {r2_threshold:.2f}：路线 A/C 可作稳健性"
            "检验；但 scale_ratio≠1，月度模型须重新标定 k_sm/k_v，不能复用逐对系数。"
            if ok else
            f"月度代理解释方差 R²={r2_txt} < {r2_threshold:.2f}：月度合成无法解释逐对"
            "变化，必须回 GEE 重导逐景 S-2 / GRD，走路线 B（插值到干涉对端点）。"
        ),
        "strategy": "A_or_C" if ok else "B",
    }


# ---------------------------------------------------------------------------
# 站点过境时刻（正确处理澳洲夏令时，替代固定 UTC 偏移）
# ---------------------------------------------------------------------------
def sar_acquisition_utc(sar_date, local_hour: float = 6.5,
                        tz_name: str = "Australia/Sydney") -> datetime:
    """把 S-1 文件名里的 **UTC 日期** 换算成 UTC 过境时刻。

    降轨当地清晨约 6:30 = 澳洲东部时（AEST UTC+10 / AEDT UTC+11），
    因此站点当地日期 = 文件名 UTC 日期 **+ 1 天**。用 zoneinfo 处理夏令时：
    AEST 下 6:30 = 当日 20:30 UTC（文件名日期），
    AEDT 下 6:30 = 当日 19:30 UTC（文件名日期）。

    返回 tz-aware UTC datetime，可直接用来从 20 min 站点序列取过境时刻值。
    固定 +10 偏移会在夏令时期间错 1 小时，不要用。
    """
    from zoneinfo import ZoneInfo

    s = _date_string(sar_date)
    if len(s) >= 8 and s[:8].isdigit() and "-" not in s[:8]:
        s = f"{s[:4]}-{s[4:6]}-{s[6:8]}"   # 兼容 YYYYMMDD 紧凑格式
    d = date.fromisoformat(s[:10])
    local_day = d + timedelta(days=1)
    hh = int(local_hour)
    mm = int(round((local_hour - hh) * 60))
    local_dt = datetime(local_day.year, local_day.month, local_day.day,
                        hh, mm, tzinfo=ZoneInfo(tz_name))
    return local_dt.astimezone(timezone.utc)
