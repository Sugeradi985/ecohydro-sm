"""时间支撑对齐模块测试：三类时间轴（逐对 / 逐景 / 月度）不能混用。"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest

from ecohydro_sm.preprocess.temporal_alignment import (
    aggregate_pairs_to_monthly,
    interpolate_series_to_dates,
    make_pair_dates,
    monthly_state_to_difference,
    pair_endpoint_difference,
    sar_acquisition_utc,
    select_cross_month_pairs,
    temporal_support_diagnostic,
)

D0 = datetime(2021, 1, 1)


def _dates(offsets):
    return [D0 + timedelta(days=int(k)) for k in offsets]


# ---------------------------------------------------------------------------
# 干涉对时间表
# ---------------------------------------------------------------------------
def test_pair_dates_dt_constant_blocks_identification_of_k_t():
    """Δt 全为 12 天时，k_t·Δt 与截距 c 完全共线，必须固定 k_t=0。"""
    r = _dates([0, 12, 24, 36])
    s = _dates([12, 24, 36, 48])
    pd = make_pair_dates(r, s)
    assert pd.n_pairs == 4
    assert np.all(pd.dt_days == 12)
    assert pd.dt_is_constant is True
    # mid_days 是"自 1970-01-01 的天数"绝对值，用相对首对比较
    assert np.allclose(pd.mid_days - pd.ref_days[0], [6, 18, 30, 42])


def test_pair_dates_variable_dt_is_identifiable():
    """nearest-neighbor depth=2 会混出 12/24 天对，k_t 才可标定。"""
    r = _dates([0, 12, 24])
    s = _dates([12, 36, 36])
    pd = make_pair_dates(r, s)
    assert np.allclose(pd.dt_days, [12, 24, 12])
    assert pd.dt_is_constant is False
    with pytest.raises(ValueError, match="非正"):
        make_pair_dates(_dates([12]), _dates([0]))


def test_cross_month_mask_matches_manual_count():
    r = _dates([0, 12, 24, 36])       # 1/1, 1/13, 1/25, 2/6
    s = _dates([12, 24, 36, 48])      # 1/13, 1/25, 2/6, 2/18
    pd = make_pair_dates(r, s)
    mask = select_cross_month_pairs(pd)
    # 只有 [1/25, 2/6] 这一对跨月
    assert mask.tolist() == [False, False, True, False]


# ---------------------------------------------------------------------------
# 路线 B：逐景插值到端点
# ---------------------------------------------------------------------------
def test_interpolate_linear_weighting():
    scenes = np.array([0.0, 1.0])[:, None]
    sd = _dates([0, 10])
    v, gap, valid = interpolate_series_to_dates(scenes, sd, _dates([5]),
                                                max_gap_days=10)
    assert np.isclose(v[0, 0], 0.5)
    assert np.allclose(gap, [5])
    assert valid.all()


def test_interpolate_beyond_gap_window_is_masked():
    scenes = np.array([0.2, 0.3])[:, None]
    sd = _dates([0, 10])
    v, gap, valid = interpolate_series_to_dates(scenes, sd, _dates([40]),
                                                max_gap_days=10)
    assert np.isnan(v[0, 0])
    assert not valid[0]
    assert gap[0] == 30


def test_interpolate_cloud_hole_falls_back_to_nearest_valid_scene():
    """云洞造成的像元级 NaN：单侧无效时退化为另一侧，双侧无效才 NaN。"""
    scenes = np.array([[0.2, 0.5], [np.nan, 0.7]])
    sd = _dates([0, 10])
    v, _, _ = interpolate_series_to_dates(scenes, sd, _dates([5]),
                                          max_gap_days=10)
    assert np.isclose(v[0, 0], 0.2)          # 右侧无效 -> 用左侧
    assert np.isclose(v[0, 1], 0.6)          # 两侧有效 -> 线性插值


def test_pair_endpoint_difference_keeps_sign_and_abs():
    """端点差分保留符号（用来赋 ΔSM 符号），同时给出 |ΔV| 供 pairwise 模型。"""
    scenes = np.array([0.10, 0.30])[:, None]
    sd = _dates([0, 12])
    out = pair_endpoint_difference(scenes, sd, _dates([0]), _dates([12]),
                                   max_gap_days=10)
    assert np.isclose(out["delta"][0, 0], 0.20)
    assert np.isclose(out["abs_delta"][0, 0], 0.20)
    assert out["valid"].all()
    assert np.allclose(out["gap_ref_days"], [0])
    assert np.allclose(out["gap_sec_days"], [0])


# ---------------------------------------------------------------------------
# 路线 A：月度降采样
# ---------------------------------------------------------------------------
def test_monthly_aggregation_uses_median_and_tracks_pair_count():
    # sec 日期 = 1/1, 1/13, 1/25, 2/6, 2/18 -> 1 月 3 对、2 月 2 对
    r = _dates([-12, 0, 12, 24, 36])
    s = _dates([0, 12, 24, 36, 48])
    pd = make_pair_dates(r, s)
    # 1 月三对里放一个异常低值，中位数应压住它
    gamma = np.array([[0.8], [0.02], [0.7], [0.6], [0.5]])
    res = aggregate_pairs_to_monthly(gamma, pd, assign_by="sec", min_pairs=2)
    assert res["n_pairs"].tolist() == [3, 2]
    assert np.isclose(res["gamma_monthly"][0, 0], 0.7)
    assert np.isclose(res["gamma_monthly"][1, 0], 0.55)
    assert np.allclose(res["dt_mean"], [12, 12])
    assert res["valid"].all()
    assert res["contiguous"] is True
    # 月内离散度：1 月含异常对，应明显大于 2 月
    assert res["within_month_std_mean"][0] > res["within_month_std_mean"][1]


def test_monthly_state_to_difference_is_not_the_state_itself():
    """最容易犯的错：把月度状态量直接喂给要 |ΔV| 的模型。"""
    v = np.array([0.1, 0.3, 0.25])[:, None]
    d = monthly_state_to_difference(v, ["2021-01-01", "2021-02-01", "2021-03-01"])
    assert np.isnan(d[0, 0])
    assert np.isclose(d[1, 0], 0.2)
    assert np.isclose(d[2, 0], -0.05)
    assert np.allclose(monthly_state_to_difference(v, absolute=True)[1:].ravel(),
                       [0.2, 0.05])


def test_monthly_difference_masks_non_contiguous_months():
    """缺月处的差分跨两个月，时间支撑不同，必须置 NaN。"""
    v = np.array([0.1, 0.3, 0.4])[:, None]
    d = monthly_state_to_difference(v, ["2021-01-01", "2021-02-01", "2021-04-01"])
    assert np.isnan(d[2, 0])
    assert np.isclose(d[1, 0], 0.2)


# ---------------------------------------------------------------------------
# 诊断：决定 A/C 还是必须 B
# ---------------------------------------------------------------------------
def _daily_pairs(n_days=365, step=12):
    refs, secs = [], []
    for k in range(0, n_days - step, step):
        refs.append(D0 + timedelta(days=k))
        secs.append(D0 + timedelta(days=k + step))
    return refs, secs


def test_diagnostic_recommends_scene_level_when_events_are_short():
    """灌溉/降水脉冲（逐日白噪声）：月均抹平事件，月度代理无解释力 -> B。"""
    rng = np.random.default_rng(7)
    times = _dates(np.arange(365))
    sm = 0.15 + 0.05 * rng.standard_normal(365)
    refs, secs = _daily_pairs()
    out = temporal_support_diagnostic(sm, times, refs, secs, r2_threshold=0.5)
    assert out["n_pairs"] >= 20
    assert out["r2"] < 0.5, out
    assert out["strategy"] == "B"
    assert "重导逐景" in out["recommendation"]


def test_diagnostic_allows_monthly_when_signal_is_smooth_seasonal():
    """纯平滑季节变化：月度降采样不丢信息 -> A/C 可作稳健性检验。"""
    t = np.arange(365)
    times = _dates(t)
    sm = 0.15 + 0.06 * np.sin(2 * np.pi * t / 365.0)
    refs, secs = _daily_pairs()
    out = temporal_support_diagnostic(sm, times, refs, secs, r2_threshold=0.5)
    assert out["r2"] > 0.8, out
    assert out["strategy"] == "A_or_C"
    # 30 天差分 vs 12 天差分，幅度必然不同 -> 月度模型必须重新标定
    assert out["scale_ratio"] > 1.5
    assert "重新标定" in out["recommendation"]


# ---------------------------------------------------------------------------
# 站点过境时刻 / 夏令时
# ---------------------------------------------------------------------------
def test_sar_acquisition_utc_handles_daylight_saving():
    """AEST(+10) 与 AEDT(+11) 下 UTC 时刻相差 1 小时，固定 +10 会错。"""
    winter = sar_acquisition_utc("2021-07-15", local_hour=6.5)
    summer = sar_acquisition_utc("2021-01-15", local_hour=6.5)
    assert winter.strftime("%Y-%m-%d %H:%M") == "2021-07-15 20:30"
    assert summer.strftime("%Y-%m-%d %H:%M") == "2021-01-15 19:30"
    assert winter.tzinfo is not None
    # 兼容 S-1 文件名里的 YYYYMMDD 紧凑格式
    assert sar_acquisition_utc("20210715").strftime("%Y-%m-%d %H:%M") == \
        "2021-07-15 20:30"
