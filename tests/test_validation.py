"""验证模块测试：指标分解、空间块 CV、三重配置。"""

import numpy as np
import pytest

from ecohydro_sm.validation.metrics import (
    anomaly_r,
    bias,
    crmse,
    full_report,
    pearson_r,
    rmse,
    stratified_report,
    ubrmse,
)
from ecohydro_sm.validation.spatial_cv import (
    SpatialBlockKFold,
    assign_spatial_blocks,
    leave_one_station_out,
)
from ecohydro_sm.validation.triple_collocation import (
    error_correlation_check,
    tc_scaled,
    triple_collocation,
)


def _pair(n=500, bias_val=0.02, noise=0.03, seed=0):
    rng = np.random.default_rng(seed)
    r = rng.uniform(0.05, 0.40, n)
    p = r + bias_val + rng.normal(0, noise, n)
    return p, r


def test_bias_separation():
    p, r = _pair(bias_val=0.02, noise=0.0)
    assert np.isclose(bias(p, r), 0.02, atol=1e-6)
    assert np.isclose(ubrmse(p, r), 0.0, atol=1e-6), "无随机误差时 ubRMSE 应为 0"
    assert rmse(p, r) > ubrmse(p, r)


def test_ubrmse_removes_bias():
    p, r = _pair(bias_val=0.05, noise=0.02)
    assert ubrmse(p, r) < rmse(p, r), "ubRMSE 必须小于 RMSE"


def test_anomaly_r_detects_dynamics():
    t = np.linspace(0, 6 * np.pi, 300)
    r = 0.2 + 0.1 * np.sin(t)
    p_good = 0.2 + 0.1 * np.sin(t) + np.random.default_rng(1).normal(0, 0.01, 300)
    p_bad = 0.2 + 0.1 * np.sin(t + np.pi)  # 相位相反
    assert anomaly_r(p_good, r) > 0.9
    assert anomaly_r(p_bad, r) < 0.0


def test_full_report_keys():
    p, r = _pair()
    rep = full_report(p, r)
    for k in ("n", "bias", "rmse", "ubrmse", "crmse", "pearson_r", "anomaly_r"):
        assert k in rep


def test_stratified_report_covers_bins():
    p, r = _pair(n=800)
    vwc = np.random.default_rng(2).uniform(0, 5, 800)
    rep = stratified_report(p, r, vwc, bins=[0, 0.5, 1.5, 3.0, 99.0])
    assert len(rep) == 4
    assert all("ubrmse" in v for v in rep.values())


def test_spatial_blocks_and_cv():
    n = 600
    x = np.random.default_rng(3).uniform(0, 30_000, n)
    y = np.random.default_rng(4).uniform(0, 30_000, n)
    blocks = assign_spatial_blocks(x, y, block_size_m=10_000)
    assert blocks.max() >= 3, "分块数过少"
    cv = SpatialBlockKFold(blocks, n_splits=4)
    folds = list(cv.split(np.zeros((n, 1))))
    assert len(folds) == 4
    for tr, te in folds:
        assert len(tr) > 0 and len(te) > 0
        assert len(np.intersect1d(blocks[tr], blocks[te])) == 0, "训练与测试块必须不重叠"


def test_leave_one_station_out():
    st = np.array([0, 0, 1, 1, 2, 2])
    folds = list(leave_one_station_out(st))
    assert len(folds) == 3
    assert all(len(te) == 2 for _, te in folds)


def test_triple_collocation_recovers_error():
    n = 4000
    rng = np.random.default_rng(5)
    truth = rng.uniform(0.05, 0.40, n)
    a = truth + rng.normal(0, 0.04, n)
    b = truth + rng.normal(0, 0.05, n)
    c = truth + rng.normal(0, 0.06, n)
    tc = triple_collocation(a, b, c)
    assert abs(tc["rmse_x"] - 0.04) < 0.01
    assert abs(tc["rmse_y"] - 0.05) < 0.01
    assert abs(tc["rmse_z"] - 0.06) < 0.01


def test_tc_independence_check_passes_for_independent_errors():
    n = 4000
    rng = np.random.default_rng(6)
    truth = rng.uniform(0.05, 0.40, n)
    a = truth + rng.normal(0, 0.04, n)
    b = truth + rng.normal(0, 0.04, n)
    c = truth + rng.normal(0, 0.04, n)
    res = error_correlation_check(a, b, c)
    assert res["independence_ok"], f"独立误差被误判：{res['deviation']}"
    # 完全独立时观测相关应接近理论期望（等方差下约 0.5，而非 0）
    obs = res["observed_vs_expected"]["corr(x−y, x−z)"]
    assert obs["observed"] == pytest.approx(obs["expected"], abs=0.08)
    assert not np.isnan(tc_scaled(a, b, c)["beta_y_vs_x"])


def test_tc_is_biased_when_errors_are_correlated():
    """记录一个已知限制：误差相关时 TC 估计会偏，且无法用三套数据自检出来。

    这是写论文时必须声明的前提风险，用测试把它固定下来。
    """
    n = 8000
    rng = np.random.default_rng(7)
    truth = rng.uniform(0.05, 0.40, n)
    shared = rng.normal(0, 0.05, n)  # a 与 b 共享的误差分量
    a = truth + shared + rng.normal(0, 0.02, n)
    b = truth + shared + rng.normal(0, 0.02, n)
    c = truth + rng.normal(0, 0.04, n)

    true_s2_a, true_s2_c = 0.05**2 + 0.02**2, 0.04**2
    tc = triple_collocation(a, b, c)
    assert tc["sigma2_x"] < 0.5 * true_s2_a, "TC 应显著低估含共享误差的那一套"
    assert tc["sigma2_z"] > 1.5 * true_s2_c, "TC 应显著高估独立误差的那一套"

    # 三套数据的自检无法发现这一点 —— 这正是限制所在
    assert "limitation" in error_correlation_check(a, b, c)


def test_tc_requires_enough_samples():
    with pytest.raises(ValueError):
        triple_collocation(np.arange(5.0), np.arange(5.0), np.arange(5.0))
