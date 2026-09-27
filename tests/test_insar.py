"""InSAR 支路与水量平衡测试：相干性去混、形变分解、储变量、闭合。"""

import numpy as np
import pytest

from ecohydro_sm.fusion.water_balance import (
    WaterBalanceConfig,
    allocate_to_root_zone,
    closure_quality_flag,
    closure_residual,
    vadose_storage_change,
)
from ecohydro_sm.insar.coherence import (
    calibrate_decay_coefficients,
    closure_phase_residual,
    invert_sm_from_coherence,
    structural_change_index,
)
from ecohydro_sm.insar.deformation import (
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


def test_structural_change_index_nonneg():
    fvc_t = np.abs(np.sin(np.linspace(0, np.pi, 40)))[:, None] * np.ones((1, 10))
    dv = structural_change_index(fvc_t)
    assert np.all(dv >= 0)
    assert np.all(np.isfinite(dv))


def test_coherence_deconvolution_improves_retrieval():
    rng = np.random.default_rng(7)
    n_t, n_px = 80, 500
    dsm = np.clip(rng.normal(0.12, 0.05, (n_t, n_px)), 0, 0.30)
    fvc = np.clip(np.outer(np.sin(np.linspace(0, np.pi, n_t)), np.ones(n_px)) * 0.6, 0, 1)
    dv = structural_change_index(fvc)
    gamma = np.clip(np.exp(-2.5 * dsm - 1.8 * dv) * (1 + rng.normal(0, 0.03, (n_t, n_px))), 0.05, 1)

    cal = calibrate_decay_coefficients(gamma, dv, dsm)
    assert cal["k_sm"] == pytest.approx(2.5, rel=0.25), f"k_sm 标定偏差过大：{cal['k_sm']}"
    assert cal["k_v"] == pytest.approx(1.8, rel=0.25), f"k_v 标定偏差过大：{cal['k_v']}"
    sm_full = invert_sm_from_coherence(gamma, dv, cal["k_sm"], cal["k_v"],
                                       log_offset=cal["log_offset"])
    sm_naive = invert_sm_from_coherence(gamma, np.zeros_like(dv), cal["k_sm"], 0.0,
                                        log_offset=cal["log_offset"])
    truth = 0.02 + dsm
    # 结构去混必须同时把偏差压下来
    assert abs(np.nanmean(sm_full - truth)) < abs(np.nanmean(sm_naive - truth))

    err_full = np.sqrt(np.nanmean((sm_full - truth) ** 2))
    err_naive = np.sqrt(np.nanmean((sm_naive - truth) ** 2))
    assert err_full < err_naive, "加入结构去混项后误差必须下降"


def test_closure_phase_residual_finite():
    trip = np.random.default_rng(8).normal(0, 0.2, (40, 3, 25))
    r = closure_phase_residual(trip)
    assert r.shape == (25,)
    assert np.all(np.isfinite(r)) and np.all(r >= 0)
    with pytest.raises(ValueError):
        closure_phase_residual(np.zeros((10, 2, 5)))


def test_decomposition_separates_trend_and_seasonal():
    rng = np.random.default_rng(9)
    n_t, n_px = 120, 50
    t = np.arange(n_t) * 6.0
    seasonal = -12.0 * np.sin(2 * np.pi * t / 365.25)[:, None] * np.ones((1, n_px))
    trend = (-6.0 * (t / 365.25))[:, None] * np.ones((1, n_px))
    disp = seasonal + trend + rng.normal(0, 1.5, (n_t, n_px))

    d = fit_decomposition(t, disp)
    assert np.allclose(d["seasonal"] + d["trend"] + d["intercept"] + d["residual"], disp, atol=1e-6)
    assert np.nanmean(d["amplitude_mm"]) == pytest.approx(12.0, abs=1.5)
    snr = seasonal_amplitude_snr(d)
    assert np.nanmedian(snr) > 3.0, "本合成情形下 SNR 应远超预检判据"


def test_lag_estimation():
    t = np.arange(120) * 6.0
    head = -np.sin(2 * np.pi * t / 365.25)
    disp = -np.sin(2 * np.pi * (t - 30) / 365.25)[:, None] * np.ones((1, 5))
    lag = estimate_lag_days(disp, head, t)
    assert abs(lag - 30) <= 12, f"滞后估计偏差过大：{lag}"


def test_ske_calibration():
    t = np.arange(120) * 6.0
    head_m = -3.0 * np.sin(2 * np.pi * t / 365.25)
    disp_mm = -15.0 * np.sin(2 * np.pi * t / 365.25)
    ske = calibrate_ske(disp_mm, head_m)
    assert ske == pytest.approx(5.0e-3, abs=1e-3)


def test_storage_paths_consistent_sign():
    ds_head = storage_change_from_head(np.array([1.0, -1.0]), sy=0.10)
    assert ds_head[0] > 0 > ds_head[1], "水位上升应对应储变量增加"
    ds_def = storage_change_from_deformation(np.array([10.0, -10.0]), ske=5e-3, sy=0.10)
    assert ds_def[0] > 0 > ds_def[1]
    with pytest.raises(ValueError):
        storage_change_from_deformation(np.array([1.0]), ske=0.0, sy=0.1)


def test_sy_sensitivity_gives_interval():
    disp = np.random.default_rng(10).normal(0, 8, (60, 30))
    span = sy_sensitivity_span(disp, 5e-3, StorageConfig(n_sensitivity_samples=100))
    assert np.all(span["p05_mm"] <= span["median_mm"])
    assert np.all(span["median_mm"] <= span["p95_mm"])
    assert np.nanstd(span["p95_mm"] - span["p05_mm"]) > 0, "必须给出非退化的区间"


def test_water_balance_closes():
    rng = np.random.default_rng(11)
    p = rng.gamma(1.2, 8.0, (50, 40))
    et = rng.gamma(2.0, 4.0, (50, 40))
    dgw = rng.normal(-2, 3, (50, 40))
    dv = vadose_storage_change(p, et, dgw)
    assert np.allclose(closure_residual(p, et, dgw, dv), 0.0, atol=1e-9)


def test_root_zone_allocation_bounds():
    rng = np.random.default_rng(12)
    dv = rng.normal(0, 20, (50, 40))
    sm_surf = np.clip(rng.normal(0.12, 0.04, (50, 40)), 0.02, 0.35)
    sm_rz = allocate_to_root_zone(dv, sm_surf, WaterBalanceConfig())
    assert sm_rz.shape == sm_surf.shape
    assert np.all((sm_rz >= 0) & (sm_rz <= 0.6))


def test_closure_quality_flag():
    p = np.array([50.0, 50.0])
    assert closure_quality_flag(p, np.array([1.0, 7.0]))[0]
    assert not closure_quality_flag(p, np.array([20.0, 1.0]))[0]
