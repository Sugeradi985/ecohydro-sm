"""干涉对模式测试：适配 ASF HyP3 INSAR_GAMMA 的逐 pair 相干性产品。

用已知参数生成合成 γ_ij，验证标定 → 反演 → 积分整条链路能还原真值。
"""

import numpy as np
import pytest

from ecohydro_sm.insar.coherence import (
    PairCoherenceSMConfig,
    assign_pair_signs,
    calibrate_decay_coefficients_pairwise,
    integrate_pair_deltas,
    invert_delta_sm_from_pairs,
)


def _synthetic_pairs(seed=7):
    rng = np.random.default_rng(seed)
    n_dates, n_pixels = 20, 4
    t = np.arange(n_dates)[:, None]
    # 每个像元一条不同相位的季节曲线，保证锚点值各不相同
    phase = np.arange(n_pixels)[None, :] * 0.4
    sm = 0.12 + 0.09 * np.sin(2 * np.pi * t / 11 + phase)
    sm = np.clip(sm, 0.02, 0.45)

    # 连续 12 天对 + 若干 24 天对，构成连通网络
    pairs = [(i, i + 1) for i in range(n_dates - 1)]
    pairs += [(i, i + 2) for i in range(0, n_dates - 2, 4)]
    ref = np.array([p[0] for p in pairs])
    sec = np.array([p[1] for p in pairs])

    dsm = sm[sec] - sm[ref]                       # 带符号
    dv = 0.02 + 0.04 * rng.random((len(pairs), n_pixels))
    dt = (sec - ref) * 12.0

    c, k_sm, k_v, k_t = -0.02, 2.5, 1.5, 0.001
    log_gamma = c - k_sm * np.abs(dsm) - k_v * np.abs(dv) - k_t * dt[:, None]
    gamma = np.exp(log_gamma)
    return dict(sm=sm, gamma=gamma, dsm=dsm, dv=dv, dt=dt, ref=ref, sec=sec,
                n_dates=n_dates, truth=(c, k_sm, k_v, k_t))


def test_calibration_recovers_known_coefficients():
    d = _synthetic_pairs()
    coef = calibrate_decay_coefficients_pairwise(d["gamma"], d["dv"], d["dsm"], d["dt"])
    c, k_sm, k_v, k_t = d["truth"]
    assert coef["log_offset"] == pytest.approx(c, abs=0.05)
    assert coef["k_sm"] == pytest.approx(k_sm, abs=0.2)
    assert coef["k_v"] == pytest.approx(k_v, abs=0.3)
    assert coef["k_t"] == pytest.approx(k_t, abs=0.0008)
    assert coef["n_samples"] == d["gamma"].size
    assert coef["rmse"] < 0.01


def test_invert_and_integrate_recover_sm():
    d = _synthetic_pairs()
    coef = calibrate_decay_coefficients_pairwise(d["gamma"], d["dv"], d["dsm"], d["dt"])
    mag = invert_delta_sm_from_pairs(d["gamma"], d["dv"], d["dt"], coef)
    assert np.all(mag >= 0)
    assert np.mean(np.abs(mag - np.abs(d["dsm"]))) < 0.01

    signs = assign_pair_signs(d["dsm"])
    sm_hat = integrate_pair_deltas(
        mag, d["ref"], d["sec"], d["n_dates"], signs=signs,
        anchor_date=0, anchor_sm=d["sm"][0, 0],
    )
    # anchor_sm 是标量，这里逐像元锚点不同，允许整体有常数差；比较去均值后的形状
    a = d["sm"] - d["sm"].mean(axis=0, keepdims=True)
    b = sm_hat - sm_hat.mean(axis=0, keepdims=True)
    assert np.mean(np.abs(a - b)) < 0.02


def test_missing_signs_produces_monotonic_increase():
    """不提供符号时序列必然单调上升——这是必须避免的失败模式，固定在此。"""
    d = _synthetic_pairs()
    coef = calibrate_decay_coefficients_pairwise(d["gamma"], d["dv"], d["dsm"], d["dt"])
    mag = invert_delta_sm_from_pairs(d["gamma"], d["dv"], d["dt"], coef)
    sm_hat = integrate_pair_deltas(mag, d["ref"], d["sec"], d["n_dates"])
    diff = np.diff(sm_hat[:, 0])
    assert np.mean(diff) > 0
    with_signs = integrate_pair_deltas(
        mag, d["ref"], d["sec"], d["n_dates"], signs=assign_pair_signs(d["dsm"])
    )
    assert np.std(np.diff(with_signs[:, 0])) > 0  # 有符号时才会出现升降


def test_scalar_and_1d_inputs():
    d = _synthetic_pairs()
    coef = calibrate_decay_coefficients_pairwise(
        d["gamma"][:, 0], d["dv"][:, 0], d["dsm"][:, 0], d["dt"]
    )
    assert coef["k_sm"] > 0
    mag = invert_delta_sm_from_pairs(d["gamma"][:, 0], d["dv"][:, 0], d["dt"], coef)
    assert mag.shape == (len(d["ref"]),)
    out = integrate_pair_deltas(mag, d["ref"], d["sec"], d["n_dates"],
                                signs=assign_pair_signs(d["dsm"][:, 0]))
    assert out.shape == (d["n_dates"],)


def test_length_mismatch_raises():
    d = _synthetic_pairs()
    with pytest.raises(ValueError):
        integrate_pair_deltas(np.zeros(3), d["ref"], d["sec"], d["n_dates"])


def test_insufficient_samples_raises():
    with pytest.raises(ValueError):
        calibrate_decay_coefficients_pairwise(
            np.full((2, 2), 0.8), np.zeros((2, 2)), np.zeros((2, 2)), np.array([12, 12])
        )


def _constant_dt_data():
    """Δt 全部 12 天：k_t·Δt 是常数，被截距吸收。"""
    d = _synthetic_pairs()
    c, k_sm, k_v, k_t = d["truth"]
    dt_const = np.full_like(d["dt"], 12.0)
    c_eff = c - k_t * 12.0
    gamma = np.exp(c_eff - k_sm * np.abs(d["dsm"]) - k_v * np.abs(d["dv"]))
    return d, dt_const, c_eff, gamma, (k_sm, k_v)


def test_constant_dt_without_fix_is_degenerate():
    """不固定 k_t 时只有共线组合 c − k_t·12 可识别，单项不可识别。"""
    d, dt_const, c_eff, gamma, _ = _constant_dt_data()
    coef = calibrate_decay_coefficients_pairwise(
        gamma, d["dv"], d["dsm"], dt_const
    )
    assert (coef["log_offset"] - coef["k_t"] * 12.0) == pytest.approx(c_eff, abs=0.05)


def test_constant_dt_fix_k_t_zero_recovers_other_coefficients():
    """固定 k_t=0 后，c / k_sm / k_v 可正常恢复，反演链路仍可用。"""
    d, dt_const, c_eff, gamma, (k_sm, k_v) = _constant_dt_data()
    cfg = PairCoherenceSMConfig(fix_k_t_zero=True)
    coef = calibrate_decay_coefficients_pairwise(
        gamma, d["dv"], d["dsm"], dt_const, cfg=cfg
    )
    assert coef["k_t"] == 0.0
    assert coef["log_offset"] == pytest.approx(c_eff, abs=0.05)
    assert coef["k_sm"] == pytest.approx(k_sm, abs=0.2)
    assert coef["k_v"] == pytest.approx(k_v, abs=0.3)
    mag = invert_delta_sm_from_pairs(gamma, d["dv"], dt_const, coef, cfg=cfg)
    assert np.mean(np.abs(mag - np.abs(d["dsm"]))) < 0.01
