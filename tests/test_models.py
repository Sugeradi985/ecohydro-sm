"""正向散射模型、光谱指数与变化检测的数值测试（不需要任何影像数据）。"""

import numpy as np
import pytest

from ecohydro_sm.preprocess.indices import fvc, ndmi, ndvi, vwc_from_ndvi
from ecohydro_sm.retrieval.change_detection import (
    ChangeDetectionConfig,
    change_detection_sm,
    dry_wet_references,
    dynamic_range_flag,
)
from ecohydro_sm.retrieval.wcm import (
    db_to_linear,
    dubois_sigma0,
    linear_to_db,
    topp_dielectric,
    topp_inverse,
    water_cloud_transmissivity,
    wcm_forward,
)


def test_topp_roundtrip():
    mv = np.linspace(0.05, 0.40, 10)
    eps = topp_dielectric(mv)
    assert np.all(np.diff(eps) > 0), "介电常数应随含水量单调上升"
    assert np.allclose(topp_inverse(eps), mv, atol=0.02), "Topp 逆变换不闭合"


def test_dubois_monotonic_in_moisture():
    mv = np.linspace(0.05, 0.35, 12)
    s0 = dubois_sigma0(mv, h_rms=1.0, theta_deg=35.0)
    assert np.all(np.diff(s0) > 0), "σ0 应随含水量单调上升"
    assert np.all(s0 > 0) and np.all(np.isfinite(s0))


def test_dubois_rejects_invalid_input():
    with pytest.raises(ValueError):
        dubois_sigma0(0.2, h_rms=0.0, theta_deg=35.0)
    with pytest.raises(ValueError):
        dubois_sigma0(0.2, h_rms=1.0, theta_deg=35.0, pol="VH")


def test_wcm_vegetation_attenuates_soil():
    soil_only = wcm_forward(0.25, 1.0, 0.0, 35.0)
    with_veg = wcm_forward(0.25, 1.0, 3.0, 35.0)
    # 植被既衰减土壤贡献、又叠加自身散射；总后向散射应高于纯土壤
    assert with_veg > soil_only
    tau = water_cloud_transmissivity(3.0, 35.0)
    assert 0 < tau < 1, "透过率应在 (0, 1)"
    assert tau < 0.5, "VWC=3 时土壤贡献应被显著衰减"


def test_db_roundtrip():
    x = np.array([0.01, 0.1, 1.0])
    assert np.allclose(db_to_linear(linear_to_db(x)), x)
    with pytest.raises(ValueError):
        linear_to_db(np.array([-1.0]))


def test_indices_bounds():
    red = np.array([0.05, 0.10, 0.30])
    nir = np.array([0.40, 0.20, 0.10])
    v = ndvi(red, nir)
    assert np.all(np.abs(v) <= 1.0)
    m = ndmi(np.array([0.4, 0.3, 0.2]), np.array([0.2, 0.25, 0.3]))
    assert np.all(np.abs(m) <= 1.0)
    f = fvc(v)
    assert np.all((f >= 0) & (f <= 1))
    assert np.all(vwc_from_ndvi(v) >= 0)


def test_change_detection_scaling():
    rng = np.random.default_rng(0)
    sm_true = np.linspace(0.02, 0.45, 50)
    s0 = -20.0 + 8.0 * (sm_true - 0.02) / 0.43 + rng.normal(0, 0.1, 50)
    # 形状必须为 (n_time, n_pixel)：干湿参考沿时间轴取
    s0 = s0[:, None]
    dry, wet = dry_wet_references(s0, ChangeDetectionConfig())
    sm = change_detection_sm(s0, dry, wet)
    assert np.all((sm >= 0) & (sm <= 0.6))
    assert np.corrcoef(sm.ravel(), sm_true)[0, 1] > 0.95


def test_dynamic_range_flag():
    assert dynamic_range_flag(np.array([-22.0]), np.array([-14.0]))[0]
    assert not dynamic_range_flag(np.array([-16.0]), np.array([-15.0]))[0]
