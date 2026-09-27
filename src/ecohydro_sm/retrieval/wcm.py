"""正向散射模型：土壤介电混合 + 裸土后向散射 + 水云模型植被层。

模型来源
--------
Topp et al. (1980)    介电常数混合模型
Dubois et al. (1995)  裸土同极化后向散射半经验模型
Attema & Ulaby (1978) 水云模型（WCM）植被层

使用前提（务必遵守）
--------------------
Dubois 模型的适用域：k·h_rms <= 2.5、入射角 >= 30°、mv <= 35%。
WCM 的 A / B / E 是**植被类型相关经验参数**，必须本地标定，
不可跨气候区直接套用文献值。正向模型的作用是生成合成训练库的先验分布，
最终产品精度由实测站点的残差校正保证。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "topp_dielectric",
    "topp_inverse",
    "dubois_sigma0",
    "water_cloud_transmissivity",
    "wcm_forward",
    "linear_to_db",
    "db_to_linear",
]


def topp_dielectric(mv: np.ndarray | float) -> np.ndarray | float:
    """Topp 模型：体积含水量 -> 土壤介电常数实部。

    参数
    ----
    mv : 体积含水量 (m3/m3)

    返回
    ----
    介电常数实部 ε'（无量纲）
    """
    mv = np.asarray(mv, dtype=float)
    return 3.03 + 9.3 * mv + 146.0 * mv**2 - 76.7 * mv**3


def topp_inverse(eps: np.ndarray | float) -> np.ndarray | float:
    """Topp 模型的逆形式：介电常数 -> 体积含水量。"""
    eps = np.asarray(eps, dtype=float)
    return -5.3e-2 + 2.92e-2 * eps - 5.5e-4 * eps**2 + 4.3e-6 * eps**3


def dubois_sigma0(
    mv: np.ndarray | float,
    h_rms: np.ndarray | float,
    theta_deg: np.ndarray | float,
    freq_ghz: float = 5.405,
    pol: str = "VV",
) -> np.ndarray | float:
    """Dubois et al. (1995) 裸土同极化后向散射系数（线性值）。

    参数
    ----
    mv : 体积含水量 (m3/m3)
    h_rms : 地表均方根高度 (cm)
    theta_deg : 入射角 (度)
    freq_ghz : 频率 (GHz)，默认 Sentinel-1 C 波段 5.405 GHz
    pol : 'VV' 或 'HH'

    返回
    ----
    线性域 σ0（非 dB）

    说明
    ----
    公式形式经四个独立来源核对一致。注意：本模型给出的 HH/VV 比值行为
    在不同粗糙度与入射角组合下需自行验证，若用于交叉极化比相关工作，
    请先与 Dubois et al. (1995) 原文核对。
    """
    pol = pol.upper()
    if pol not in ("VV", "HH"):
        raise ValueError(f"pol 必须是 'VV' 或 'HH'，收到 {pol!r}")

    mv = np.asarray(mv, dtype=float)
    h_rms = np.asarray(h_rms, dtype=float)
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))

    eps = topp_dielectric(mv)
    lam_cm = 30.0 / freq_ghz  # 波长，cm（c = 3e10 cm/s）
    k = 2.0 * np.pi / lam_cm  # 波数，1/cm

    ks_sin = k * h_rms * np.sin(theta)
    if np.any(ks_sin <= 0):
        raise ValueError("k·h_rms·sin(theta) 必须为正，请检查 h_rms 与入射角。")

    cos_t, sin_t, tan_t = np.cos(theta), np.sin(theta), np.tan(theta)

    if pol == "VV":
        return (
            10.0 ** (-2.35)
            * (cos_t**3 / sin_t**3)
            * 10.0 ** (0.046 * eps * tan_t)
            * ks_sin**1.1
            * lam_cm**0.7
        )
    return (
        10.0 ** (-2.75)
        * (cos_t**1.5 / sin_t**5)
        * 10.0 ** (0.028 * eps * tan_t)
        * ks_sin**1.4
        * lam_cm**0.7
    )


def water_cloud_transmissivity(vwc: np.ndarray | float, theta_deg, b_param: float = 0.12) -> np.ndarray | float:
    """水云模型的植被双向透过率 τ² = exp(−2·B·VWC/cosθ)。

    参数
    ----
    vwc : 植被含水量 (kg/m2)
    b_param : 衰减系数 B，**必须本地标定**
    """
    vwc = np.asarray(vwc, dtype=float)
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))
    return np.exp(-2.0 * b_param * vwc / np.cos(theta))


def wcm_forward(
    mv,
    h_rms,
    vwc,
    theta_deg,
    freq_ghz: float = 5.405,
    pol: str = "VV",
    a_param: float = 0.35,
    b_param: float = 0.12,
    e_param: float = 1.0,
) -> np.ndarray:
    """水云模型总后向散射系数（线性域）。

        σ0_total = σ0_veg + τ² · σ0_soil
        σ0_veg   = A · VWC^E · cosθ · (1 − τ²)
        τ²       = exp(−2 B · VWC / cosθ)

    参数
    ----
    mv, h_rms, vwc, theta_deg : 见 dubois_sigma0
    a_param, b_param, e_param : WCM 经验参数 A / B / E，**必须本地标定**

    返回
    ----
    线性域总 σ0
    """
    tau2 = water_cloud_transmissivity(vwc, theta_deg, b_param)
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))
    sigma_soil = dubois_sigma0(mv, h_rms, theta_deg, freq_ghz=freq_ghz, pol=pol)
    vwc_arr = np.asarray(vwc, dtype=float)
    sigma_veg = a_param * vwc_arr**e_param * np.cos(theta) * (1.0 - tau2)
    return np.asarray(sigma_veg + tau2 * sigma_soil, dtype=float)


def linear_to_db(x) -> np.ndarray:
    """线性功率 -> dB。"""
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        raise ValueError("线性域 σ0 必须为正才能转 dB。")
    return 10.0 * np.log10(x)


def db_to_linear(x) -> np.ndarray:
    """dB -> 线性功率。"""
    return 10.0 ** (np.asarray(x, dtype=float) / 10.0)
