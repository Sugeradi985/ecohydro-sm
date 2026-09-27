"""土壤水反演：合成训练库 -> 神经网络反演 -> 实测残差二级校正。

为什么是"物理骨架 + NN + 残差校正"的三段式，而不是纯黑箱 ML
------------------------------------------------------------
1. 物理模型（WCM + Dubois）保证外推时的物理一致性：σ0 随含水量单调、
   随植被衰减，NN 不会学到违背散射机理的映射；
2. NN 承担非线性反演，避免逐像元迭代求解的算力负担；
3. 第二级残差校正吸收本地系统性偏差（土壤质地、定标、RTC 残差），
   这是跨区迁移能力的来源——纯黑箱模型换研究区就会崩。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .wcm import linear_to_db, wcm_forward

__all__ = [
    "SyntheticDBConfig",
    "build_synthetic_db",
    "InversionModel",
    "ResidualCorrector",
]


@dataclass
class SyntheticDBConfig:
    """合成训练库的采样范围。范围要覆盖研究区实测分布并留余量。"""

    mv_min: float = 0.02
    mv_max: float = 0.50
    mv_step: float = 0.01
    h_rms_min: float = 0.3
    h_rms_max: float = 3.0
    h_rms_step: float = 0.15
    vwc_min: float = 0.0
    vwc_max: float = 6.0
    vwc_step: float = 0.25
    theta_min: float = 30.0
    theta_max: float = 45.0
    theta_step: float = 2.5
    freq_ghz: float = 5.405
    pol: str = "VV"
    a_param: float = 0.35
    b_param: float = 0.12
    e_param: float = 1.0
    noise_db: float = 0.35  # 斑点与定标噪声（dB），防止 NN 学到零噪声的理想映射
    random_state: int = 42


def build_synthetic_db(cfg: SyntheticDBConfig | None = None):
    """用 WCM 正向模型生成合成训练库。

    返回
    ----
    X : (n, 4) 特征矩阵 [σ0_VV(dB), σ0_VH(dB), 入射角(deg), NDVI]
    y : (n,) 体积含水量 (m3/m3)
    meta : dict，采样网格信息

    说明
    ----
    NDVI 与 VWC 的换算用线性关系（需本地实测生物量标定）；
    VH 通道由 VV 加一个与 VWC 相关的体散射项近似得到，
    这是简化处理——若研究区有极化定标数据，应替换为本区经验关系。
    """
    cfg = cfg or SyntheticDBConfig()
    rng = np.random.default_rng(cfg.random_state)

    mv = np.arange(cfg.mv_min, cfg.mv_max + 1e-9, cfg.mv_step)
    h_rms = np.arange(cfg.h_rms_min, cfg.h_rms_max + 1e-9, cfg.h_rms_step)
    vwc = np.arange(cfg.vwc_min, cfg.vwc_max + 1e-9, cfg.vwc_step)
    theta = np.arange(cfg.theta_min, cfg.theta_max + 1e-9, cfg.theta_step)

    # 全网格遍历在 (>40 × 19 × 25 × 7) 规模下过大，改为随机采样组合
    n_samples = 60_000
    mv_s = rng.choice(mv, n_samples)
    h_s = rng.choice(h_rms, n_samples)
    vwc_s = rng.choice(vwc, n_samples)
    th_s = rng.choice(theta, n_samples)

    s0_vv = wcm_forward(
        mv_s, h_s, vwc_s, th_s,
        freq_ghz=cfg.freq_ghz, pol="VV",
        a_param=cfg.a_param, b_param=cfg.b_param, e_param=cfg.e_param,
    )
    s0_vh = wcm_forward(
        mv_s, h_s, vwc_s, th_s,
        freq_ghz=cfg.freq_ghz, pol="HH",
        a_param=cfg.a_param, b_param=cfg.b_param, e_param=cfg.e_param,
    ) * (0.4 + 0.1 * vwc_s)  # 交叉极化近似：以 HH 为基并按 VWC 调整

    noise = rng.normal(0.0, cfg.noise_db, n_samples)
    s0_vv_db = linear_to_db(s0_vv) + noise
    s0_vh_db = linear_to_db(np.maximum(s0_vh, 1e-6)) + noise

    ndvi = np.clip(vwc_s / 6.0, 0.0, 1.0)  # VWC -> NDVI 的粗略映射，需标定

    X = np.column_stack([s0_vv_db, s0_vh_db, th_s, ndvi])
    y = mv_s
    meta = {
        "n_samples": int(n_samples),
        "mv_range": (float(mv_s.min()), float(mv_s.max())),
        "vwc_range": (float(vwc_s.min()), float(vwc_s.max())),
        "theta_range": (float(th_s.min()), float(th_s.max())),
        "config": cfg,
    }
    return X, y, meta


class InversionModel:
    """土壤水反演器（合成库预训练 + 可选实测微调）。

    默认用 MLP；`estimator="rf"` 可切随机森林（小样本、强非线性时更稳）。
    """

    def __init__(self, estimator: str = "mlp", random_state: int = 42, **kwargs):
        if estimator == "mlp":
            base = MLPRegressor(
                hidden_layer_sizes=kwargs.pop("hidden_layer_sizes", (64, 64, 32)),
                activation="relu",
                max_iter=kwargs.pop("max_iter", 800),
                early_stopping=True,
                random_state=random_state,
            )
        elif estimator == "rf":
            base = RandomForestRegressor(
                n_estimators=kwargs.pop("n_estimators", 300),
                min_samples_leaf=kwargs.pop("min_samples_leaf", 5),
                random_state=random_state,
                n_jobs=-1,
            )
        else:
            raise ValueError(f"estimator 必须是 'mlp' 或 'rf'，收到 {estimator!r}")

        self.pipeline = Pipeline([("scaler", StandardScaler()), ("model", base)])
        self.feature_names = ["s0_vv_db", "s0_vh_db", "theta_deg", "ndvi"]

    def fit(self, X, y):
        self.pipeline.fit(np.asarray(X), np.asarray(y))
        return self

    def predict(self, X) -> np.ndarray:
        return np.asarray(self.pipeline.predict(np.asarray(X)), dtype=float)


class ResidualCorrector:
    """第二级残差校正：用实测站点残差拟合一个轻量校正器。

    输入特征应为 [预测 SM, NDVI, 砂含量%, 黏含量%] 等本地协变量。
    站点数很少（<20）时，退化为全局 bias 校正，避免过拟合。
    """

    def __init__(self, min_samples_for_model: int = 20, random_state: int = 42):
        self.min_samples_for_model = min_samples_for_model
        self.random_state = random_state
        self._model: Pipeline | None = None
        self._bias: float = 0.0

    def fit(self, X_res, residual):
        """residual = 实测 SM − 模型预测 SM。"""
        X_res = np.asarray(X_res, dtype=float)
        residual = np.asarray(residual, dtype=float)
        self._bias = float(np.mean(residual))
        if len(residual) >= self.min_samples_for_model:
            self._model = Pipeline([
                ("scaler", StandardScaler()),
                ("model", RandomForestRegressor(
                    n_estimators=200, min_samples_leaf=3,
                    random_state=self.random_state, n_jobs=-1)),
            ])
            self._model.fit(X_res, residual)
        else:
            self._model = None
        return self

    def apply(self, X_res, prediction) -> np.ndarray:
        prediction = np.asarray(prediction, dtype=float)
        if self._model is None:
            return prediction + self._bias
        return prediction + np.asarray(self._model.predict(np.asarray(X_res, dtype=float)), dtype=float)
