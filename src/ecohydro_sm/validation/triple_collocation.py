"""三重配置（Triple Collocation）。

用途
----
站点稀疏区（如西北内陆）没有足够的地面真值，TC 能在**不需要真值**的条件下
估计三套独立产品各自的随机误差方差。这是本方案在缺站点时的保底验证手段。

前提假设（必须检验，不可默认成立）
----------------------------------
1. 三个数据集的误差相互独立；
2. 误差与真值正交（无系统偏差耦合）；
3. 误差零均值、方差稳定。

SMAP 与 ERA5-Land 都依赖陆面模式，独立性可能不足——
使用前务必做误差相关性诊断，否则 TC 结果不可信。
"""

from __future__ import annotations

import numpy as np

__all__ = ["triple_collocation", "tc_scaled", "error_correlation_check"]


def triple_collocation(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """标准三重配置：估计三套数据各自的误差方差。

    参数
    ----
    x, y, z : 三套独立的土壤水估计（同形，单位一致）

    返回
    ----
    dict，含 sigma2_x / y / z（误差方差）与 rmse_x / y / z
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    if len(x) < 10:
        raise ValueError("三重配置至少需要 10 个有效三元组。")

    # 以 x 为基准做相对定标，消除尺度差异
    s2_x = np.var(x, ddof=1) - np.cov(x, y)[0, 1] * np.cov(x, z)[0, 1] / np.cov(y, z)[0, 1]
    s2_y = np.var(y, ddof=1) - np.cov(y, x)[0, 1] * np.cov(y, z)[0, 1] / np.cov(x, z)[0, 1]
    s2_z = np.var(z, ddof=1) - np.cov(z, x)[0, 1] * np.cov(z, y)[0, 1] / np.cov(x, y)[0, 1]

    def _safe(v):
        return float(v) if v > 0 else float("nan")

    return {
        "n": int(len(x)),
        "sigma2_x": _safe(s2_x), "sigma2_y": _safe(s2_y), "sigma2_z": _safe(s2_z),
        "rmse_x": float(np.sqrt(max(s2_x, 0))),
        "rmse_y": float(np.sqrt(max(s2_y, 0))),
        "rmse_z": float(np.sqrt(max(s2_z, 0))),
    }


def tc_scaled(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """带尺度因子的三重配置（Stoffelen 形式），同时给出相对真值的缩放系数。

    当三套数据的动态范围差异较大时用这个版本，比标准 TC 更稳。
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]

    c_xy = np.cov(x, y)[0, 1]
    c_xz = np.cov(x, z)[0, 1]
    c_yz = np.cov(y, z)[0, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_yx = c_xy / np.var(x, ddof=1)
        beta_zx = c_xz / np.var(x, ddof=1)
    s2_x = np.var(x, ddof=1) - c_xy * c_xz / c_yz if c_yz != 0 else np.nan
    return {
        "n": int(len(x)),
        "beta_y_vs_x": float(beta_yx),
        "beta_z_vs_x": float(beta_zx),
        "sigma2_x": float(s2_x) if np.isfinite(s2_x) and s2_x > 0 else float("nan"),
        "rmse_x": float(np.sqrt(max(s2_x, 0))) if np.isfinite(s2_x) else float("nan"),
    }


def error_correlation_check(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, threshold: float = 0.2
) -> dict:
    """TC 前提诊断：误差独立性是否成立。

    原理
    ----
    若三套数据的误差相互独立，则差分序列的相关系数有确定期望：

        E[corr(x−y, x−z)] = σ²_x / sqrt((σ²_x + σ²_y)(σ²_x + σ²_z))

    **注意 1**：即使误差完全独立，该相关系数也不为 0（误差方差相近时约 0.5）。
    因此不能简单地"相关系数高 = 不独立"，必须拿观测值跟独立假设下的期望值比。

    **注意 2（重要限制）**：仅有三套数据时，误差相关性在数学上是**不可识别**的
    ——TC 的期望值本身就由 TC 估计的方差算出，因此本检查是**必要条件而非充分条件**。
    实测表明：当 x 与 y 共享误差分量时，TC 会低估 σ²_x、高估 σ²_z，
    但本检查的偏差量仍然很小，**无法可靠检出**。

    若独立性存疑（例如 SMAP 与 ERA5-Land 都依赖陆面模式），必须：
    1) 引入第四套独立数据；或 2) 用地面站点做绝对定标；
    本函数的 `independence_ok` 只表示"未发现内部矛盾"，不代表已验证独立。

    参数
    ----
    threshold : 观测与期望相关系数的允许偏差
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]

    est = triple_collocation(x, y, z)
    s2x, s2y, s2z = est["sigma2_x"], est["sigma2_y"], est["sigma2_z"]

    def corr(a, b):
        return float(np.corrcoef(a, b)[0, 1]) if len(a) > 3 else float("nan")

    def expected(shared, o1, o2):
        denom = np.sqrt((shared + o1) * (shared + o2))
        return float(shared / denom) if denom > 0 else float("nan")

    checks = {
        "corr(x−y, x−z)": (corr(x - y, x - z), expected(s2x, s2y, s2z)),
        "corr(y−x, y−z)": (corr(y - x, y - z), expected(s2y, s2x, s2z)),
        "corr(z−x, z−y)": (corr(z - x, z - y), expected(s2z, s2x, s2y)),
    }
    deviations = {
        k: (obs - exp) for k, (obs, exp) in checks.items()
        if np.isfinite(obs) and np.isfinite(exp)
    }
    variances_admissible = all(
        np.isfinite(v) and v > 0 for v in (s2x, s2y, s2z)
    )
    ok = bool(variances_admissible and deviations and
              all(abs(v) < threshold for v in deviations.values()))
    return {
        "observed_vs_expected": {k: {"observed": o, "expected": e} for k, (o, e) in checks.items()},
        "deviation": deviations,
        "variances_admissible": bool(variances_admissible),
        "independence_ok": ok,
        "threshold": threshold,
        "limitation": "三套数据无法识别误差相关性；independence_ok 仅为必要条件。"
                      "独立性存疑时请引入第四套数据或地面站点。",
    }
