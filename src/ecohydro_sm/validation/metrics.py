"""精度指标。

为什么必须拆开报
----------------
只报一个 RMSE 会把三类完全不同的误差混在一起：
  bias    —— 系统性偏差（定标、模型参数），可校正；
  ubRMSE  —— 去偏后的随机误差，才是**真正的方法精度**；
  cRMSE   —— 相对气候态的误差，反映对异常的刻画能力。

anomaly R（去季节气候态后的相关）才是生态水文应用真正关心的量：
对气候态拟合得好不代表能捕捉干旱事件。
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "bias",
    "rmse",
    "ubrmse",
    "crmse",
    "pearson_r",
    "anomaly_r",
    "stratified_report",
    "full_report",
]


def _clean(pred, ref):
    p = np.asarray(pred, dtype=float).ravel()
    r = np.asarray(ref, dtype=float).ravel()
    m = np.isfinite(p) & np.isfinite(r)
    return p[m], r[m]


def bias(pred, ref) -> float:
    p, r = _clean(pred, ref)
    return float(np.mean(p - r)) if p.size else float("nan")


def rmse(pred, ref) -> float:
    p, r = _clean(pred, ref)
    return float(np.sqrt(np.mean((p - r) ** 2))) if p.size else float("nan")


def ubrmse(pred, ref) -> float:
    """去偏 RMSE：随机误差。方法精度应以此为准。"""
    p, r = _clean(pred, ref)
    if not p.size:
        return float("nan")
    return float(np.sqrt(np.mean(((p - np.mean(p)) - (r - np.mean(r))) ** 2)))


def crmse(pred, ref) -> float:
    """相对气候态的 RMSE。"""
    p, r = _clean(pred, ref)
    if not p.size:
        return float("nan")
    clim = np.mean(r)
    return float(np.sqrt(np.mean(((p - clim) - (r - clim)) ** 2)))


def pearson_r(pred, ref) -> float:
    p, r = _clean(pred, ref)
    if p.size < 3:
        return float("nan")
    return float(np.corrcoef(p, r)[0, 1])


def anomaly_r(pred, ref, climatology=None) -> float:
    """异常相关：剔除气候态后的相关系数。

    参数
    ----
    climatology : 参考序列的季节气候态（与 ref 同长）。
        缺省用 ref 的均值（样本不足时至少去均值）。
    """
    p, r = _clean(pred, ref)
    if p.size < 3:
        return float("nan")
    clim = np.mean(r) if climatology is None else np.asarray(climatology, dtype=float).ravel()
    if np.ndim(clim) == 0:
        clim = np.full(r.shape, float(clim))
    pa, ra = p - clim, r - clim
    if np.std(pa) < 1e-12 or np.std(ra) < 1e-12:
        return float("nan")
    return float(np.corrcoef(pa, ra)[0, 1])


def full_report(pred, ref, climatology=None) -> dict:
    """一次性给出全部核心指标。"""
    return {
        "n": int(np.sum(np.isfinite(np.asarray(pred, dtype=float)) & np.isfinite(np.asarray(ref, dtype=float)))),
        "bias": bias(pred, ref),
        "rmse": rmse(pred, ref),
        "ubrmse": ubrmse(pred, ref),
        "crmse": crmse(pred, ref),
        "pearson_r": pearson_r(pred, ref),
        "anomaly_r": anomaly_r(pred, ref, climatology),
    }


def stratified_report(pred, ref, stratify_by, bins, bin_labels=None, climatology=None) -> dict:
    """分层报告精度。

    **禁止只报一个全样本数值**——按 VWC 区间分层是硬性要求，
    因为 C 波段土壤水的精度退化主要发生在高植被区。

    参数
    ----
    stratify_by : 分层变量（如 VWC）
    bins : 分层的边界序列，如 [0, 0.5, 1.5, 3.0, np.inf]
    """
    pred = np.asarray(pred, dtype=float).ravel()
    ref = np.asarray(ref, dtype=float).ravel()
    s = np.asarray(stratify_by, dtype=float).ravel()
    m = np.isfinite(pred) & np.isfinite(ref) & np.isfinite(s)
    pred, ref, s = pred[m], ref[m], s[m]

    out = {}
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        sel = (s >= lo) & (s < hi)
        label = bin_labels[i] if bin_labels else f"[{lo}, {hi})"
        out[label] = full_report(pred[sel], ref[sel], climatology) if sel.sum() > 2 else {
            "n": int(sel.sum()), "bias": float("nan"), "rmse": float("nan"),
            "ubrmse": float("nan"), "crmse": float("nan"),
            "pearson_r": float("nan"), "anomaly_r": float("nan"),
        }
    return out
