# Yanco 试验区：从已下载数据到第一批结果的完整处理链

> 适用前提：已收集公开土壤湿度观测、已用 ASF HyP3 处理 S-1 干涉数据、
> 已用 GEE 处理同范围 S-1 GRD 与 S-2 并下载到本地。
> 配套：`configs/yanco.yaml`、`scripts/00_prepare_yanco.py`、
> `experiment_runbook.md`（通用判据）、`station_requirements.md`。

---

## 0. 先确认三件事（30 分钟，别跳过）

### 0.1 地理范围：Yanco 在澳大利亚**东南部**，不是西南部

Yanco 位于新南威尔士州 Murrumbidgee 流域西部平原，约 **145.8–146.5°E / 34.7–35.2°S**，
UTM 带 **55S（EPSG:32755）**，地势平坦（117–150 m），年降水约 394 mm。
如果 GEE 的 bbox 落在西澳（约 115–118°E），下载的数据整批作废。

**动作**：打开任意一景 GeoTIFF，核对 `bounds` 是否落在上述经纬度范围内。

### 0.2 「土壤湿度观测产品」是站点还是栅格？

这是决定后续能不能训练的分叉：

| 你手上的数据 | 能做什么 | 不能做什么 |
|---|---|---|
| **OzNet / Yanco 站点观测**（oznet.org.au，20 min，0–5 cm；Y1–Y13 剖面 0–90 cm，YA 13 站 / YB 11 站表层簇） | 标签、残差校正、完整精度验证 | — |
| SMAP / SMOS / ESA CCI / ERA5-Land **栅格产品** | 三重配置、相对比较、基准 | **不能当标签训练**，不能报「与实测的 RMSE」 |

**动作**：确认文件里有没有 `station_id` 和 `depth_cm` 字段。没有 → 走无真值路线。

### 0.3 时区：会整体错开一天

S-1 文件名日期是 **UTC**；OzNet 站点时间戳多为澳洲东部时（AEST UTC+10 / AEDT UTC+11）。
**降轨约当地清晨 6–7 点 = 前一日 19–21 UTC**，所以按文件名日期直接匹配站点当地日期会整体错位。

**动作**：把站点时间统一换算为 UTC（用固定偏移，不要用系统本地时区），
再取 S-1 **过境时刻**的值，而不是日平均。日平均会平滑掉降水/灌溉响应。

---

## 1. 数据清点（脚本 00，只读）

```bash
python scripts/00_prepare_yanco.py --root . --config configs/yanco.yaml
python scripts/00_prepare_yanco.py --selftest   # 不依赖真实数据，先验证工具本身
```

产出 `outputs/yanco/` 下：

| 文件 | 内容 |
|---|---|
| `manifest_grd.csv` | 每景 GRD 的日期、CRS、分辨率、边界、数值范围 |
| `manifest_s2.csv` | S-2 同上 |
| `manifest_insar_pairs.csv` | 每个干涉对的参考/副影像日期、Δt、γ 统计 |
| `match_table.csv` | S-1 日期 ↔ 最近有效 S-2、间隔天数、覆盖该日期的干涉对 |
| `inventory_report.md` | PASS / WARN / FAIL / NOTE 检查清单 |

**必须清零所有 FAIL 项**再往下走，尤其是网格一致性。

---

## 2. 辐射口径统一（最容易埋雷的一步）

| 来源 | 给的是什么 |
|---|---|
| GEE `COPERNICUS/S1_GRD` | **σ⁰（dB）**，含 `angle` 波段（椭球入射角）；已做轨道、热噪声去除、辐射定标、距离-多普勒地形校正；**未做斑点滤波、未做角度辐射坡度校正** |
| HyP3 `INSAR_GAMMA` | `amplitude` = 参考景的 **σ⁰ 线性功率**（仅几何地形校正，无辐射地形校正）；`corr` = 0–1 相干性；`unw_phase`、`water_mask`；像元 **40 或 80 m** |
| HyP3 `RTC_GAMMA` | **γ⁰**（`radiometry=gamma0`），与上面两种都不同口径 |

**关键结论**：仓库 `retrieval/wcm.py::dubois_sigma0` 预测的正是 **σ⁰ 线性功率**，
所以 GEE GRD 转线性后**恰好与合成库一致**，比 RTC γ⁰ 更顺。
若混用 RTC γ⁰，需 γ⁰ ≈ σ⁰ / cosθ 转换（Yanco 平坦，误差小但口径必须统一）。

**动作清单**：

1. `10 ** (dB / 10)` 把 GRD 转成线性功率；
2. 用 `preprocess/speckle.py::lee_sigma` 做斑点滤波（L1 支路必须做）；
3. **相干性支路不要复用被滤波的幅度**——两条支路保持独立，否则交叉验证失去意义；
4. 把 10 m GRD/S-2 与 40/80 m 相干性重投影到**同一张参考网格**（EPSG:32755），
   判据：同名点残差 **< 0.5 像元**；
5. 标称分辨率 = 实际聚合后分辨率，不可虚标。

---

## 3. S-2 指数与掩膜（脚本 01）

- 计算 NDVI / NDMI / FVC / VWC；SCL 云掩膜 + 腐蚀缓冲
- Yanco 基本无冻融 → 可跳过 `freeze_thaw_mask`（仍用土壤温度复核一次）
- 入射角掩膜 30–45°；过境前后 1 天降水 > 1 mm 剔除
- **VWC 域掩膜**：YA 灌溉作物（水稻/大麦/玉米）VWC 常 > 2 kg/m²，
  L2 在这些像元**打标降级，不承诺精度**
- **不要做月度最大值合成**（引入物候偏差）；
  用「最近有效光学观测 + 距 SAR 过境天数」作为独立特征，让模型学时间错配
- `match_table.csv` 里的 `gap_days` 就是这一列的来源

---

## 4. 站点表准备

- 20 min → 重采样到 S-1 过境时刻（UTC），只用 **0–5 cm** 层
- 站点 → 像元抽表：用 S-2 NDVI 的**像元内变异系数**筛代表性，异质像元剔除
- **YA / YB 必须分开**，按 `landcover × irrigation` 分层报告，不可合并成一个精度数字
- **灌溉记录**：YA 是灌溉作物。没有它，"降水后 1–3 天上升"的事件响应检验会把灌溉误判成降水
- 站点数：`residual_correction.min_samples_for_model = 20`；
  低于 20 时 `ResidualCorrector` 自动退化为全局 bias，不做协变量建模

---

## 5. L1 基线（脚本 02 → 03 → 04）

```
02_build_synthetic_db.py   Topp → Dubois σ0 → WCM 正演，生成合成训练库
03_train_inversion.py      MLP / RF 训练 + 实测残差二级校正
04_retrieve_sm.py          生产 + 聚合
```

- WCM 的 **A / B / E 必须本地标定**；Yanco 已有文献里的 VWC–NDVI 关系只能当**初值**
- 合成库自洽 ubRMSE < 0.05 m³/m³ 才继续
- 升/降轨分别反演再融合，同期差 **< 0.05 m³/m³**；超出说明入射角或定标有问题
- 验收：ubRMSE ≤ 0.04、|bias| ≤ 0.02、anomaly R ≥ 0.6、干湿动态范围 ≥ 3 dB

---

## 6. L2：HyP3 干涉对与仓库模型的**结构性不匹配**（重点）

仓库原有的 `invert_sm_from_coherence` 要的是 **逐时相 γ(t)**（每个日期相对干参考的相干性），
而 HyP3 `INSAR_GAMMA` 给的是 **逐干涉对 γ_ij**。
**把 pair 堆成 `(n_time, ...)` 直接喂进去，得到的是语义错误的"时序"，反演结果没有物理意义。**

两条路：

### 路线 A（推荐先用）：pairwise 模式，不重跑数据

模型改为对 pair 建模（已在 `insar/coherence.py` 实现并有 6 项测试）：

```
ln γ_ij = c − k_sm·|ΔSM_ij| − k_v·|ΔV_ij| − k_t·Δt_ij
```

| 函数 | 作用 |
|---|---|
| `calibrate_decay_coefficients_pairwise` | 同时标定 c、k_sm、k_v、**k_t** |
| `invert_delta_sm_from_pairs` | 反演每个 pair 的 **\|ΔSM\|** |
| `assign_pair_signs` | 由站点/L1/ERA5-Land 的 ΔSM 给符号 |
| `integrate_pair_deltas` | 最小二乘从干锚点积分成逐日期 SM 序列 |

```python
from ecohydro_sm.insar.coherence import (
    calibrate_decay_coefficients_pairwise, invert_delta_sm_from_pairs,
    assign_pair_signs, integrate_pair_deltas,
)

coef = calibrate_decay_coefficients_pairwise(gamma_pairs, dv_pairs, dsm_ref_pairs, dt_days)
mag = invert_delta_sm_from_pairs(gamma_pairs, dv_pairs, dt_days, coef)
sm = integrate_pair_deltas(mag, ref_idx, sec_idx, n_dates,
                           signs=assign_pair_signs(dsm_ref_pairs),
                           anchor_date=dry_idx, anchor_sm=0.02)
```

**三个不能省的点**：

1. **k_t 项**（时间去相干）：12 天对与 24 天对的基线不同，不建模会把时间去相干吸进 k_sm；
2. **符号必须由外部参考给**：反演只给幅度。不提供符号会得到单调上升的 SM 序列——
   `test_missing_signs_produces_monotonic_increase` 把这个失败模式固定下来了；
3. **闭合相位质控本轮跳过**：单对产品不含闭合相位（需三景三元组），
   论文里要明确声明这一限制。

### 路线 B：重跑 star 网络

以固定「干参考景」重新提交 HyP3，得到 per-date γ，就回到仓库原有的
`calibrate_decay_coefficients` / `invert_sm_from_coherence`。
但长基线会引入时间去相干，仍需在模型里加 Δt 项，且**重跑成本高**。

### L2 判据

拟合 R² ≥ 0.5、k_sm/k_v/k_t 均为正；干旱/裸土区平均 γ > 0.4；
L1 与 L2 差 < 0.06 m³/m³（两者观测机制独立，这个交叉很有说服力）。
比较「是否加结构去混」时看 **RMSE 与 bias**，不要看 ubRMSE（它先剔除了偏差）。

---

## 7. 验证（脚本 08）

- **空间块 CV**（块 ≥ 10 km）+ 随机 CV **两个都报**，R² 差 > 0.15 说明泄漏严重
- leave-one-station-out
- 按 VWC / 土地覆被 / 灌溉制度**分层报告**
- 指标分离：bias / ubRMSE / cRMSE / anomaly R
- 三重配置：L1 / ERA5-Land / SMAP。注意 SMAP 与 ERA5-Land 都依赖陆面模式，
  误差可能相关，TC 有已知偏倚（代码里已显式记录）
- **Yanco 的最大优势**：它是 SMAP/SMOS 定标验证的经典区，有大量已发表精度基准
  （SMAPEx、PRISM、OzNet 系列），可以直接做国际可比的公平对比

---

## 8. 本轮不要做 L3 / L4

`configs/yanco.yaml` 已把两者关闭：

- HyP3 单对干涉产品**不构成时序形变**，无法做谐波+趋势分解；
- Yanco 不是地面沉降区，也没有配套地下水位井；
- 缺灌溉与抽水记录时，水量平衡闭合残差会系统性偏大。

**把目标收缩到 L1 + L2，把验证纪律做扎实**，这一轮的产出才干净。

---

## 9. Yanco 在整体课题里的定位

| | Yanco | 民勤盆地 |
|---|---|---|
| 角色 | **方法学验证区** | **科学问题区** |
| 优势 | 数据全公开、20 min 长期、密集站网、国际基准可比 | 超采—治理反向自然实验，Han 2020 vs Tian 2026 争议落到 ΔS 缺口 |
| 产出 | L1/L2 方法学论文、可复现的精度数字 | L4/L3 生态水文机理、地下水闭合 |

Yanco 比黑河更容易被国际审稿人接受，正好接替此前「双研究区」里方法学区的角色。
但**它不承载科学故事**：碳-水权衡与地下水闭合仍要回到民勤。

---

## 10. 执行顺序速查

```
[ ] 00_prepare_yanco.py --selftest           工具自检
[ ] 00_prepare_yanco.py --root .            清点，清零所有 FAIL
[ ] 核对 bbox / 时区 / 真值类型 / 辐射口径
[ ] 统一网格（EPSG:32755），同名点残差 < 0.5 像元
[ ] GRD dB → 线性 + 斑点滤波（相干性支路不复用）
[ ] 站点表：UTC 过境时刻、0-5 cm、YA/YB 分层、代表性筛选
[ ] 01_prepare_s2.py                        指数与掩膜
[ ] 02 → 03 → 04                            L1 基线
[ ] L2 pairwise：标定 → 反演 → 符号 → 积分
[ ] 08_validate.py                          空间块 CV + 分层 + 三重配置
[ ] 与 Yanco 已发表基准对比，出第一张图
```
