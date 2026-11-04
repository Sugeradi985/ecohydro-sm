# 时间支撑对齐：逐干涉对相干性 ↔ 月度光学 / GRD

> 适用：ASF HyP3 `INSAR_GAMMA` 逐干涉对产品 + GEE 导出的**月度** S-2 / S-1 GRD。
> 配套：`configs/yanco.yaml` 的 `temporal_alignment` 段、
> `src/ecohydro_sm/preprocess/temporal_alignment.py`、
> `scripts/00b_align_temporal.py`、`docs/yanco_workflow.md`。
>
> 命名说明：本文用 A/B/C 指**时间支撑对齐**路线；`yanco_workflow.md` §6 里
> 另有一组「pairwise vs star 网络」的路线 A/B，两者无关，故那边记作 T-A/T-B/T-C。

---

## 1. 三类数据不在同一根时间轴上

| 数据 | 观测量 | 时间支撑 | 物理含义 |
|---|---|---|---|
| HyP3 `corr` | γ_ij（每个干涉对一份） | **[t_i, t_j]**，通常 12 天 | **差分**：这段时间里散射体变了多少 |
| GEE S-2 月值 | NDVI / NDMI / FVC / VWC | ≈30 天合成 | **状态**：本月平均或最大冠层状态 |
| GEE GRD 月值 | σ⁰（dB） | ≈30 天合成 | **状态**：本月平均后向散射 |
| OzNet 站点 | 体积含水量 | 20 min | 近瞬时，可当任意时刻的**真值锚** |

相干性衰减模型（仓库 `insar/coherence.py` 的 pairwise 形式）是：

```
ln γ_ij = c − k_sm·|ΔSM_ij| − k_v·|ΔV_ij| − k_t·Δt_ij
```

其中 **ΔV_ij = V(t_j) − V(t_i)** 是**两个端点日期上的差分**。月值给不出这两个端点。

---

## 2. 直接按「干涉对所在月份」取月值，会犯三个错

**① 量纲错**。把 V_month（状态）当成 ΔV（差分）喂进模型。二者不仅数值尺度不同，连符号和物理含义都不同；k_v 会被拟合成一个没有意义的数。

**② 解释变量重复 → 回归稀释**。S-1 12 天重访，一个月约 2–3 个干涉对。月内对的两端落在同一个月，取到同一个月值，`V(t_j) − V(t_i)` 被算成 **0**，而 γ_ij 却各不相同。于是样本里出现大量「同一个 x、不同 y」，最小二乘只能把 k_sm / k_v 往 0 压，残差里塞满未建模的月内变化。**这是最隐蔽也最致命的一条**：结果看上去能跑通，系数却是系统性低估的。

**③ 事件被平滑掉**。Yanco 是灌溉农业区，灌溉 / 降水脉冲的时间尺度是几天。月度合成（尤其最大值合成 MVC）把这些脉冲抹平，而 12 天相干性恰恰对它们最敏感 —— 等于把 L2 相对 L1 的增量信息全部丢掉。MVC 还有额外问题：它偏向月内最"绿"的一天，**对应的日期不确定**，与确定的 S-1 采集日之间可能有 2–3 周的物候偏移。

> **结论**：月值不是"低分辨率的逐景数据"，它是另一种观测量。不能升采样冒充端点，只能降采样响应变量去配合它。

---

## 3. 动手前先确认两件事（各 10 分钟）

### 3.1 月度合成用的是什么 reducer？

| reducer | 能否用作端点代理 | 说明 |
|---|---|---|
| `mean` / `median` | 勉强可作月度状态 | 需注意月内事件已被平均掉 |
| `max`（NDVI MVC） | **禁止** | 值系统性偏高、对应日期不确定、物候偏差最大 |
| 未知 | **按最坏情况处理** | 查 GEE 脚本里的 `reducer`，或看导出文件的 `system:time_start` 是否全是月首 |

配置里 `input_spec.monthly_composite` 目前是 `unknown`，**必须确认后填实**。

### 3.2 Δt 是否可变（决定 k_t 能不能标定）

ASF HyP3 **不替你选干涉对**：配对由 Vertex 的 Baseline / SBAS 工具或 SDK 的
`nearest-neighbor`（`depth` 最多 2）决定。因此：

- 若所有对都是 12 天 → **k_t·Δt 与截距 c 完全共线，k_t 不可识别**。
  此时 `least_squares` 会把任意值分给 k_t 而让 c 反向补偿，系数没有物理意义。
  **必须固定 k_t = 0**：

```python
from ecohydro_sm.insar.coherence import (
    PairCoherenceSMConfig, calibrate_decay_coefficients_pairwise,
)
pd = make_pair_dates(ref_dates, sec_dates)
cfg = PairCoherenceSMConfig(fix_k_t_zero=pd.dt_is_constant)
coef = calibrate_decay_coefficients_pairwise(gamma_pairs, dv, dsm_ref,
                                             pd.dt_days, cfg=cfg)
```

  `fix_k_t_zero=True` 会把参数向量降为 3 维、彻底去掉 Δt 项；
  **不是**把 k_t 的上下界都设成 0（那样仍留下退化自由度，结果依赖初值）。
  单元测试 `test_constant_dt_without_fix_is_degenerate` 固定了共线组合
  `c − k_t·12` 可识别、单项不可识别这一事实。
- 若混有 12 / 24 天对 → k_t 可识别，保留。

`00b_align_temporal.py` 的报告会自动给这个判定；`PairDates.dt_is_constant`
也可在代码里直接查。HyP3 产品名里的 `nnn` 字段就是时间 separation in days，
比解析两个日期更可靠。

---

## 4. 三条路线

| | **B 逐景插值到端点**（主结果） | **A 月度降采样**（稳健性） | **C 仅跨月对**（敏感性） |
|---|---|---|---|
| 响应变量 | 保持逐对 γ_ij | γ_ij → 月度 γ_M | 逐对 γ_ij |
| 解释变量 | 逐景光学插值出 V(t_i)、V(t_j) 后取端点差分 | **月际差分** ΔV_M = V_M − V_{M−1} | 月际差分近似 ΔV_ij |
| 时间支撑 | 一致（12 天 ↔ 12 天） | 一致（30 天 ↔ 30 天） | **仍错配**（30 天 ↔ 12 天） |
| 样本量 | n_pairs | n_months（2 年仅 24 个） | 约 n_pairs × 40% |
| 事件响应 | 保留 | 丢失 | 部分丢失 |
| 额外成本 | 需在 GEE 重导**逐景**产品 | 无 | 无 |
| 定位 | 主结果 | 稳健性检验 / 敏感度分析 | 补充说明，**不能作主结果** |

**推荐组合**：B 作主结果、A 作稳健性检验并列报告。两者之差本身就是"时间支撑敏感度"的定量结果，可以写进论文。C 只在完全无法重导数据时作为补充，且必须明写 limitation。

**不要试图"把月值插值成逐日"再喂给逐对模型** —— 那只是给 30 天的平均值套上一条平滑曲线，制造伪精度，不增加任何信息。

---

## 5. 用站点诊断决定 A/C 是否可用（先做这一步）

不要拍脑袋选路线。OzNet 是 20 min 观测，用它量化月度合成丢了多少信息：

```bash
python scripts/00b_align_temporal.py --root . \
    --station-csv data_raw/yanco/oznet/YA3_0-5cm.csv \
    --time-col timestamp --value-col sm_0_5cm
```

诊断逻辑（`temporal_support_diagnostic()`）：

1. 由站点算每个干涉对的真实 ΔSM_ij；
2. 把同一序列聚成月均值，取月际差分作「月度代理」；
3. 报告二者的**相关系数平方**（尺度无关的解释方差）、`scale_ratio`、`rmse`。

判据：

- **R² < 0.5** → 月度合成解释不了逐对变化，**必须走 B**；
- **R² ≥ 0.5** → A/C 可作稳健性检验，但 `scale_ratio ≠ 1` 说明月度差分与逐对差分幅度不同，**月度模型必须重新标定 k_sm / k_v，不能复用逐对系数**。

为什么用相关系数平方而不是"直接拿月度代理当 ΔSM"的 R²：30 天差分与 12 天差分幅度本就不同，固定尺度比较会把**尺度失配**误判成**信息丢失**。

这个诊断的数字（"月度合成仅解释了 X% 的逐对土壤水变化"）本身就是一个可写进论文的定量论证。

---

## 6. 路线 B 详细操作

### 6.1 在 GEE 重导逐景产品

**不要**再做 `calendarRange + mean/max`。原来的脚本里把这段去掉，改成逐景导出：

```javascript
var roi = ee.Geometry.Rectangle([145.80, -35.25, 146.55, -34.65]);

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate('2021-01-01', '2023-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .map(function(img) {
    var qa = img.select('SCL');
    var ok = qa.neq(3).and(qa.neq(8)).and(qa.neq(9))
               .and(qa.neq(10)).and(qa.neq(11));
    var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
    var ndmi = img.normalizedDifference(['B8', 'B11']).rename('NDMI');
    return img.addBands([ndvi, ndmi]).updateMask(ok)
      .set('scene_date', img.date().format('YYYY-MM-dd'));
  });

// 逐景导出：把日期列表拉到客户端后循环提交任务
s2.aggregate_array('scene_date').evaluate(function(dates) {
  dates.forEach(function(d) {
    var img = s2.filter(ee.Filter.eq('scene_date', d)).first()
                .select(['NDVI', 'NDMI']).clip(roi);
    Export.image.toDrive({
      image: img, description: 'S2_scene_' + d, folder: 'yanco_s2_scene',
      fileNamePrefix: 'S2_' + d.replace(/-/g, ''),
      region: roi, scale: 10, crs: 'EPSG:32755', maxPixels: 1e9
    });
  });
});
```

GRD 同理，但**两个额外要求**：

```javascript
var grd = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi).filterDate('2021-01-01', '2023-12-31')
  // 同一轨道、同一升/降轨、同一极化，否则月内混合入射角会引入伪变化
  .filter(ee.Filter.eq('relativeOrbitNumber', 55))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisations', 'VV'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'));
// 保留 VV + angle 两个波段，不做均值合成，逐景导出
```

导出文件名带 `YYYYMMDD`，仓库 `00_prepare_yanco.py` 会识别为 `granularity=scene`。

### 6.2 插值到干涉对端点

```python
from ecohydro_sm.preprocess.temporal_alignment import (
    make_pair_dates, pair_endpoint_difference, sar_acquisition_utc,
)

pd = make_pair_dates(ref_dates, sec_dates, dt_days=dt_days)
fvc = ...  # (n_scenes, n_pixels) 逐景 FVC / NDVI，第一维是时间
out = pair_endpoint_difference(fvc, scene_dates, ref_dates, sec_dates,
                               max_gap_days=10)
delta_v = out["abs_delta"]        # -> calibrate_decay_coefficients_pairwise
valid   = out["valid"]            # 剔除端点离最近有效光学超过 10 天的对
```

要点：

- **线性插值 + 窗口约束**。端点前后 ±10 天内没有有效光学观测就标记该对无效；
  云多的月份宁可丢样本，不要用 20 天外的景硬插。
- **记录 gap_days 并作为特征**（`out["gap_ref_days"]` / `gap_sec_days`），
  让模型看到时间错配，也便于做敏感性分析。
- 像元级云洞：插值时单侧无效退化为另一侧最近值，双侧无效置 NaN；
  pairwise 标定函数内部已有 `np.isfinite` 屏蔽，不必额外填充。
- **站点取过境时刻值**，不要用日平均：

```python
t_utc = sar_acquisition_utc("2021-01-15", local_hour=6.5,
                            tz_name="Australia/Sydney")
# 文件名日期是 UTC 日期，站点当地日期 = 它 + 1 天。
# 1 月为 AEDT(UTC+11)：当地 1/16 06:30 = 1/15 19:30 UTC
# 7 月为 AEST(UTC+10)：当地 7/16 06:30 = 7/15 20:30 UTC
# 用 ZoneInfo 按日期自动判定，别用固定 +10（夏令时会错 1 小时）
```

### 6.3 L1 / L2 比较前必须统一时间尺度

上一轮定的判据「L1 与 L2 差 < 0.06 m³/m³」在**两者时间尺度一致时**才成立。
若 L1 用的是月度 GRD 产品（月尺度 SM）、L2 是逐对反演积分出的逐日期 SM，
直接相减没有意义。二选一：

- 把 L2 的逐日期序列**聚合到月**（中位数）再与 L1 比；或
- L1 也改用**逐景** GRD 反演（推荐，与 B 路线一致，一次重导解决）。

---

## 7. 路线 A 详细操作

```python
from ecohydro_sm.preprocess.temporal_alignment import (
    aggregate_pairs_to_monthly, monthly_state_to_difference,
)

agg = aggregate_pairs_to_monthly(gamma_pairs, pd,
                                 assign_by="sec", reducer="median", min_pairs=2)
gamma_m = agg["gamma_monthly"]                       # (n_months, n_pixels)
dt_bar  = agg["dt_mean"]                             # 必须作为协变量保留
dv_m = monthly_state_to_difference(vwc_monthly, agg["months"], absolute=True)
dsm_m = monthly_state_to_difference(sm_monthly_ref, agg["months"], absolute=True)
```

六条规则：

1. **归属按副影像月份**（γ 反映"到 t_j 为止"的变化）；`assign_by` 可改但要全链一致。
2. **用中位数**聚合 γ，它对长基线 / 雨日异常对更稳健。
3. **解释变量必须做月际差分**：月度合成是状态量，直接喂给要 |ΔV| 的模型是量纲错误。
4. **保留 Δt̄ 协变量**：各月干涉对的平均时间基线不同（有的月 2 对 12 天、有的月混了 24 天），
   不建模会把时间去相干误判成含水量变化。若全部 Δt 相同则固定 k_t=0（见 §3.2）。
5. **缺月处的差分置 NaN**：跨两个月的差分时间支撑与相邻月不同（`require_contiguous_months`）。
6. **重新标定 k_sm / k_v**，不要复用逐对系数；月际差分的幅度与逐对差分不同。

**样本量警告**：2 年只有 24 个月，月度样本的有效独立数约等于月数本身。
必须报告置信区间，**不能报"逐日精度"或逐日 anomaly R**；A 路线的结论只能到月尺度。
`agg["gamma_std"]`（月内 γ 离散度）是现成的诊断量：值大说明月度聚合把显著的月内变化平滑掉了。

---

## 8. 路线 C 的限制（写进论文）

只保留 `ref_month != sec_month` 的对（`select_cross_month_pairs`），
用月际差分近似 ΔV_ij。这至少保证 ΔV 是真实的差分量，但：

- 12 天重访下跨月对约占 40%，样本量减半；
- 端点是月中任意日期 vs 月度状态，时间支撑仍是 30 天 对 12 天；
- Δt 与端点月份不对齐。

只能作敏感性分析，且必须声明。若诊断 R² < 0.5，C 路线没有意义。

---

## 9. 月值数据仍然有用的地方（别浪费）

月值 GRD / S-2 **不能支撑逐对相干性**，但可以用于：

- **L1 月尺度基线**：月值 σ⁰ 已相当于多景平均、做了斑点抑制，做月度 SM 反演反而更稳；
- **WCM 参数标定的月度约束**、干湿参考端元、季节物候曲线；
- 与 ERA5-Land / SMAP 月值产品做三重配置（时间尺度天然一致）；
- 路线 A 的月度 γ 建模。

一句话：**月值归月值，逐对归逐对，各走一条链，最后在同一时间尺度上比较。**

---

## 10. 常见错误清单

| # | 错误 | 后果 |
|---|---|---|
| 1 | 月值状态量直接当 ΔV_ij | 量纲错，k_v 无意义 |
| 2 | 月内对算 ΔV=0 | 回归稀释，k_sm/k_v 系统性偏低 |
| 3 | MVC 月度 NDVI 当端点值 | 物候偏移 2–3 周，值系统性偏高 |
| 4 | 月度产品时间戳取月首 | 凭空引入最多 30 天偏差（统一取月中 15 日） |
| 5 | Δt 全相同仍拟合 k_t | k_t 与 c 共线，系数不可识别 |
| 6 | 把月值插成逐日曲线冒充逐景 | 伪精度，不增加信息 |
| 7 | 用日平均站点值配 S-1 过境 | 平滑掉灌溉/降水响应 |
| 8 | 固定 UTC+10 换算澳洲时间 | 夏令时期间错 1 小时 |
| 9 | L1 月值与 L2 逐日期直接相减 | 尺度不一致，"差 < 0.06" 判据失效 |
| 10 | A 路线复用 B 路线的 k 系数 | 月际差分幅度不同，反演尺度错误 |

---

## 11. 代码接口速查

| 函数 / 脚本 | 作用 |
|---|---|
| `make_pair_dates()` | 由 ref/sec 日期建干涉对时间表，含 `dt_is_constant`、`cross_month_mask` |
| `interpolate_series_to_dates()` | 逐景线性插值到任意目标日期 + gap / 有效性 |
| `pair_endpoint_difference()` | 路线 B：求每个 pair 的端点差分（带符号 + 绝对值） |
| `select_cross_month_pairs()` | 路线 C：跨月对掩码 |
| `aggregate_pairs_to_monthly()` | 路线 A：γ_ij → γ_M，含 `n_pairs` / `gamma_std` / `dt_mean` |
| `monthly_state_to_difference()` | 路线 A/C：月度状态量 → 月际变化量，缺月置 NaN |
| `temporal_support_diagnostic()` | 用高频站点定量决定 A/C 是否可用 |
| `sar_acquisition_utc()` | S-1 文件名 UTC 日期 → UTC 过境时刻，ZoneInfo 处理 DST |
| `scripts/00b_align_temporal.py` | 读 manifest 产出 `pair_alignment.csv` + `alignment_report.md`，`--selftest` 可空跑 |
| `scripts/00_prepare_yanco.py` | 清点；已支持 `YYYYMM` / `YYYY_MM` 月度命名，自动报时间支撑与 k_t 可识别性 |

---

## 12. 论文里怎么写

建议把时间支撑问题正面写成方法贡献，而不是藏起来：

> 逐干涉对相干性（12 天差分）与光学 / SAR 幅度月度合成（30 天状态）具有不同时间支撑。
> 我们首先在 Yanco 用 20 min 站点观测量化了该错配：月度合成仅能解释逐对土壤水变化的
> R² = X（scale ratio = Y），据此在 GEE 重导逐景产品并插值到干涉对端点；
> 同时给出月度降采样结果作为时间支撑敏感度的对照，两者在 [指标] 上相差 Z。

这比"我们做了 S-1/S-2 融合"更能回应审稿人对时间错配的质疑，也把
INSAR4SM 单站验证未讨论的时间尺度问题变成了自己的增量。
