# 工作流说明

## 0. 凭据管理（重要）

**任何凭据都不进仓库。** `.gitignore` 已屏蔽 `.netrc`、`.dodsrc`、`.env` 等。

| 服务 | 方式 |
|---|---|
| ASF / Earthdata（HyP3、Sentinel-1、NISAR） | `~/.netrc` 写入 `machine urs.earthdata.nasa.gov login <user> password <pass>` |
| Copernicus Data Space | 环境变量 `CDSE_USER` / `CDSE_PASSWORD` |
| CDS（ERA5-Land） | `~/.cdsapirc` |
| GitHub | `gh auth login`，不要写 token 到文件 |

## 1. 数据获取与预处理

### Sentinel-1 GRD → RTC

推荐直接用 **ASF HyP3 的 `RTC_GAMMA`**，与既有 MintPy 工作流一致：

```
job 参数：radiometry=gamma0, resolution=30, dem_matching=True,
         speckle_filter=True, include_dem=True, include_inc_map=True
```

本地 SNAP 处理的顺序不可颠倒：
`Apply Orbit File → Thermal Noise Removal → Calibration → GRD Border Noise Removal → Terrain Flattening → Terrain Correction → Speckle Filter → 分贝化`

### Sentinel-1 SLC → 干涉

```
HyP3 (insar_gburst / insar_isce) → prep_hyp3 → smallbaselineApp.py
```

- 干涉对只用**短时空基线**（≤ 12 天 / ≤ 150 m）
- 大气校正 ERA5 + pyaps3：**形变支路必须做，相干性支路不做**
- 相干性支路用大窗口多视（~200 m）降噪

### NISAR（可选，L 波段分支）

```
ASF 下载 GUNW → MintPy prep_nisar（≥1.6.4）→ smallbaselineApp.py
备选：ARIA-tools
```

NISAR 当前为 **Provisional**，2026 Q4 将重处理。**全流程必须写成一键重跑脚本**，
结论不能依赖单版本数据。

## 2. 执行顺序

```
01_prepare_s2.py              S2 指数 + 云/冻融/降水掩膜
02_build_synthetic_db.py      WCM 合成训练库
03_train_inversion.py         NN 反演器 + 实测残差二级校正
04_retrieve_sm.py             L1 表层土壤水
05_insar_coherence_sm.py      L2 相干性独立约束
06_deformation_to_storage.py  L4 形变 → 储变量
07_water_balance_closure.py   L3 水量平衡闭合 → 根区土壤水
08_validate.py                验证（不可跳过）
```

每个脚本都支持 `--selftest`，可在**没有任何真实数据**的情况下
用合成数据验证算法链路：

```bash
for f in scripts/0*.py; do python $f --config configs/minqin.yaml --selftest; done
```

## 3. 用 GEE 处理 Sentinel-2（推荐）

脚本：`scripts/gee_s2_indices.js`（GEE Code Editor）。

### 哪些环节适合 GEE，哪些不适合

| 环节 | 是否上 GEE | 原因 |
|---|---|---|
| S-2 云掩膜、NDVI/NDMI/FVC | ✅ **强烈推荐** | 免去几十 GB 下载；逐像元分位端元在 GEE 上几秒完成 |
| 站点抽表（构建训练/验证表） | ✅ **强烈推荐** | `reduceRegions` 直接出表，完全不用导出栅格，这是效率提升最大的一步 |
| Sentinel-1 RTC | ❌ 不推荐 | GEE 的 S-1 只有 GRD 且**无地形辐射校正**，生产请用 HyP3 RTC |
| **InSAR 相干性与形变（L2/L4）** | ❌ **不可能** | GEE **没有 Sentinel-1 SLC**，相位信息不可用，必须本地 HyP3 + ISCE2/MintPy |

**所以整体是混合架构**：GEE 负责 S-2，本地负责 S-1 与全部 InSAR。

### 三个必须注意的点

1. **用 `COPERNICUS/S2_SR_HARMONIZED`，不要用 `COPERNICUS/S2_SR`**。
   后者在 2022-01 前后有 +1000 DN 的处理基线偏移，混用会让 NDVI 时序出现虚假跳变。
2. **云掩膜要 SCL + s2cloudless 双保险**。SCL 会低估薄卷云与云影，单靠它会在
   高反射地表上留下大量污染像元。
3. **导出前必须与 S-1 RTC 网格对齐**（`crs` + `crsTransform` 从参考 RTC 影像读取）。
   不对齐会在融合阶段引入系统性配准误差，且这种误差很难在精度指标上被发现。

### 效率与代价

| | 本地 | GEE |
|---|---|---|
| 指数与掩膜计算 | 需先下载原始数据（几十 GB） | 分钟级，无需下载 |
| 逐像元 NDVI 分位端元 | 需自己写分块统计 | `ee.Reducer.percentile` 一行 |
| 站点抽表 | 需先导出栅格再采样 | 直接出 CSV，秒级 |
| 栅格导出 | — | **瓶颈**：任务排队，几百景会很慢 |

**代价**：GEE 免费额度限非商业用途；数据集会随处理基线更新而变化，
投稿前必须记录 Collection ID 与访问日期，并把导出产品归档。

## 4. 本地标定是强制项

以下参数**没有通用值**，必须用在地观测标定，否则结果不可用：

| 参数 | 标定数据 |
|---|---|
| WCM 的 A / B / E | 生物量样方 + 同步 σ0 |
| VWC–NDVI 关系 | 地上生物量实测 |
| NDVI 端元（ndvi_soil / ndvi_veg） | 研究区 NDVI 直方图 5%/95% 分位 |
| S_ke | 季节形变幅度 ÷ 水位幅度 |
| S_y | 抽水量 / GRACE / 岩性，并给区间 |
| k_sm / k_v | 参考土壤水时序联合标定 |

## 5. 常见错误

| 症状 | 原因 |
|---|---|
| R² 高达 0.8 但换区就崩 | 随机 CV 的空间泄漏；改空间块 CV |
| 生长季误差突然变大 | 结构去相干未被分离；启用 L2 的结构去混 |
| 把总形变当储变量 | 未分离永久压实；先做谐波+趋势分解 |
| ΔS 量级离谱 | 直接把 Δd 当 ΔS；用 S_y/S_ke 并给区间 |
| 产品"1 km"但细节模糊 | 标称分辨率与实际聚合不一致 |
