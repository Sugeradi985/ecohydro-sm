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

## 3. 本地标定是强制项

以下参数**没有通用值**，必须用在地观测标定，否则结果不可用：

| 参数 | 标定数据 |
|---|---|
| WCM 的 A / B / E | 生物量样方 + 同步 σ0 |
| VWC–NDVI 关系 | 地上生物量实测 |
| NDVI 端元（ndvi_soil / ndvi_veg） | 研究区 NDVI 直方图 5%/95% 分位 |
| S_ke | 季节形变幅度 ÷ 水位幅度 |
| S_y | 抽水量 / GRACE / 岩性，并给区间 |
| k_sm / k_v | 参考土壤水时序联合标定 |

## 4. 常见错误

| 症状 | 原因 |
|---|---|
| R² 高达 0.8 但换区就崩 | 随机 CV 的空间泄漏；改空间块 CV |
| 生长季误差突然变大 | 结构去相干未被分离；启用 L2 的结构去混 |
| 把总形变当储变量 | 未分离永久压实；先做谐波+趋势分解 |
| ΔS 量级离谱 | 直接把 Δd 当 ΔS；用 S_y/S_ke 并给区间 |
| 产品"1 km"但细节模糊 | 标称分辨率与实际聚合不一致 |
