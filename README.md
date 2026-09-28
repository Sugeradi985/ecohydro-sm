# ecohydro-sm

SAR + 光学融合的生态水文土壤水生产管线，并用 InSAR 相干性与形变把产品从表层（0–5 cm）推进到根区与含水层储变量。

## 设计立场

不做"第 N 个 Sentinel-1/2 融合土壤水产品"。表层融合已经很成熟，增量有限。
本仓库的目标是**表层—根区—含水层三层水分与储量闭合**：

| 层 | 产品 | 手段 | 代码位置 |
|---|---|---|---|
| L1 | 表层土壤水 0–5 cm | HyP3 RTC γ⁰ + S-2 植被 → WCM 合成库 + NN 反演 + 残差校正 | `src/ecohydro_sm/retrieval/` |
| L2 | 相干性独立表层约束 | γ 时序 + S-2 结构变化去混 + 闭合相位质控 | `src/ecohydro_sm/insar/` |
| L3 | 根区土壤水 0–1 m | 水量平衡残差沿深度分配 | `src/ecohydro_sm/fusion/water_balance.py` |
| L4 | 含水层储变量 ΔS | 形变 → S_ke/S_y → ΔS | `src/ecohydro_sm/insar/storage.py` |

## 安装

```bash
pip install -r requirements.txt
# 或
pip install -e .
```

栅格 I/O 依赖 `rasterio`/`xarray`，为可选依赖：只在 `io/` 模块中导入，**核心算法模块只依赖 numpy / scipy / scikit-learn**，便于纯数值单元测试与合成数据试验。

## 快速开始

```bash
# 1. 由 Sentinel-2 生成植被指数与掩膜
python scripts/01_prepare_s2.py --config configs/minqin.yaml

# 2. 用 WCM 正向模型生成合成训练库
python scripts/02_build_synthetic_db.py --config configs/minqin.yaml

# 3. 训练 NN 反演器（并用实测站点做残差校正）
python scripts/03_train_inversion.py --config configs/minqin.yaml

# 4. 生产 L1 表层土壤水
python scripts/04_retrieve_sm.py --config configs/minqin.yaml

# 5. L2 相干性支路
python scripts/05_insar_coherence_sm.py --config configs/minqin.yaml

# 6. L4 形变 → 储变量
python scripts/06_deformation_to_storage.py --config configs/minqin.yaml

# 7. L3 水量平衡闭合 → 根区土壤水
python scripts/07_water_balance_closure.py --config configs/minqin.yaml

# 8. 验证（空间块交叉验证 + 三重配置 + 分层报告）
python scripts/08_validate.py --config configs/minqin.yaml
```

## 验证纪律（写死在管线里，不是可选项）

现有大量文献用随机划分做交叉验证，存在空间自相关泄漏，R² 普遍虚高（严格分块后常从 0.7–0.8 掉到 ~0.5）。
本仓库强制：

- **空间块交叉验证**（`validation/spatial_cv.py`），块大小 ≥ 10 km 或按轨道/站点分块；
- **leave-one-station-out**；
- **三重配置**（`validation/triple_collocation.py`），用于站点稀疏区的无真值精度评估；
- **按 VWC 区间分层报告**，禁止只报一个全样本数值；
- 指标分离：`bias` / `ubRMSE` / `cRMSE` / `R` / `anomaly R`。

## 适用域声明（必须写在每个产品上）

- L1：全覆被可用，但 VWC > 3 kg/m² 时精度显著退化，须打标。
- L2：**仅适用于 VWC < 2 kg/m² 的裸土、稀疏植被与干旱半干旱区**。农田绿洲内部去相干，不承诺精度。
- L4：需要水位井或分层标做标定；`S_ke`/`S_y` 跨数量级，**结果必须给区间 + 敏感性分析**，禁止单点估计。

## 目录结构

```
configs/            配置（研究区、参数、路径）
src/ecohydro_sm/
    io/             栅格读写（可选依赖 rasterio）
    preprocess/     斑点滤波、光谱指数、掩膜（云/冻融/入射角）
    retrieval/      WCM 正向模型、NN 反演、变化检测、残差校正
    insar/          相干性、闭合相位、形变分解、储变量
    fusion/         水量平衡闭合与根区分配
    validation/     指标、空间块 CV、三重配置
    utils/          地理与绘图工具
scripts/            按执行顺序编号的生产脚本
scripts/gee_s2_indices.js  Google Earth Engine 上的 S-2 指数与掩膜处理
tests/              纯数值单元测试（不需要影像数据）
docs/               工作流说明与实验手册
docs/station_requirements.md  地面站点的要素、数量、代表性与获取途径
docs/in_situ_data_sources.md  已有研究案例用了什么站点数据、各自能否公开获取
```

## 混合架构：GEE + 本地

| 环节 | 平台 | 说明 |
|---|---|---|
| Sentinel-2 云掩膜、NDVI/NDMI/FVC、站点抽表 | **GEE** | 见 `scripts/gee_s2_indices.js`，免下载、分钟级 |
| Sentinel-1 RTC | 本地（HyP3） | GEE 的 S-1 无地形辐射校正 |
| InSAR 相干性与形变（L2/L4） | **必须本地** | GEE 没有 S-1 SLC，无相位信息 |

## 引用

若使用本仓库，请同时引用方法学来源：Rahmati et al. (2026, RSE) Sentinel-1 土壤水综述；
Karamvasis & Karathanassi (2023, Computers & Geosciences) INSAR4SM；
De Zan & Gomba (2018, RSE) 闭合相位反演。详见 `docs/references.md`。

## 许可

MIT
