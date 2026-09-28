# 参考文献

写代码前回原文核对数值；数值用于论文时务必引用原文而非二次转述。

## 综述（优先读）

- **Rahmati et al., 2026**, *Remote Sensing of Environment* —《Soil moisture retrieval from Sentinel-1: Lessons learned after more than a decade in orbit》, doi:10.1016/j.rse.2025.115146。当前最权威的 Sentinel-1 土壤水综述；明确指出干涉可观测（相干性/闭合相位）在干旱半干旱区有潜力，并把多源融合（Sentinel-2、SMAP、NISAR）列为优先方向。
- Brocca et al., 2024 — "the rush to high resolution"：指出多数号称 1 km 的产品并非真 1 km，只有 Sentinel-1 系的 1 km 是真实分辨率。

## Sentinel-1/2 融合反演

- Veloso et al., 2017 — S-1 + S-2 植被参数耦合 WCM。
- Zribi et al., 2019 — WCM 分离植被贡献，R² 0.65–0.75。
- El Hajj et al., 2017 — S-1 + S-2 神经网络反演。
- Lozac'h et al., 2020 (M2GARSS) — **S²MP** 业务化产品：输入 = VV σ0 + 入射角 + NDVI，Orfeo Toolbox 流水线。
- Bazzi et al., 2019/2020 — 回归树与随机森林融合。
- Bauer-Marschallinger et al., 2018 — TU Wien / ESA CCI Sentinel-1 1 km 表层土壤水。
- Faridani et al., 2025, *Agricultural Water Management* — 变化检测在葡萄园的 20 m 应用。
- Massart et al., 2025, *Agricultural Water Management* — 莫桑比克 500 m 变化检测产品。

## InSAR 相干性与闭合相位

- **De Zan & Gomba, 2018**, *RSE* — 由 SAR 闭合相位反演植被与土壤水。
- **Karamvasis & Karathanassi, 2023**, *Computers & Geosciences* 178:105410 — **INSAR4SM**：干旱区 RMSE 0.029 m³/m³、R 0.78、~250 m，配套开源工具箱。
- De Zan, Filippucci & Brocca, 2026, *RSE* — 干涉法多年多站验证：高相干区相关性强，低相干/植被/雪盖区退化。
- Walker et al., 2025, *Remote Sensing* 17(21) — 苏格兰泥炭地 InSAR 相干性与实测土壤水，春夏季 R² 达 0.83。

## 形变 → 地下水储变量

- Hoffmann et al., 2001 — 弹性骨架储水系数与 InSAR 形变的经典工作。
- Chaussard et al., 2014 — Santa Clara Valley 的形变—水位联合分析。
- Zhang X. et al., 2025, *RSE* — 华北平原 InSAR-VSM：2 km 地下水亏损，S_ke 2×10⁻⁴–2.1×10⁻²。
- Shi X. et al., 2025 — 唐山 2012–2021 十年形变约束地下水体积损失。
- Hengshui 案例（2025）— 弹性骨架储水系数与滞后校正。

## 正向散射模型

- Topp et al., 1980 — 介电常数混合模型。
- Dubois et al., 1995 — 裸土同极化后向散射半经验模型（适用域：k·h_rms ≤ 2.5、θ ≥ 30°、mv ≤ 35%）。
- Attema & Ulaby, 1978 — 水云模型（WCM）。

## 数据产品

- NISAR L3 SME2 — 200 m 土壤水，DSG/PMI/TSR 三算法平均，精度目标 0.06 m³/m³，VWC > 5 kg/m² 需打标。文档：nisar-docs.asf.alaska.edu/sme2
- ESA CCI Soil Moisture、SMAP L3、ERA5-Land、GLEAM / PML_V2 ET、SoilGrids 2.0。

## 地面观测数据与可获取性（详见 docs/in_situ_data_sources.md）

- **ISMN**（ismn.earth，注册下载）— 全球聚合；SMOSMANIA（法国南部 21 站剖面）、Banizombou（Niger）、Ouémé（Benin）均在其中。
- **Madelon et al., 2023**, *HESS* 27:1221–1242, doi:10.5194/hess-27-1221-2023 — S²MP 1 km 验证。**Merguellil 数据需发邮件向 mehrez.zribi@ird.fr 索取；S²MP 产品与代码需向 nicolas.baghdadi@inrae.fr 索取**（论文 Data/Code availability 原文）。
- **Bousbih et al., 2018**, *Remote Sensing* 10(12):1953, doi:10.3390/rs10121953 — Kairouan 平原，**20 个谷物参考田块**的地面同步观测（土壤水 + 粗糙度 + LAI），短期野外试验，不进公开库。
- *Water* 12(3):866, 2020 — 500 m 融合，四站：Occitanie（SMOSMANIA）、Merguellil（CESBIO OSR: osr-cesbio.ups-tlse.fr）、Banizombou、Ouémé。
- **INSAR4SM 的精度来自 1 个 ISMN 站点**（原文 "a station"，单数）；预印本 arXiv:2210.10665v1 报 RMSE 0.027/R 0.88，正式版为 0.029/0.78、ERA5-Land 0.049/0.62。
- **黑河上游八宝河**生态水文传感器网络逐时土壤水分（2013–2017，40 节点，4/20 cm）— 国家青藏高原科学数据中心，**开放获取**，doi:10.11888/Hydro.tpdc.271137。
- **HiWATER WATERNET** 黑河中游（2012，50 节点，4/10 cm，10 min）— 国家青藏高原科学数据中心，**开放获取**，doi:10.3972/hiwater.118.2013.db。
- **黑河流域地表过程综合观测网**（11–15 站，含宇宙射线仪 CRNS 区域土壤水、土壤温湿廓线至 160 cm）— 逐年发布，申请获取。
- **民勤站 2021–2023 年水环境要素日尺度数据集**（含土壤含水量）— 国家生态科学数据中心，doi:10.12199/nesdc.ecodb.2021YFF0703900.mcsos.2025.11，**保护期至 2028-03-28**，协议共享。
- **石羊河尾闾柽柳降水改变试验样地土壤湿度数据集（2020，民勤县，5TM，日尺度）** — 国家冰川冻土沙漠科学数据中心，doi:10.12072/ncdc.nieer.db3924.2023，申请获取。
- **2002 年甘肃省石羊河流域土壤特征参数数据集**（30 采样点，民勤 12 点；质地/容重/饱和含水率/水分特征曲线）— 国家生态科学数据中心，doi:10.12199/nesdc.ecodb.wwa.mon.018，**公开共享**，可直接用于本地参数化。
- **中国气象数据网**（data.cma.cn）— 土壤水分自动站**日值/逐时值不开放**，只开放**土壤墒情旬值及以上产品**。
