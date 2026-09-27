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
