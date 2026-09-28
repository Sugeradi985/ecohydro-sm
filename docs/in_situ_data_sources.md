# 已有研究案例用了什么站点数据？能不能公开拿到？

> 配套文档：`station_requirements.md`（本站点的要求与布设）、`experiment_runbook.md`（处理步骤与判据）

---

## 1. 一句话结论

**能公开拿到的，基本都是"长期业务观测网"（ISMN、SMOSMANIA、黑河观测网）；
而绝大多数高分论文的验证数据来自"短期密集野外试验（ground campaign）"，
这类数据通常不在任何公开库里，只能发邮件向作者/承担机构索取。**

所以文献里报的 RMSE，你未必能复现——因为**你拿不到它的验证数据**。
这不是文献不严谨，而是土壤水遥感领域的结构性现状。

---

## 2. 案例 × 站点数据 × 可获取性（主表）

| 研究 | 站点数据 | 规模 | 可获取性 |
|---|---|---|---|
| **INSAR4SM**<br>Karamvasis &amp; Karathanassi 2023<br>*Comput. Geosci.* 178:105410 | ISMN 站点（California/Arizona 干旱区） | **1 个站点** | ✅ ISMN 公开（ismn.earth 注册下载） |
| 同上（辅助对比） | ERA5-Land 再分析 | 栅格 | ✅ 公开 |
| **S²MP 1 km 验证**<br>Madelon et al. 2023<br>*HESS* 27:1221–1242 | ISMN 站点 | 多站 | ✅ ISMN 公开 |
| 同上 | **MERGUELLIL（突尼斯 Kairouan）** | thetaprobe 网，5 / 40 cm，裸土 | ❌ **需发邮件向作者索取**（mehrez.zribi@ird.fr） |
| 同上（被检验的产品） | **S²MP 产品本身 + 生产代码** | 1 km | ❌ **需发邮件索取**（nicolas.baghdadi@inrae.fr） |
| **500 m 融合反演**<br>*Water* 12(3):866, 2020 | **Occitanie（法国南部）**：SMOSMANIA 21 个自动站剖面（5/10/20/30 cm），文中只用 5 cm 的 3 个站 | 21 站中取 3 | ✅ ISMN 公开 |
| 同上 | **Merguellil（Kairouan 平原）**：连续 thetaprobe 站，裸土，5 / 40 cm，重量法标定 | 站点网 | ⚠️ 走 CESBIO OSR 门户（osr-cesbio.ups-tlse.fr），需注册/申请 |
| 同上 | **Banizombou（Niger）**：5 cm 连续观测 | 单站 | ✅ ISMN 公开 |
| 同上 | **Ouémé（Benin）**：AMMA-CATCH 观测网，5 cm | 站点网 | ✅ ISMN 公开（AMMA-CATCH 也可单独申请） |
| **灌溉与土壤水制图**<br>Bousbih et al. 2018<br>*Remote Sens.* 10(12):1953 | **20 个谷物参考田块**的地面同步观测（土壤水 + 粗糙度 + LAI） | 20 田块 × 2015–2017 | ❌ **短期野外试验，不在公开库**，需向作者/IRSTEA-INRAE 索取 |
| **Hydroterra+ / 4DMED-DEMETRAS**<br>De Zan, Filippucci &amp; Brocca 2026<br>*RSE*（doi 10.1016/j.rse.2026.115266） | 7 个试验区，与陆面模式、其他遥感产品及 in situ 比较 | 7 区 | ⚠️ 部分公开，in situ 需按区索取 |

---

## 3. 三个必须看清的事实

### 3.1 INSAR4SM 的 RMSE 0.029 / R 0.78，来自 1 个站点

原文写的是 *"...assessed with independent SSM observations from **a station** of the ISMN"*——**单数**。

| 版本 | in situ RMSE | in situ R | ERA5-Land RMSE | ERA5-Land R |
|---|---|---|---|---|
| arXiv 预印本（2210.10665v1） | 0.027 | 0.88 | 0.035 | 0.71 |
| 正式刊出版（2023） | **0.029** | **0.78** | **0.049** | **0.62** |

**怎么读这个数字**：
- 它是"该方法可行"的**存在性证明**，不是"精度达 0.029"的普适结论；
- 单站验证无法做空间交叉验证，也无法做 VWC 分层，泛化能力未被检验；
- 预印本 → 正式版 R 从 0.88 掉到 0.78，本身就是"换/扩充验证数据后精度回落"的典型表现。

**对本课题的意义**：这恰恰是本方案的差异化空间——
用 **≥ 10 站 × ≥ 2 年 + 空间块 CV + 按 VWC 分层** 给出相干性路线的精度，
是可以站得住的增量贡献。但前提是你得先有站点。

### 3.2 "高分论文"的验证数据大多是短期的、不公开的

Zribi / Baghdadi / El Hajj / Bousbih 这一系列的共同模式：

- 研究期 **2–3 年**，地面观测是**与卫星过境同步的密集野外观测**（查墒 + 粗糙度 + LAI）；
- 这类数据有明确的**项目资助期**（CHAAMS、TOSCA/CNES、ESA-ESTEC ITT），
  项目结束后数据不进 ISMN，也不做长期维护；
- 2023 年的 HESS 论文里，连 **S²MP 产品本身和代码**都只能通过"发邮件给 Baghdadi"获取——
  这说明该领域的数据共享基础设施确实薄弱。

**结论**：想复现这些精度，不现实；想与之对比，只能用**同一研究区的公开子集**做公平对比，
或者干脆把对比对象换成公开产品（SMAP/SMOS/ESA CCI/CGLS）。

### 3.3 中国公开源的口径与你想的不一样

| 源 | 口径 | 现实情况 |
|---|---|---|
| **中国气象数据网**（data.cma.cn） | 土壤水分自动站 | **日值/逐时值不对外开放**。官方答复：因自动站稳定性与数据质量问题，"主要向科学数据用户开放**土壤墒情旬值及以上数据产品**"。旬值（约 10 天一次，10 cm）勉强匹配 S-1 重访，但**事件级响应检验做不了** |
| **国家生态科学数据中心** | 民勤站 2021–2023 水环境要素日尺度数据集（含土壤含水量） | ⚠️ **保护期至 2028-03-28**，之后才可获取；当前为"协议共享"（签协议后下载） |
| **国家冰川冻土沙漠科学数据中心** | 石羊河尾闾柽柳降水改变试验样地土壤湿度（2020，民勤县，5TM 传感器，日尺度） | ⚠️ 申请获取；数据量小（302 KiB），单点对照试验性质 |
| 同上 | 民勤县西渠镇试验池数据集（2019–2023，土壤水热盐） | ⚠️ 申请获取 |
| **国家生态科学数据中心** | 2002 年石羊河流域土壤特征参数（30 采样点，民勤 12 点） | ✅ **公开共享**，填用途直接下载；土壤质地/饱和含水率/容重/水分特征曲线 → 可用于 WCM 与 Dubois 的本地参数化 |
| **国家青藏高原科学数据中心** | 黑河上游八宝河生态水文传感器网（**40 节点，4/20 cm，逐时，2013–2017**） | ✅ **开放获取** |
| 同上 | HiWATER WATERNET 中游（**50 节点，4/10 cm，10 min，2012**） | ✅ **开放获取** |
| 同上 | 黑河流域地表过程综合观测网（11–15 站，含**宇宙射线仪 CRNS** 区域土壤水、土壤温湿廓线至 160 cm） | ⚠️ 申请获取；2013–2026 已发布数据集 564 个 |

---

## 4. 可获取性分级（照这个顺序去找）

**A 级 — 注册/开放即下载，优先用尽**
- ISMN（ismn.earth）：全球聚合，含 SMOSMANIA、Banizombou、Ouémé、黑河相关网络（CTP_SMTMN 56 站、MAQU 20 站）
- 国家青藏高原科学数据中心：黑河 WATERNET / 八宝河传感器网（**开放获取**）
- 国家生态科学数据中心：石羊河 2002 土壤参数
- ERA5-Land、SMAP/L3、SMOS、ESA CCI、CGLS（作为 TC 的第三方与对比基准）

**B 级 — 需申请/签协议，提前 1–2 个月启动**
- 黑河综合观测网年度数据集（含 CRNS）
- 民勤站水环境要素（**保护期到 2028-03-28**；或走协议共享谈判）
- 石羊河尾闾柽柳样地、西渠镇试验池（ncdc.ac.cn 申请）

**C 级 — 发邮件向作者索取，成功率不确定**
- Merguellil：mehrez.zribi@ird.fr
- S²MP 产品与代码：nicolas.baghdadi@inrae.fr
- Bousbih 2018 的 20 个参考田块观测：经 IRD / INRAE 转询

**D 级 — 靠自己**
- 民勤盆地自建剖面站（见 `station_requirements.md` §5）
- 石羊河流域管理局 / 甘肃省水文部门的地下水位井与灌溉记录（**L4 标定绕不开**）

---

## 5. 对民勤课题的现实建议：双研究区

民勤盆地公开土壤水站点**不足以支撑投稿级验证**（需要 15–20 站 × 2 年）。
最务实的布局是把"方法学验证"和"科学应用"拆到两个区：

| | **方法学区：黑河中游（张掖 / 大满灌区）** | **目标区：民勤盆地** |
|---|---|---|
| 目的 | 验证 L1/L2 反演与验证纪律本身 | 回答"绿化—地下水—沉降"科学问题 |
| 地面数据 | WATERNET 50 节点、八宝河 40 节点（开放获取）+ CRNS 区域土壤水（申请） | 自建 3–5 个剖面站 + 气象局旬值 + 申请到的民勤站数据 |
| 优势 | 同为河西走廊沙漠绿洲灌溉农业，地表过程与民勤高度可比；**数据现成、密度够** | 有超采—治理的反向自然实验，科学故事强 |
| 论文表述 | "方法在河西走廊绿洲灌溉区经密集观测网验证" | "应用于民勤盆地，闭合 P−ET−R−ΔS" |

**这样做的好处**：
1. 不用等民勤站点建好就能发方法学论文（黑河数据现在就能下载）；
2. 民勤站的 3–5 个自建点，其定位变成"**目标区的一致性检验**"而非"训练/验证样本"，
   数量门槛从 15–20 降到 3–5，成本可控；
3. 黑河的 CRNS 数据正好用来解决**尺度代表性**问题（CRNS 半径 130–240 m，与 100–500 m 产品匹配最好），
   这是本方案 §4 里最容易被审稿人挑的点。

**兜底**（若两端数据都拿不到）：三重配置（Triple Collocation）不依赖真值。
但要记住它的根本限制——**三套数据无法识别误差相关性**，
已在 `validation/triple_collocation.py` 中显式记录并写进测试。

---

## 6. 索取/申请的操作要点

1. **先确认对方是否还在维护**：项目结束 5 年以上的数据，作者多半已换方向，回复率低；
   优先联系仍在运行的观测网（黑河、CESBIO OSR）。
2. **邮件里说清三件事**：用途（论文/基金）、需要的时间范围与深度层、是否接受共同署名；
   主动提出在致谢中标注 DOI/CSTR，会显著提高回复率。
3. **国内平台申请**要预留时间：ncdc.ac.cn / nesdc.org.cn / tpdc.ac.cn 的"申请获取"通常 1–5 个工作日，
   但**保护期内的数据集无法提前获取**（民勤站数据集即属此类，2028-03-28 后开放）。
4. **拿到后第一件事**：核对深度层定义、单位（m³/m³ vs %）、是否做过土壤标定、
   以及时间戳是地方时还是 UTC。这四项错一个，后续全部作废。
5. **引用与致谢**：国内数据集通常要求在成果中标注"数据来源于 ×× 数据中心"并将成果反馈至指定邮箱，
   写论文时不要漏。

---

## 7. 关键文献（补充至 `references.md`）

- Karamvasis, K., Karathanassi, V. (2023). Soil moisture estimation from Sentinel-1 interferometric observations over arid regions. *Computers &amp; Geosciences*, 178, 105410. doi:10.1016/j.cageo.2023.105410（预印本 arXiv:2210.10665）
- Madelon, R., Rodríguez-Fernández, N.J., Bazzi, H., Baghdadi, N., Albergel, C., Dorigo, W., Zribi, M. (2023). Soil moisture estimates at 1 km resolution making a synergistic use of Sentinel data. *HESS*, 27, 1221–1242. doi:10.5194/hess-27-1221-2023
- Bousbih, S., Zribi, M., El Hajj, M., Baghdadi, N., Lili-Chabaane, Z., Gao, Q., Fanise, P. (2018). Soil moisture and irrigation mapping in a semi-arid region, based on the synergetic use of Sentinel-1 and Sentinel-2 data. *Remote Sensing*, 10(12), 1953. doi:10.3390/rs10121953
- El Hajj, M., Baghdadi, N., Zribi, M., Bazzi, H. (2017). Synergic use of Sentinel-1 and Sentinel-2 images for operational soil moisture mapping at high spatial resolution over agricultural areas. *Remote Sensing*, 9(12), 1292.
- Estimating 500-m resolution soil moisture using Sentinel-1 and optical data synergy. *Water*, 12(3), 866, 2020.（Occitanie / Merguellil / Banizombou / Ouémé 四站）
- De Zan, F., Filippucci, P., Brocca, L. (2026). *Remote Sensing of Environment*. doi:10.1016/j.rse.2026.115266
- 黑河上游生态水文传感器网络逐时土壤水分观测数据集（2013–2017）. 国家青藏高原科学数据中心. doi:10.11888/Hydro.tpdc.271137
- HiWATER: 黑河流域中游生态水文无线传感器网络 WATERNET 观测数据集（2012）. 国家青藏高原科学数据中心. doi:10.3972/hiwater.118.2013.db
- 民勤站 2021–2023 年水环境要素日尺度数据集. 国家生态科学数据中心. doi:10.12199/nesdc.ecodb.2021YFF0703900.mcsos.2025.11（保护期至 2028-03-28）
- 石羊河尾闾荒漠绿洲过渡带柽柳降水改变试验样地土壤湿度数据集（2020 年）. 国家冰川冻土沙漠科学数据中心. doi:10.12072/ncdc.nieer.db3924.2023
- 2002 年甘肃省石羊河流域土壤特征参数数据集. 国家生态科学数据中心. doi:10.12199/nesdc.ecodb.wwa.mon.018（**公开共享**）
