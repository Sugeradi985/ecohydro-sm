// ============================================================
// Sentinel-2 指数与掩膜批量处理（Google Earth Engine Code Editor）
// ============================================================
// 用途：为本仓库 L1/L2 支路生成分析就绪的 S-2 产品（NDVI / NDMI / FVC / 有效掩膜）
//
// 使用前必须修改的 4 处：
//   1. AOI          —— 研究区范围
//   2. START / END  —— 时间范围（建议 >= 2 年）
//   3. TARGET_CRS / TARGET_TRANSFORM —— **必须与你的 Sentinel-1 RTC 参考网格一致**
//   4. STATIONS     —— 地面站点（用 reduceRegions 直接抽表，无需导出栅格）
//
// 关键设计：
//   - 用 S2_SR_HARMONIZED 而非 S2_SR：后者在 2022-01 前后存在 +1000 DN 偏移，
//     直接混用会让 NDVI 时序出现虚假跳变。
//   - 云掩膜用 SCL + s2cloudless 云概率**双保险**：SCL 本身会低估薄卷云与云影。
//   - NDVI 端元（裸土/全植被）按像元逐点取 5%/95% 分位，不做全研究区统一取值，
//     这正是实验手册里"必须按研究区直方图重取端元"的要求。
// ============================================================

var AOI = ee.Geometry.Rectangle([102.9, 38.5, 103.9, 39.3], 'EPSG:4326', false);
var START = '2020-01-01';
var END = '2024-12-31';

// ---- 与 Sentinel-1 RTC 参考网格对齐（务必替换为你自己的值）----
// 做法：把一景 HyP3 RTC 产品上传到 GEE Asset，用
//   ee.Image(asset).projection() 读取 crs 与 transform，填到下面。
var TARGET_CRS = 'EPSG:32647';
var TARGET_TRANSFORM = [20, 0, 400000, 0, -20, 4300000]; // [xRes, 0, xMin, 0, -yRes, yMax]
var SCALE = 20;

// ---- 地面站点（可选，但强烈建议配置）----
// 有了它就可以直接在 GEE 里抽训练/验证表，完全不用导出栅格。
var STATIONS = ee.FeatureCollection([]); // 替换为上传的站点 Asset

// ============================================================
// 1. 数据准备
// ============================================================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(AOI)
  .filterDate(START, END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 85))
  .select(['B4', 'B8', 'B8A', 'B11', 'SCL']);

var cloudProb = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(AOI)
  .filterDate(START, END);

// ============================================================
// 2. 逐景计算指数与掩膜
// ============================================================
function addIndices(img) {
  var scl = img.select('SCL');

  // SCL 掩膜：剔除缺陷/云影/云/薄卷云/雪
  var valid = scl.neq(1).and(scl.neq(3)).and(scl.neq(8))
                 .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));

  // 用 s2cloudless 云概率补充（SCL 对薄卷云与云影漏检较多）
  var idx = img.get('system:index');
  var prob = ee.Image(cloudProb.filter(ee.Filter.eq('system:index', idx)).first())
               .select('probability');
  valid = valid.and(prob.lt(40)).rename('VALID');

  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndmi = img.normalizedDifference(['B8A', 'B11']).rename('NDMI');

  return ee.Image.cat([ndvi, ndmi, valid])
    .updateMask(valid)
    .copyProperties(img, ['system:time_start', 'system:index']);
}

var processed = s2.map(addIndices);

// ============================================================
// 3. 逐像元 NDVI 分位端元（用于 FVC，必须本地化，不要用默认 0.05/0.85）
// ============================================================
var ndviCol = processed.select('NDVI');
var percentiles = ndviCol.reduce(ee.Reducer.percentile([5, 95], ['p5', 'p95']));
var ndviSoil = percentiles.select('NDVI_p5');
var ndviVeg = percentiles.select('NDVI_p95');

function addFVC(img) {
  var ndvi = img.select('NDVI');
  var fvc = ndvi.subtract(ndviSoil).divide(ndviVeg.subtract(ndviSoil))
              .clamp(0, 1).rename('FVC');
  return img.addBands(fvc);
}

var final = processed.map(addFVC).select(['NDVI', 'NDMI', 'FVC', 'VALID']);

// ============================================================
// 4. 导出方式一（推荐）：站点抽表 —— 不导出栅格，直接拿训练/验证表
//    这是效率提升最大的一步：秒级完成，且不受导出配额限制。
// ============================================================
if (STATIONS.size().getInfo() > 0) {
  var list = final.toList(final.size());
  var table = ee.FeatureCollection(list.map(function (img) {
    img = ee.Image(img);
    var sampled = img.reduceRegions({
      collection: STATIONS,
      reducer: ee.Reducer.mean(),
      scale: SCALE,
      crs: TARGET_CRS,
      crsTransform: TARGET_TRANSFORM
    });
    return sampled.map(function (f) {
      return f.set('date', ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'))
              .set('image_id', img.get('system:index'));
    });
  }).flatten());

  Export.table.toDrive({
    collection: table,
    description: 'S2_station_samples',
    fileFormat: 'CSV'
  });
}

// ============================================================
// 5. 导出方式二：逐景导出栅格（本地做 S-1/S-2 融合时才需要）
//    注意：导出是瓶颈，几百景会触发任务排队；建议只导出与 S-1 过境匹配的日期。
// ============================================================
var dates = final.aggregate_array('system:time_start').getInfo();
var imageList = final.toList(final.size()).getInfo();

// 按需打开；一次性全导会占满任务队列
var DO_EXPORT = false;
if (DO_EXPORT) {
  imageList.forEach(function (imgInfo, i) {
    var img = ee.Image(imageList[i].id);
    var dateStr = ee.Date(dates[i]).format('YYYYMMdd').getInfo();
    Export.image.toDrive({
      image: img.toFloat(),
      description: 'S2_' + dateStr,
      folder: 'ecohydro_sm_s2',
      region: AOI,
      scale: SCALE,
      crs: TARGET_CRS,
      crsTransform: TARGET_TRANSFORM,
      maxPixels: 1e10,
      fileFormat: 'GeoTIFF'
    });
  });
}

// ============================================================
// 6. 快速质检（在 Console 里看，不用导出）
// ============================================================
var validRate = final.select('VALID').mean().multiply(100);
print('有效观测率均值 (%)', validRate.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: AOI, scale: 500, maxPixels: 1e10
}));
print('S-2 景数', final.size());
Map.centerObject(AOI, 9);
Map.addLayer(final.select('NDVI').median(), {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green']}, 'NDVI 中位');
Map.addLayer(validRate, {min: 20, max: 100, palette: ['red', 'yellow', 'green']}, '有效观测率');

// ============================================================
// 复现性提醒（写论文时必须记录）
// ============================================================
// GEE 数据集会随处理基线更新而变化。投稿前务必记录：
//   - ImageCollection ID 与访问日期
//   - 导出产品的版本号（建议把导出结果归档并给 DOI/永久链接）
//
// 硬边界：GEE **没有 Sentinel-1 SLC**，相位信息不可用，
// 因此 InSAR 相干性与形变支路（L2/L4）**不能**在 GEE 上做，
// 必须走本地 HyP3 + ISCE2/MintPy。
