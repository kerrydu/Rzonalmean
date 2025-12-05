// Google Earth Engine script: 计算湖南各市平均灯光亮度（与本地R/Python脚本等效）
// 使用方式：把本脚本粘贴到 https://code.earthengine.google.com/ 的 Code Editor 里运行
// 先将你的数据上传为资产（Assets）：
// 1) 右侧 Assets 面板 -> NEW -> Table upload（上传 hunan.shp 的压缩包 .zip）得到一个 FeatureCollection 资产
// 2) NEW -> Image upload（上传 dmsp-like2020.tif）得到一个 Image 资产
// 然后在下方设置 ASSET_IMAGE 与 ASSET_HUNAN

// ========== 用户需要设置的资产ID与参数 ==========
var ASSET_IMAGE = 'users/your_username/dmsp_like2020';  // 你的 TIF 资产ID
var ASSET_HUNAN = 'users/your_username/hunan';          // 你的 Hunan shapefile 资产ID（FeatureCollection）

var BAND_MODE = 'last';   // 可选: 数字N, 'first', 'last'(默认), 'avg'(各band像元平均为单层), 'all'(保留全部band)
var CITY_FIELD = null;    // 若知道城市字段名（如 'ShiName'），可填；为空时自动探测
var EXPORT_TO_DRIVE = true;  // 是否导出CSV到Google Drive
var EXPORT_FILE_NAME = 'hunan_light_mean_gee';
var EXPORT_FOLDER = null; // 指定Drive文件夹名称，null表示根目录
// ============================================

// 读取资产
var img = ee.Image(ASSET_IMAGE);
var hunan = ee.FeatureCollection(ASSET_HUNAN);

// 打印基础元信息（避免重计算，简单且稳定）
print('Image band names:', img.bandNames());
print('Image projection:', img.projection());
print('Image nominal scale (m):', img.projection().nominalScale());
print('Hunan feature count:', hunan.size());

// 自动识别城市字段
function detectCityField(fc) {
  var first = ee.Feature(fc.first());
  var keys = ee.List(first.toDictionary().keys());
  var candidates = ee.List(['ShiName', 'CITY', 'CITY_NAME', 'NAME', 'NAME_1', 'NAME_2',
                            'ADMIN_NAME', 'ADM2_NAME', 'COUNTY', 'PREF', 'PREFECTURE',
                            'CITY_CN', 'NAME_CH', 'NAME_ZH']);
  var matches = candidates.filter(function(k){ return keys.contains(k); });
  var blacklist = ee.List(['system:index']);
  var chosen = ee.Algorithms.If(
    ee.Number(matches.size()).gt(0),
    matches.get(0),
    keys.filter(function(k){ return blacklist.contains(k).not(); }).get(0)
  );
  return ee.String(chosen);
}

var cityField = CITY_FIELD ? ee.String(CITY_FIELD) : detectCityField(hunan);
print('Using city field:', cityField);

// 根据 BAND_MODE 选择/聚合 band
function selectByMode(image, mode) {
  var bandNames = image.bandNames();
  var nBands = bandNames.size();
  var m = (mode || 'last').toString().toLowerCase();

  // 是否为纯数字（1-based索引）
  var isDigits = m.match('^\\d+$');
  isDigits = isDigits ? true : false;

  if (isDigits) {
    var idx0 = ee.Number.parse(m).subtract(1);  // 转0基
    var bname = ee.String(bandNames.get(idx0));
    return image.select([bname]);
  }
  if (m === 'first') {
    var b0 = ee.String(bandNames.get(0));
    return image.select([b0]);
  }
  if (m === 'last') {
    var bl = ee.String(bandNames.get(nBands.subtract(1)));
    return image.select([bl]);
  }
  if (m === 'avg') {
    // 跨band求均值，得到单层Image，命名为 'mean'
    return image.reduce(ee.Reducer.mean()).rename('mean');
  }
  if (m === 'all') {
    return image; // 保留所有band
  }
  // 默认
  var bl2 = ee.String(bandNames.get(nBands.subtract(1)));
  return image.select([bl2]);
}

var selectedImg = selectByMode(img, BAND_MODE);
var scale = img.projection().nominalScale();

// 统一附加一个稳定的ID用于后续Join
var hunanWithId = hunan.map(function(f){ return f.set('uid', ee.String(f.id())); });

// 结果计算
var result;
if (BAND_MODE === 'all') {
  // 按band分别计算均值
  var perBand = selectedImg.reduceRegions({
    collection: hunanWithId,
    reducer: ee.Reducer.mean(),
    scale: scale,
    maxPixels: 1e13
  });
  // 同时再对band求像元平均得到单层，再计算均值，作为 mean_light
  var avgImg = img.reduce(ee.Reducer.mean()).rename('mean');
  var avgRes = avgImg.reduceRegions({
    collection: hunanWithId,
    reducer: ee.Reducer.mean(),
    scale: scale,
    maxPixels: 1e13
  });
  // 按uid Join，把 avgRes 的 'mean' 加为 mean_light
  var join = ee.Join.inner();
  var filt = ee.Filter.equals({leftField: 'uid', rightField: 'uid'});
  var joined = join.apply(perBand, avgRes, filt);
  result = ee.FeatureCollection(joined.map(function(pair){
    var left = ee.Feature(pair.get('primary'));
    var right = ee.Feature(pair.get('secondary'));
    return left.set('mean_light', right.get('mean'));
  }));
} else {
  // 单层：selectedImg只有1个band
  var single = selectedImg.reduceRegions({
    collection: hunanWithId,
    reducer: ee.Reducer.mean(),
    scale: scale,
    maxPixels: 1e13
  });
  var propName = ee.String(selectedImg.bandNames().get(0));
  // 统一输出为 mean_light
  result = single.map(function(f){
    return ee.Feature(f.geometry())
      .copyProperties(f)
      .set('mean_light', f.get(propName));
  });
}

// 只保留需要的字段：city, mean_light（若 all 模式则还会包含各band列）
var bandCols = selectedImg.bandNames();
var baseProps = ee.List(['uid', 'system:index']);
var cityColName = cityField;

// 添加 city 字段并选择输出列
var selectedResult = result.map(function(f){
  var cityName = f.get(cityColName);
  return f.set('city', cityName);
});

// 打印前10条
print('Preview (first 10):', selectedResult.limit(10));

// 导出CSV到Google Drive
if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection: selectedResult,
    description: EXPORT_FILE_NAME,
    fileNamePrefix: EXPORT_FILE_NAME,
    folder: EXPORT_FOLDER,
    fileFormat: 'CSV'
  });
}
