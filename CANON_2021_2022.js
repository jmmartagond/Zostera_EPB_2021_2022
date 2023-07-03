// Función para unir NIR y DSM y enmascarar los datos alpha=0
var nir_dsm = function(nir, dsm){
  return nir.addBands(dsm).updateMask(nir.select('b4').neq(0));
} ;

// NDVI CANON
var NDVI_cn = function(image) {
  var ndvi = image.normalizedDifference(['NIR', 'R']).rename('NDVI');
  return image.addBands(ndvi);
};

// EVI2 CANON c/nombres
var EVI2_cn = function(image) {
  var NIR_cn = image.select('NIR'); 
  var RED_cn = image.select('R');
  var EVI2_cn = NIR_cn.subtract(RED_cn).divide(NIR_cn.add(RED_cn).add(1)).multiply(2.4).rename('EVI2');
  return image.addBands(EVI2_cn);
};

var evi2Vis = {
  bands:['EVI2'],
  min: -1,
  max: 1,
  palette: ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
          '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
          '012E01', '011D01', '011301']
};

// Clasificacion

// Clasificacion No Supervisada
var classificationKmeans = function(image){
  var training = image.sample({
  region: pix_train,
  scale: 0.5,
  projection:'EPSG:4326',
  numPixels: 5000
  });
  var clusterer = ee.Clusterer.wekaKMeans(10).train(training);
  return image.cluster(clusterer);
};

// Clasificación supervisada
// Grounds Truth Points: Definir poligonos de pixeles donde se ha identificado las cubiertas
// var imports_PB= require('users/martagond/Punta_Banda:PB');


var classificationRForest = function(image,t_points){
  var training = image.sampleRegions({
    collection: t_points,
    properties: ['clase'],
    scale: 0.5
  });
  var trained = ee.Classifier.smileRandomForest(10).train({
    features: training, 
    classProperty: 'clase', 
    });
  return image.classify(trained);
};

// Matriz de confusión
var conf_matriz = function (image, v_points){
  var muestreo = image.sampleRegions({
  collection: v_points,
  properties: ['clase'],
  scale: 0.5,
  // projection:'EPSG:32611'
});
  return muestreo.errorMatrix('clase', 'classification');
};

// Extraer la capa clasificada como Zostera
var zostera_layer = function(image){
  var patch = image.eq(1); // Extraer la clase zostera, haciendo 1 donde se encuentra la clase zostera
  return image.updateMask(patch.eq(1)); // Enmascarar los lugares donde no hay zostera
};

// Calcular el area de la clase zostera
var area_t = function (image){
  var patch = image.eq(1);
  var area_pix = patch.multiply(ee.Image.pixelArea());
  var area = area_pix.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: imports_wtr_poly.max_poly_B8.geometry(),
  scale: 0.5,
  crs:'EPSG:32611',
  maxPixels: 1e12
  });
  var metros_area = ee.Number(area.get('classification'));
  var area_Km = metros_area.divide(1e6);
  return print('El area de la clase Zostera marina en km es', area_Km);
};

var class_color =
'<RasterSymbolizer>' +
  '<ChannelSelection>' + //used when image has more than one band (to specify which band in which channel).
    '<GrayChannel>' + 
        '<SourceChannelName>1</SourceChannelName>' +
    '</GrayChannel>' +
  '</ChannelSelection>' +
  '<ColorMap type="values">' +
    '<ColorMapEntry color="#188205" quantity="1" />' +
    '<ColorMapEntry color="#1a7004" quantity="2" />' +
    '<ColorMapEntry color="#98ff00" quantity="3" />' +
    '<ColorMapEntry color="#00ca22" quantity="4" />' +
    '<ColorMapEntry color="#5cd9ff" quantity="5" />' +
    '<ColorMapEntry color="#cebd9e" quantity="6" />' +
    '<ColorMapEntry color="#ab6550" quantity="7" />' +
    '<ColorMapEntry color="#0d38d6" quantity="8" />' +
    '<ColorMapEntry color="#ceff2e" quantity="9" />' +
    '<ColorMapEntry color="#d63000" quantity="10" />' +
  '</ColorMap>' +
'</RasterSymbolizer>';


// Objeto para importar script de polígono de cuerpo de agua
var imports_wtr_poly= require('users/martagond/Punta_Banda:PB_Planet_Polígono_de_Agua');



// 03 de Diciembre de 2021

//La imagen del 03 de Diciembre de 2021 incluye un segmento con algunos valores de 0, lo cual confunde demasiado
// a la clasificacion, arrojando un gran area como cierta clase, por ello se recorta a un poligono manualmente
NIR031221 = NIR031221.clip(cut_1221);
// NIR031221 = NIR031221.clip(canal_031221);
NIR031221 = NIR031221.select(['b1','b2','b3'],['NIR','R','G']); // Renombrar las bandas de NIR
NIR031221 = NDVI_cn(NIR031221); // Calculo de NDVI
NIR031221 = EVI2_cn(NIR031221); // Calculo de EVI2
var NIR031221_clip = NIR031221.clip(imports_wtr_poly.max_poly_B8);
// var NIR031221_clip_NS_Kmeans = classificationKmeans(NIR031221_clip);

// var ground_t_diciembre=zostera.merge(zostera_m).merge(spartina).merge(arena).merge(agua);
// ground_t_diciembre = ground_t_diciembre.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
// var puntos_entrenamiento_diciembre = ground_t_diciembre.filter(ee.Filter.lt('random', 0.6));
// var puntos_validacion_diciembre = ground_t_diciembre.filter(ee.Filter.gte('random', 0.6));

var ground_t_diciembre = zostera_122021.merge(ulva_122021).merge(spartina_122021).merge(agua_122021).merge(arena_122021);
ground_t_diciembre = ground_t_diciembre.randomColumn();
var puntos_entrenamiento_diciembre = ground_t_diciembre.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion_diciembre = ground_t_diciembre.filter(ee.Filter.gte('random', 0.6));
var NIR031221_clip_class_RForest = classificationRForest(NIR031221_clip, puntos_entrenamiento_diciembre);
var conf_mat_NIR031221_clip_RForest = conf_matriz(NIR031221_clip_class_RForest, puntos_validacion_diciembre);
// print( 'El valor de overall accuracy de la clasificacion del 03/Diciembre/2021 es', conf_mat_NIR031221_clip_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 03/Diciembre/2021 es', conf_mat_NIR031221_clip_RForest.kappa());
var NIR031221_clip_class_RForest_layer = zostera_layer(NIR031221_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_NIR031221 = area_t(NIR031221_clip_class_RForest); // Calcular el área de la clase zostera
// Map.addLayer(NIR031221, {}, 'Dron 03/12/2021 CANON 10 cm Georef', false);
// Map.addLayer(NIR031221_clip, {}, 'Dron 03/12/2021 CANON 10 cm Georef', false);
// Map.addLayer(NIR031221, {bands:'NDVI'}, 'Dron 03/12/2021 CANON 10 cm Georef NDVI ', false);
// Map.addLayer(NIR031221, evi2Vis, 'Dron 03/12/2021 CANON 10 cm Georef EVI2 ', false);
// Map.addLayer(NIR031221_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest NIR 03/12/2021');
  // Map.addLayer(NIR031221_clip_class_RForest_layer, {palette: [ 'green']}, 'Parche de Zostera 03/Diciembre/2022');
exports.NIR031221_clip_class_RForest_layer = NIR031221_clip_class_RForest_layer;

// NIR031221_clip_class_RForest_layer = NIR031221_clip_class_RForest_layer.toByte()

// Export.image.toDrive({
//   image: NIR031221_clip_class_RForest_layer,
//   description: 'Raster_Zostera_marina_03Diciembre2021_CANON',
//   folder:'Productos_Finales',
//   region:cut_1221,
//   scale:0.5,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });

// 19 y 20 de Mayo de 2022

// Aplicar funcion para unir cada tile de nir y dsm
NIR1905_1 = nir_dsm(NIR1905_1, NIR1905_DSM_1);
NIR1905_2 = nir_dsm(NIR1905_2, NIR1905_DSM_2);
NIR1905_3 = nir_dsm(NIR1905_3, NIR1905_DSM_3);
NIR1905_4 = nir_dsm(NIR1905_4, NIR1905_DSM_4);
NIR1905_5 = nir_dsm(NIR1905_5, NIR1905_DSM_5);
NIR1905_6 = nir_dsm(NIR1905_6, NIR1905_DSM_6);
NIR1905_7 = nir_dsm(NIR1905_7, NIR1905_DSM_7);
NIR1905_8 = nir_dsm(NIR1905_8, NIR1905_DSM_8);
NIR1905_9 = nir_dsm(NIR1905_9, NIR1905_DSM_9);
NIR1905_10 = nir_dsm(NIR1905_10, NIR1905_DSM_10);
NIR2005 = nir_dsm(NIR2005, NIR2005_DSM);

var nir_0522 = ee.ImageCollection([
  NIR1905_1,
  NIR1905_2,
  NIR1905_3,
  NIR1905_4,
  NIR1905_5,
  NIR1905_6,
  NIR1905_7,
  NIR1905_8,
  NIR1905_9,
  NIR1905_10,
  NIR2005
  ]).mosaic();



nir_0522 = nir_0522.select(['b1','b2','b3','b4','b1_1'],['NIR','R','G','Alpha','DSM']); //Renombrar las bandas
nir_0522 = NDVI_cn(nir_0522);
nir_0522 = EVI2_cn(nir_0522);
var nir_0522_clip = nir_0522.clip(imports_wtr_poly.max_poly_B8);
// var nir_0522_clip = nir_0522.clip(poly_max_0522_mod);
var nir_0522_clip_NS_Kmeans = classificationKmeans(nir_0522_clip);
var ground_t=zostera.merge(zostera_m).merge(spartina).merge(arena).merge(agua);
ground_t = ground_t.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
var puntos_entrenamiento = ground_t.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion = ground_t.filter(ee.Filter.gte('random', 0.6));
var nir_0522_clip_class_RForest = classificationRForest(nir_0522_clip, puntos_entrenamiento);
var conf_mat_nir_0522_clip_class_RForest = conf_matriz(nir_0522_clip_class_RForest, puntos_validacion);
// print( 'El valor de overall accuracy de la clasificacion del 19-20/Mayo/2022 es', conf_mat_nir_0522_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 19-20/Mayo/2022 es', conf_mat_nir_0522_clip_class_RForest.kappa());
var nir_0522_clip_class_RForest_layer = zostera_layer(nir_0522_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_nir_0522 = area_t(nir_0522_clip_class_RForest); // Calcular el área de la clase zostera
// Map.addLayer(nir_0522, {bands:['NIR','R','G']}, 'Dron 05/2022 CANON NIR');
// Map.addLayer(nir_0522_clip, {bands:['NIR','R','G']}, 'Dron 05/2022 CANON NIR Cliped');
// Map.addLayer(nir_0522, {bands:'NDVI'}, 'Dron 05/2022 CANON 10 cm Georef NDVI ', false);
// Map.addLayer(nir_0522, evi2Vis, 'Dron 05/2022 CANON 10 cm Georef EVI2 ', false);
// Map.addLayer(nir_0522_clip_NS_Kmeans.randomVisualizer(), {}, 'Clusters KMeans NIR 0522', false);
// Map.addLayer(nir_0522_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest NIR 05 2022');
// Map.addLayer(nir_0522_clip_class_RForest_layer, {palette: [ 'green']}, 'Parche de Zostera Mayo/2022');
exports.nir_0522_clip_class_RForest_layer = nir_0522_clip_class_RForest_layer;

// nir_0522_clip_class_RForest_layer = nir_0522_clip_class_RForest_layer.toByte();

// Export.image.toDrive({
//   image: nir_0522_clip_class_RForest_layer,
//   description: 'Raster_Zostera_marina_Mayo2022_CANON',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:0.5,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });
