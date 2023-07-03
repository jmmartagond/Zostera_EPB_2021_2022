

// NDVI
var NDVI_ms = function(image) {
  var ndvi = image.normalizedDifference(['NIR', 'R']).rename('NDVI');
  return image.addBands(ndvi);
};

// EVI2
var EVI2_ms = function(image) {
  var NIR_ms = image.select('NIR'); 
  var RED_ms = image.select('R');
  var EVI2_ms = NIR_ms.subtract(RED_ms).divide(NIR_ms.add(RED_ms).add(1)).multiply(2.4).rename('EVI2');
  return image.addBands(EVI2_ms);
};

// Parámetros para visualizar EVI2
var evi2Vis = {
  bands:['EVI2'],
  min: -1,
  max: 1,
  palette: ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
          '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
          '012E01', '011D01', '011301']
};

// // Definir ROI donde el algoritmo va a realizar la clasificacion
// // Objeto para importar script de Planet con el cuerpo de agua
var imports_wtr_poly= require('users/martagond/Punta_Banda:PB_Planet_Polígono_de_Agua');

// // Clasificacion
// // Clasificacion No Supervisada
var classificationKmeans = function(image){
  var training = image.sample({
  region: pix_train,
  scale: 0.1,
  projection:'EPSG:32611',
  numPixels: 5000
  });
  var clusterer = ee.Clusterer.wekaKMeans(10).train(training);
  return image.cluster(clusterer);
};

// Clasificación supervisada
// Grounds Truth Points: Definir poligonos de pixeles donde se ha identificado las cubiertas
// Une todos los puntos de ground truth
var ground_t=zostera.merge(zostera_m).merge(spartina).merge(arena).merge(agua).merge(ulva);
// var ground_t=zostera_122021.merge(ulva_122021).merge(spartina_122021).merge(arena_122021).merge(agua_122021);


// Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
ground_t = ground_t.randomColumn();

var puntos_entrenamiento = ground_t.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion = ground_t.filter(ee.Filter.gte('random', 0.6));

var classificationRForest = function(image){
  var training = image.sampleRegions({
    collection: puntos_entrenamiento,
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
  tileScale: 2,
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
  // crs:'EPSG:32611',
  maxPixels: 1e12,
  bestEffort: true
  // tileScale: 16
  });
  var metros_area = ee.Number(area.get('classification'));
  var area_Km = metros_area.divide(1e6);
  return print('El area de la clase Zostera marina en km es', area_Km);
};


// Parámetros para colorear las clases
// var class_color =
// '<RasterSymbolizer>' +
//   '<ChannelSelection>' + //used when image has more than one band (to specify which band in which channel).
//     '<GrayChannel>' + 
//         '<SourceChannelName>1</SourceChannelName>' +
//     '</GrayChannel>' +
//   '</ChannelSelection>' +
//   '<ColorMap type="values">' +
//     '<ColorMapEntry color="#188205" quantity="0" />' +
//     '<ColorMapEntry color="#1a7004" quantity="1" />' +
//     '<ColorMapEntry color="#98ff00" quantity="2" />' +
//     '<ColorMapEntry color="#00ca22" quantity="3" />' +
//     '<ColorMapEntry color="#5cd9ff" quantity="4" />' +
//     '<ColorMapEntry color="#cebd9e" quantity="5" />' +
//     '<ColorMapEntry color="#ab6550" quantity="6" />' +
//     '<ColorMapEntry color="#0d38d6" quantity="7" />' +
//   '</ColorMap>' +
// '</RasterSymbolizer>';

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


// Unir las imágenes en RGB del 02/Diciembre/2021
// var RGB0212 = ee.ImageCollection([
//   RGB0212_1,
//   RGB0212_2
//   ]).mosaic();

// // Renombrar las bandas
// RGB0212 = RGB0212.select(['b1','b2','b3'],['R','G','B']);

// var RGB0212_clip = RGB0212.clip(imports_wtr_poly.max_poly_B8);


// Unir las imágenes en MS del 02/Diciembre/2021
var MS_0212 = ee.ImageCollection([
  MS0212_1,
  MS0212_2,
  MS0212_3,
  MS0212_4
  ]).mosaic();

// Renombrar las bandas
MS_0212 = MS_0212.select(['b1','b2','b3','b4'],['G','R','RE','NIR']);

// Calculo de indices
MS_0212 = NDVI_ms(MS_0212);
MS_0212 = EVI2_ms(MS_0212);

// Definir el poligono de marea maxima como objeto tipo ee.Geometry
var poly_max_B8 = imports_wtr_poly.max_poly_B8.geometry()

// Recortar la imagen al poligono
var MS_0212_clip = MS_0212.clip(poly_max_B8);
// var MS_0212_clip = MS_0212.clip(imports_wtr_poly.max_poly_B8);
// var MS_0212_clip = MS_0212.clip(poly_230422_PL);

// var area_image = MS_0212_clip.pixelArea();

// Calcular el area de la imagen recortada
// var area_image = ee.Image.pixelArea().mask(MS_0212_clip.select('G'));
// print(area_image)

// Map.addLayer(area_image, {}, 'Area');

// var unmaskedArea = area_image.reduceRegion({
//   reducer: ee.Reducer.sum(),
//   geometry: poly_max_B8,
//   scale: 0.22,
//   maxPixels: 1e13,
//   bestEffort: true
// }); 

// print('Area de la imagen recortada:', unmaskedArea);

// var unmaskedArea = area_image.reduceRegion({
//   reducer: ee.Reducer.sum(),
//   geometry: poly_max_B8,
//   scale: MS_0212.projection().nominalScale(),
//   maxPixels: 1e13,
//   bestEffort: true
// }); 
// print('Area of unmasked pixels:', unmaskedArea);


// Calcular área de imagen recortada
// var area_clip = MS_0212.reduceRegion({
//   reducer: ee.Reducer.sum(),
//   geometry: poly_max_B8,
//   scale: 0.2,
//   maxPixels: 1e13
// })
// print('Area de la imagen recortada', area_clip)

// // Calcular área de imagen recortada
// var area_clip = MS_0212_clip.reduceRegion({
//   reducer: ee.Reducer.sum(),
//   geometry: imports_wtr_poly.max_poly_B8,
//   scale: 0.5,
//   maxPixels: 1e13
// })
// print('Area de la imagen recortada', area_clip)

// var MS_0212_clip_NS_Kmeans = classificationKmeans(MS_0212_clip);
var MS_0212_clip_class_RForest = classificationRForest(MS_0212_clip);
var MS_0212_clip_class_RForest_layer = zostera_layer(MS_0212_clip_class_RForest);

var conf_mat_MS_0212_clip_class_RForest = conf_matriz(MS_0212_clip_class_RForest, puntos_validacion);

// var area_z_m_MS_0212_clip_class_RForest = area_t(MS_0212_clip_class_RForest);



// NIR0312_10 = NIR0312_10.select(['b1','b2','b3'],['NIR','R','G']);

// var NIR0312_clip = NIR0312_10.clip(imports_wtr_poly.max_poly_B8);


// print( 'El valor de overall accuracy de la clasificacion del 02/Diciembre/2021 es', conf_mat_MS_0212_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 02/Diciembre/2021 es', conf_mat_MS_0212_clip_class_RForest.kappa());
// print( 'El orden de las clases en la Matriz de oOnfusión del 02/Diciembre/2021 es', conf_mat_MS_0212_clip_class_RForest.getInfo());
// print( 'El valor de la Exactitud del  Productor de la clasificacion del 02/Diciembre/2021 es', conf_mat_MS_0212_clip_class_RForest.producersAccuracy());
// print( 'El valor de la Exactitud del  Usuario de la clasificacion del 02/Diciembre/2021 es', conf_mat_MS_0212_clip_class_RForest.consumersAccuracy());
// Map.addLayer(RGB0212_clip, {bands:['R','G','B']}, 'Sequoia 02/12/2022 RGB Cliped');
// Map.addLayer(NIR0312_clip, {bands:['NIR','R','G']}, 'Canon 03/12/2022 NIR Cliped');
// Map.addLayer(MS_0212_clip, {bands:['NIR','R','G'], min: -0.003, max: 0.15}, 'Sequoia 02/12/2022 MS Cliped');
// Map.addLayer(MS_0212_clip_NS_Kmeans.randomVisualizer(), {}, 'Clusters KMeans MS 02 12 2021', false);
// Map.addLayer(MS_0212_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest MS 02/12/2021', false);
// Map.addLayer(MS_0212_clip_class_RForest_layer, {palette: [ 'green']}, 'Parche de Zostera 02/Diciembre/2021');
exports.MS_0212_clip_class_RForest_layer = MS_0212_clip_class_RForest_layer;

// Export.image.toDrive({
//   image: MS_0212_clip_class_RForest_layer,
//   description: 'Capa_Zostera_marina_del_02Diciembre2021_Sequoia',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:0.5,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });
