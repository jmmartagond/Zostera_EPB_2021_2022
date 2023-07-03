
// Calculo del indice NDVI para Planet 4-Bandas
var NDVI_4 = function(image) {
  var ndvi = image.normalizedDifference(['b4', 'b3']).rename('NDVI');
return image.addBands(ndvi);
};

// Calculo del indice NDVI para Planet 8-Bandas
var NDVI = function(image) {
  var ndvi = image.normalizedDifference(['b8', 'b6']).rename('NDVI');
return image.addBands(ndvi);
};

// Calculo del indice EVI2 para Planet 4-Bandas
var EVI2_4 = function(image) {
  var NIR = image.select('b4');
  var RED = image.select('b3');
  var EVI2 = NIR.subtract(RED).divide(NIR.add(RED).add(1)).multiply(2.4).rename('EVI2');
return image.addBands(EVI2);
};

// Calculo del indice EVI2 para Planet 8-Bandas
var EVI2 = function(image) {
  var NIR = image.select('b8');
  var RED = image.select('b6');
  var EVI2 = NIR.subtract(RED).divide(NIR.add(RED).add(1)).multiply(2.4).rename('EVI2');
return image.addBands(EVI2);
};

// Importar el polígono para realizar el recorte de la Zona de Interés
var imports_wtr_poly= require('users/martagond/Punta_Banda:PB_Planet_Polígono_de_Agua');

// Clasificación
// Clasificación No Supervisada
var clus_kMeans = function(image){
  var training = image.sample({
    region: pix_train,
    scale: 3,
    // projection:'EPSG:32611',
    numPixels: 5000});
  var clusterer = ee.Clusterer.wekaKMeans(10).train(training);
  return image.cluster(clusterer);
};

// Clasificación Supervisada
// Función de entrenamiento y clasificación
var classificationRForest = function(image,t_points){
  var training = image.sampleRegions({
    collection: t_points,
    properties: ['clase'],
    scale: 3,
    // projection:'EPSG:32611'
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
  scale: 3,
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
  scale: 3,
  crs:'EPSG:32611'
  // maxPixels: 1e12
  });
  var metros_area = ee.Number(area.get('classification'));
  var area_Km = metros_area.divide(1e6);
  return print('El area de la clase Zostera marina en km es', area_Km);
};

// Parámetros de visualización en RGB 8 bandas
var visParamsRGB_8 = {bands: ['b6', 'b4', 'b2'], min: 100, max: 8000 ,gamma:3.0};

// Parámetros de visualización en Falso Infrarrojo 8 bandas
var visParamsInfra_8 = {bands: ['b8', 'b6', 'b4'], min: 100, max: 8000 ,gamma:3.0};

// Parámetros para visualizar EVI2
var evi2Vis = {
  bands:['EVI2'],
  min: -1,
  max: 1,
  palette: ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
          '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
          '012E01', '011D01', '011301']
};

// Parámetros para colorear las clases
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
  '</ColorMap>' +
'</RasterSymbolizer>';


// 31 de Octubre de 2021
// Unir imagenes de 31 de Octubre de 2021
var pb20211031_Pl = ee.ImageCollection([
  pb20211031_Pl_1,
  pb20211031_Pl_2
  ]).mosaic();
  
pb20211031_Pl= NDVI_4(pb20211031_Pl);
pb20211031_Pl = EVI2_4(pb20211031_Pl);
var pb20211031_Pl_clip = pb20211031_Pl.clip(imports_wtr_poly.max_poly_B8);
var pb20211031_Pl_clip_NS_Kmeans = clus_kMeans(pb20211031_Pl_clip);
var ground_t=zostera.merge(spartina).merge(arena).merge(agua);// Une todos los puntos de ground truth
ground_t = ground_t.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
var puntos_entrenamiento = ground_t.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion = ground_t.filter(ee.Filter.gte('random', 0.6));
var pb20211031_Pl_clip_class_RForest = classificationRForest(pb20211031_Pl_clip, puntos_entrenamiento);
var pb20211031_Pl_clip_class_RForest_layer = zostera_layer(pb20211031_Pl_clip_class_RForest);
// var area_z_m_pb20211031_Pl = area_t(pb20211031_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20211031_Pl_clip_class_RForest = conf_matriz(pb20211031_Pl_clip_class_RForest, puntos_validacion); //Calcula la matriz de confusión
// print( 'El valor de overall accuracy de la clasificacion del 31/Octubre/2021 es', conf_mat_pb20211031_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 31/Octubre/2021 es', conf_mat_pb20211031_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20211031_Pl, {bands: ['b3', 'b2', 'b1'], min: 100, max: 8000 ,gamma:3.0}, 'RGB 31/Octubre/2021 Planet', false);
// Map.addLayer(pb20211031_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'KMeans 31/Octubre/2021  Planet', false);
// Map.addLayer(pb20211031_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 31/Octubre/2021 Planet');
// Map.addLayer(pb20211031_Pl_clip_class_RForest_layer, {palette: [ 'green']}, 'Parche de Zostera 31/Octubre/2021');
exports.pb20211031_Pl_clip_class_RForest_layer = pb20211031_Pl_clip_class_RForest_layer;

// Filtro Mayoria
// var pb20211031_Pl_clip_class_RForest_majorityf = pb20211031_Pl_clip_class_RForest.reduceNeighborhood({
//   reducer:ee.Reducer.mode(),
//   kernel:ee.Kernel.square(5)
// });
// Map.addLayer(pb20211031_Pl_clip_class_RForest_majorityf.sldStyle(class_color), {}, 'R F 31/Octubre/2021 Planet_Majority');

// pb20211031_Pl_clip_class_RForest_layer = pb20211031_Pl_clip_class_RForest_layer.toByte();

// Export.image.toDrive({
//   image: pb20211031_Pl_clip_class_RForest_layer,
//   description: 'Raster_Zostera_marina_31Octubre2021_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });


// 26 de Enero de 2022
// Unir imagenes de 26 de Enero de 2022
var pb20220126_Pl = ee.ImageCollection([
  pb20220126_Pl_1,
  pb20220126_Pl_2
  ]).mosaic();
  
pb20220126_Pl =NDVI_4(pb20220126_Pl); //Calcular Índice NDVI
pb20220126_Pl = EVI2_4(pb20220126_Pl); // Calcular índice EVI2
var pb20220126_Pl_clip = pb20220126_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
var pb20220126_Pl_clip_NS_Kmeans = clus_kMeans(pb20220126_Pl_clip); //Realiza la clasificación no supervisada de kMeans
var ground_t_enero=zostera_enero.merge(ulva_enero).merge(spartina).merge(arena).merge(agua);// Une todos los puntos de ground truth
ground_t_enero = ground_t_enero.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
var puntos_entrenamiento_enero = ground_t_enero.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion_enero = ground_t_enero.filter(ee.Filter.gte('random', 0.6));
var pb20220126_Pl_clip_class_RForest = classificationRForest(pb20220126_Pl_clip, puntos_entrenamiento_enero); //Realiza la clasificación Supervisada
var pb20220126_Pl_clip_class_RForest_layer = zostera_layer(pb20220126_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20220126_Pl = area_t(pb20220126_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20220126_Pl_clip_class_RForest = conf_matriz(pb20220126_Pl_clip_class_RForest, puntos_validacion_enero); //Calcula la matriz de confusión
// print( 'El valor de overall accuracy de la clasificacion del 26/Enero/2022 es', conf_mat_pb20220126_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 26/Enero/2022 es', conf_mat_pb20220126_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20220126_Pl, {bands: ['b3', 'b2', 'b1'], min: 100, max: 8000 ,gamma:3.0}, 'RGB 26/Enero/2022 Planet');
// Map.addLayer(pb20220126_Pl, {bands:'NDVI'}, 'NDVI 26/Enero/2022', false);
// Map.addLayer(pb20220126_Pl, evi2Vis, 'EVI2 26/Enero/2022', false);
// Map.addLayer(pb20220126_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'KMeans 26/Enero/2022 Planet', false);
// Map.addLayer(pb20220126_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 26/Enero/2022 Planet');
// Map.addLayer(pb20220126_Pl_clip_class_RForest_layer, {palette: [ 'green']}, 'Parche de Zostera 26/Enero/2022');
exports.pb20220126_Pl_clip_class_RForest_layer = pb20220126_Pl_clip_class_RForest_layer;


//Cambia la proyección de la Clasificacion
// pb20220126_Pl_clip_class_RForest=pb20220126_Pl_clip_class_RForest.reproject({
//   crs: 'EPSG:32611',
//   scale: 3
// });

// Export.image.toDrive({
//   image: pb20220126_Pl_clip_class_RForest_layer,
//   description: 'Capa_Zostera_marina_del_26Enero2022_Planet',
//   folder:'Tesis',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });

// 26 de Febrero de 2022
pb20220226_Pl = NDVI(pb20220226_Pl); //Calcular Índice NDVI
pb20220226_Pl = EVI2(pb20220226_Pl); // Calcular índice EVI2
var pb20220226_Pl_clip = pb20220226_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
var pb20220226_Pl_clip_NS_Kmeans = clus_kMeans(pb20220226_Pl_clip); //Realiza la clasificación no supervisada de kMeans
var pb20220226_Pl_clip_class_RForest = classificationRForest(pb20220226_Pl_clip, puntos_entrenamiento); //Realiza la clasificación Supervisada
var pb20220226_Pl_clip_class_RForest_layer = zostera_layer(pb20220226_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20220226_Pl = area_t(pb20220226_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20220226_Pl_clip_class_RForest = conf_matriz(pb20220226_Pl_clip_class_RForest, puntos_validacion);
// print( 'El valor de overall accuracy de la clasificacion del 26/Febrero/2022 es', conf_mat_pb20220226_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 26/Febrero/2022 es', conf_mat_pb20220226_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20220226_Pl, {bands: ['b3', 'b2', 'b1'], min: 100, max: 8000 ,gamma:3.0}, 'RGB 26/Febrero/2022 Planet', false);
// Map.addLayer(pb20220226_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'KMeans 26/Febrero/2022 Planet', false);
// Map.addLayer(pb20220226_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 26/Febrero/2022 Planet');
// Map.addLayer(pb20220226_Pl_clip_class_RForest_layer, {palette: [ 'green']}, 'Parche de Zostera 26/Febrero/2022');
exports.pb20220226_Pl_clip_class_RForest_layer = pb20220226_Pl_clip_class_RForest_layer;

// Export.image.toDrive({
//   image: pb20220226_Pl_clip_class_RForest_layer,
//   description: 'Capa_Zostera_marina_del_26Febrero2022_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });

// 23 de Abril de 2022
pb20220423_Pl = NDVI(pb20220423_Pl); //Calcular Índice NDVI
pb20220423_Pl = EVI2(pb20220423_Pl); // Calcular índice EVI2
// var pb20220423_Pl_clip = pb20220423_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
var pb20220423_Pl_clip = pb20220423_Pl.clip(poly_230422_PL); // Recorte de la Zona de Interés
var pb20220423_Pl_clip_NS_Kmeans = clus_kMeans(pb20220423_Pl_clip); //Realiza la clasificación no supervisada de kMeans
var ground_t_mayo=zostera_mayo.merge(ulva_mayo).merge(spartina).merge(arena).merge(agua); // Une todos los puntos de ground truth
ground_t_mayo = ground_t_mayo.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
var puntos_entrenamiento_mayo = ground_t_mayo.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion_mayo = ground_t_mayo.filter(ee.Filter.gte('random', 0.6));
var pb20220423_Pl_clip_class_RForest = classificationRForest(pb20220423_Pl_clip, puntos_entrenamiento_mayo); //Realiza la clasificación Supervisada
var pb20220423_Pl_clip_class_RForest_layer = zostera_layer(pb20220423_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20220423_Pl = area_t(pb20220423_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20220423_Pl_clip_class_RForest = conf_matriz(pb20220423_Pl_clip_class_RForest, puntos_validacion_mayo);
// print( 'El valor de overall accuracy de la clasificacion del 23/Abril/2022 es', conf_mat_pb20220423_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 23/Abril/2022 es', conf_mat_pb20220423_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20220423_Pl, visParamsRGB_8, 'RGB 23/Abril/2022', false);
// Map.addLayer(pb20220423_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'K Means 23/Abril/2022', false);
// Map.addLayer(pb20220423_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 23/Abril/2022', false);
// Map.addLayer(pb20220423_Pl_clip_class_RForest_layer ,{palette: [ 'green']}, 'Parche de Zostera 23 de Abril de 2022');
exports.pb20220423_Pl_clip_class_RForest_layer = pb20220423_Pl_clip_class_RForest_layer; 
 
// Export.image.toDrive({
//   image: pb20220423_Pl_clip_class_RForest_layer,
//   description: 'Capa_Zostera_marina_del_23Abril2022_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });

 
// 06 de Junio de 2022
var pb20220606_Pl = NDVI(pb20220606_Pl_1); //Calcular Índice NDVI
pb20220606_Pl= EVI2(pb20220606_Pl); // Calcular índice EVI2
pb20220606_Pl_1=NDVI(pb20220606_Pl_1)
pb20220606_Pl_1=EVI2(pb20220606_Pl_1)
pb20220606_Pl_2=NDVI(pb20220606_Pl_2)
pb20220606_Pl_2=EVI2(pb20220606_Pl_2)
var pb20220606_Pl_clip = pb20220606_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
var pb20220606_Pl_clip_NS_Kmeans = clus_kMeans(pb20220606_Pl_clip); //Realiza la clasificación no supervisada de kMeans
var ground_t_junio=zostera_junio.merge(ulva_junio).merge(spartina).merge(arena).merge(agua);// Une todos los puntos de ground truth
ground_t_junio = ground_t.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
var puntos_entrenamiento_junio = ground_t_junio.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion_junio = ground_t_junio.filter(ee.Filter.gte('random', 0.6));
var pb20220606_Pl_clip_class_RForest = classificationRForest(pb20220606_Pl_clip, puntos_entrenamiento_mayo); //Realiza la clasificación supervisada
var pb20220606_Pl_clip_class_RForest_layer = zostera_layer(pb20220606_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20220606_Pl = area_t(pb20220606_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20220606_Pl_clip_class_RForest =conf_matriz(pb20220606_Pl_clip_class_RForest, puntos_validacion_mayo);
// print( 'El valor de overall accuracy de la clasificacion del 06/Junio/2022 es', conf_mat_pb20220606_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 06/Junio/2022 es', conf_mat_pb20220606_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20220606_Pl_1, visParamsRGB_8, 'RGB 06/Junio/2022_1', false);
// Map.addLayer(pb20220606_Pl_1, {"opacity":1,"bands":["b8"],"min":493.2,"max":2622.8,"gamma":1} , 'B8 06/Junio/2022_1', false);
// Map.addLayer(pb20220606_Pl_1, {"opacity":1,"bands":["EVI2"],"min":-0.6373138961306357,"max":2.015847522162645,"gamma":1} , 'EVI2 06/Junio/2022_1', false);
// Map.addLayer(pb20220606_Pl_2, visParamsRGB_8, 'RGB 06/Junio/2022_2', false);
// Map.addLayer(pb20220606_Pl_2, {"opacity":1,"bands":["b8"],"min":585.7,"max":2975.3,"gamma":1} , 'B8 06/Junio/2022_2', false);
// Map.addLayer(pb20220606_Pl_2, {"opacity":1,"bands":["EVI2"],"min":-0.6373138961306357,"max":2.015847522162645,"gamma":1} , 'EVI2 06/Junio/2022_2', false);
// Map.addLayer(pb20220606_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'K Means 06/Junio/2022', false);
// Map.addLayer(pb20220606_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 06/Junio/2022', false);
// Map.addLayer(pb20220606_Pl_clip_class_RForest_layer ,{palette: [ 'green']}, 'Parche de Zostera 06 de Junio de 2022');
exports.pb20220606_Pl_clip_class_RForest_layer = pb20220606_Pl_clip_class_RForest_layer; 

// Export.image.toDrive({
//   image: pb20220606_Pl_clip_class_RForest_layer,
//   description: 'Capa_Zostera_marina_del_06Junio2022_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });

//19 de Junio de 2022
pb20220619_Pl = NDVI(pb20220619_Pl); //Calcular Índice NDVI
pb20220619_Pl = EVI2(pb20220619_Pl); // Calcular índice EVI2
var pb20220619_Pl_clip = pb20220619_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
var pb20220619_Pl_clip_NS_Kmeans = clus_kMeans(pb20220619_Pl_clip); //Realiza la clasificación no supervisada de kMeans
var pb20220619_Pl_clip_class_RForest = classificationRForest(pb20220619_Pl_clip, puntos_entrenamiento_mayo); //Realiza la clasificación supervisada
var pb20220619_Pl_clip_class_RForest_layer = zostera_layer(pb20220619_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20220619_Pl = area_t(pb20220619_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20220619_Pl_clip_class_RForest = conf_matriz(pb20220619_Pl_clip_class_RForest, puntos_validacion_mayo);
// print( 'El valor de overall accuracy de la clasificacion del 19/Junio/2022 es', conf_mat_pb20220619_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 19/Junio/2022 es', conf_mat_pb20220619_Pl_clip_class_RForest.kappa());
// print( 'El valor de producers accuracy de la clasificacion del 19/Junio/2022 es', conf_mat_pb20220619_Pl_clip_class_RForest.producersAccuracy());
// Map.addLayer(pb20220619_Pl, visParamsRGB_8, 'RGB 19/Junio/2022', false);
// Map.addLayer(pb20220619_Pl, {"opacity":1,"bands":["b8"],"min":495.40000000000003,"max":3298.6,"gamma":1}, 'B8 19/Junio/2022', false);
// Map.addLayer(pb20220619_Pl, {"opacity":1,"bands":["EVI2"],"min":-0.43868286727173766,"max":1.7384783983760168,"gamma":1}, 'EVI2 19/Junio/2022', false);
// Map.addLayer(pb20220619_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'KMeans 19/Junio/2022', false);
// Map.addLayer(pb20220619_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 19/Junio/2022', false);
// Map.addLayer(pb20220619_Pl_clip_class_RForest_layer ,{palette: [ 'green']}, 'Parche de Zostera 19 de Junio de 2022');

// Filtro Mayoria
// var pb20220619_Pl_clip_class_RForest_majorityf = pb20220619_Pl_clip_class_RForest.reduceNeighborhood({
//   reducer:ee.Reducer.mode(),
//   kernel:ee.Kernel.square(3)
// });
// Map.addLayer(pb20220619_Pl_clip_class_RForest_majorityf.sldStyle(class_color), {}, 'R F 19/Junio/2022 Planet_Majority');


// pb20220619_Pl_clip_class_RForest_layer = pb20220619_Pl_clip_class_RForest_layer.toByte();
// Export.image.toDrive({
//   image: pb20220619_Pl_clip_class_RForest_layer,
//   description: 'Raster_Zostera_marina_19Junio2022_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// });


//15 de Febrero de 2023
// pb20230215_Pl = NDVI(pb20230215_Pl); //Calcular Índice NDVI
// pb20230215_Pl = EVI2(pb20230215_Pl); // Calcular índice EVI2
// var pb20230215_Pl_clip = pb20230215_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
// var pb20230215_Pl_clip_NS_Kmeans = clus_kMeans(pb20230215_Pl_clip); //Realiza la clasificación no supervisada de kMeans
// var pb20230215_Pl_clip_class_RForest = classificationRForest(pb20230215_Pl_clip, puntos_entrenamiento_mayo); //Realiza la clasificación supervisada
// var pb20230215_Pl_clip_class_RForest_layer = zostera_layer(pb20230215_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20230216_Pl = area_t(pb20230216_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20230215_Pl_clip_class_RForest = conf_matriz(pb20230215_Pl_clip_class_RForest, puntos_validacion_mayo);
// print( 'El valor de overall accuracy de la clasificacion del 15/Febrero/2023 es', conf_mat_pb20230215_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 15/Febrero/2023 es', conf_mat_pb20230215_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20230215_Pl, visParamsRGB_8, 'RGB 15/Febrero/2023', false);
// Map.addLayer(pb20230215_Pl, {"opacity":1,"bands":["b8"],"min":495.40000000000003,"max":3298.6,"gamma":1}, 'B8 15/Febrero/2023', false);
// Map.addLayer(pb20230215_Pl, {"opacity":1,"bands":["EVI2"],"min":-0.43868286727173766,"max":1.7384783983760168,"gamma":1}, 'EVI2 15/Febrero/2023', false);
// Map.addLayer(pb20230215_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'KMeans 15/Febrero/2023', false);
// Map.addLayer(pb20230215_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 15/Febrero/2023', false);
// Map.addLayer(pb202230215_Pl_clip_class_RForest_layer ,{palette: [ 'green']}, 'Parche de Zostera 15/Febrero/2023');

//16 de Febrero de 2023
// pb20230216_Pl = NDVI(pb20230216_Pl); //Calcular Índice NDVI
// pb20230216_Pl = EVI2(pb20230216_Pl); // Calcular índice EVI2
// var pb20230216_Pl_clip = pb20230216_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
// var pb20230216_Pl_clip_NS_Kmeans = clus_kMeans(pb20230216_Pl_clip); //Realiza la clasificación no supervisada de kMeans
// var pb20230216_Pl_clip_class_RForest = classificationRForest(pb20230216_Pl_clip, puntos_entrenamiento_mayo); //Realiza la clasificación supervisada
// var pb20230216_Pl_clip_class_RForest_layer = zostera_layer(pb20230216_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20230216_Pl = area_t(pb20230216_Pl_clip_class_RForest); // Calcular el área de la clase zostera
// var conf_mat_pb20230216_Pl_clip_class_RForest = conf_matriz(pb20230216_Pl_clip_class_RForest, puntos_validacion_mayo);
// print( 'El valor de overall accuracy de la clasificacion del 16/Febrero/2023 es', conf_mat_pb20230216_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 16/Febrero/2023 es', conf_mat_pb20230216_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20230216_Pl, visParamsRGB_8, 'RGB 16/Febrero/2023', false);
// Map.addLayer(pb20230216_Pl, {"opacity":1,"bands":["b8"],"min":495.40000000000003,"max":3298.6,"gamma":1}, 'B8 16/Febrero/2023', false);
// Map.addLayer(pb20230216_Pl, {"opacity":1,"bands":["EVI2"],"min":-0.43868286727173766,"max":1.7384783983760168,"gamma":1}, 'EVI2 16/Febrero/2023', false);
// Map.addLayer(pb20230216_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'KMeans 16/Febrero/2023', false);
// Map.addLayer(pb20230216_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 16/Febrero/2023', false);
// Map.addLayer(pb202230216_Pl_clip_class_RForest_layer ,{palette: [ 'green']}, 'Parche de Zostera 16/Febrero/2023');


// 29 de Marzo de 2023
pb20230329_Pl = NDVI(pb20230329_Pl); //Calcular Índice NDVI
pb20230329_Pl = EVI2(pb20230329_Pl); // Calcular índice EVI2
var pb20230329_Pl_clip = pb20230329_Pl.clip(imports_wtr_poly.max_poly_B8); // Recorte de la Zona de Interés
// var pb20230329_Pl_clip = pb20230329_Pl.clip(poly_230422_PL); // Recorte de la Zona de Interés
var pb20230329_Pl_clip_NS_Kmeans = clus_kMeans(pb20230329_Pl_clip); //Realiza la clasificación no supervisada de kMeans
var ground_t_mayo=zostera_mayo.merge(ulva_mayo).merge(spartina).merge(arena).merge(agua); // Une todos los puntos de ground truth
ground_t_mayo = ground_t_mayo.randomColumn(); // Agregar una columna aleatoria y dividir los puntos en dos conjuntos de entrenamiento y validación
var puntos_entrenamiento_mayo = ground_t_mayo.filter(ee.Filter.lt('random', 0.6));
var puntos_validacion_mayo = ground_t_mayo.filter(ee.Filter.gte('random', 0.6));
var pb20230329_Pl_clip_class_RForest = classificationRForest(pb20230329_Pl_clip, puntos_entrenamiento_mayo); //Realiza la clasificación Supervisada
var pb20230329_Pl_clip_class_RForest_layer = zostera_layer(pb20230329_Pl_clip_class_RForest); //Crea una capa de los pixeles clasificados como zostera
// var area_z_m_pb20230329_Pl = area_t(pb20230329_Pl_clip_class_RForest); // Calcular el área de la clase zostera
var conf_mat_pb20230329_Pl_clip_class_RForest = conf_matriz(pb20230329_Pl_clip_class_RForest, puntos_validacion_mayo);
// print( 'El valor de overall accuracy de la clasificacion del 29/Marzo/2023 es', conf_mat_pb20230329_Pl_clip_class_RForest.accuracy());
// print( 'El valor kappa de la clasificacion del 29/Marzo/2023 es', conf_mat_pb20230329_Pl_clip_class_RForest.kappa());
// Map.addLayer(pb20230329_Pl, visParamsRGB_8, 'RGB 29/Marzo/2023', false);
// Map.addLayer(pb20230329_Pl, {bands:'NDVI'}, 'NDVI 29/Marzo/2023', false);
// Map.addLayer(pb20230329_Pl, evi2Vis, 'EVI2 29/Marzo/2023', false);
// Map.addLayer(pb20230329_Pl_clip_NS_Kmeans.randomVisualizer(), {}, 'K Means 29/Marzo/2023', false);
// Map.addLayer(pb20230329_Pl_clip_class_RForest.sldStyle(class_color), {}, 'Random Forest 29/Marzo/2023', false);
// Map.addLayer(pb20230329_Pl_clip_class_RForest_layer ,{palette: [ 'green']}, 'Parche de Zostera 29/Marzo/2023');
exports.pb20230329_Pl_clip_class_RForest_layer = pb20230329_Pl_clip_class_RForest_layer; 


// Mapa de Probabilidad

// 2022
// Coleccion de capas aisladas de zostera de 2022
var zostera_2022_Pl = ee.ImageCollection([
  pb20220126_Pl_clip_class_RForest_layer,
  pb20220226_Pl_clip_class_RForest_layer,
  pb20220423_Pl_clip_class_RForest_layer,
  pb20220606_Pl_clip_class_RForest_layer,
  pb20220619_Pl_clip_class_RForest_layer
  ]);
  
// Sumar la clase zostera
var suma_zostera_2022_Pl = zostera_2022_Pl.sum();


// Map.addLayer(suma_zostera_2022_Pl ,{palette: [ 'green']}, 'Suma de Zostera 2022');
exports.suma_zostera_2022_Pl = suma_zostera_2022_Pl;

// suma_zostera_2022_Pl = suma_zostera_2022_Pl.toByte();
// Export.image.toDrive({
//   image: suma_zostera_2022_Pl,
//   description: 'Probabilidad_Zostera_2022_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// }); 


// Coleccion completa de Planet

var zostera_total_Pl = ee.ImageCollection([
  pb20211031_Pl_clip_class_RForest_layer,
  pb20220126_Pl_clip_class_RForest_layer,
  pb20220226_Pl_clip_class_RForest_layer,
  pb20220423_Pl_clip_class_RForest_layer,
  pb20220606_Pl_clip_class_RForest_layer,
  pb20220619_Pl_clip_class_RForest_layer
  ]);

var suma_zostera_t_Pl = zostera_total_Pl.sum();
// Map.addLayer(suma_zostera_t_Pl ,{palette: [ 'green']}, 'Suma de Zostera Total');
exports.suma_zostera_t_Pl = suma_zostera_t_Pl;

// suma_zostera_t_Pl = suma_zostera_t_Pl.toByte();
// Export.image.toDrive({
//   image: suma_zostera_t_Pl,
//   description: 'Probabilidad_Zostera_Total_Planet',
//   folder:'Productos_Finales',
//   region:imports_wtr_poly.max_poly_B8,
//   scale:3,
//   crs:'EPSG:32611',
//   maxPixels:1e12
// }); 
