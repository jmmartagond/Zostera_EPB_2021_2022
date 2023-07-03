//------ Definir polígono de cuerpo de agua------//

// Unir imagenes de 04 de Enero de 2022, usada para calcular el polígono de cobertura de agua en marea máxima
var pb20220104_Pl = ee.ImageCollection([
  pb20220104_Pl_1,
  pb20220104_Pl_2
  ]).mosaic();

// Map.addLayer(pb20220105_Pl, visParamsRGB, 'PB 20220105'); // Imagen en RGB
// Map.addLayer(pb20220619_Pl, visParamsRGB, 'PB 20220619'); // Imagen en RGB

// -- Método banda b8 (NIR)--//
// Map.addLayer(pb20220105_Pl.select('b8'), {min: 100, max: 8000 , gamma:3.0}, 'PB 20220105 B8', false); //Solo banda B8
// Map.addLayer(pb20220104_Pl.select('b4'), {min: 100, max: 8000 , gamma:3.0}, 'PB 20220104 B8', false);

//Definir umbral de la banda B8: 
// Cuando la banda 8 toma valores menores a 950, se hacen 1
var B8_thres = pb20220105_Pl.select(['b8'],['NIR']).lt(950);
// var B8_thres = pb20220104_Pl.select(['b4'],['NIR']).lt(950);
// Enmasacara todos los valores de la banda B8
B8_thres = B8_thres.updateMask(B8_thres.neq(0));

// Forma un conjunto de vectores (FeatureCollection) de los pixeles dentro del umbral de B8
var vectors_B8 = B8_thres.addBands(B8_thres.select('NIR')).reduceToVectors({
  geometry: crp_agua,
  // crs: B8_thres.projection(),
  scale: 3,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'zone',
  reducer: ee.Reducer.mean()
});

// Muestra los poligonos en el mapa
// var displaycB8 = ee.Image(0).updateMask(0).paint(vectors_B8, '000000', 1);
// Map.addLayer(displaycB8, {palette: '000000'}, 'Vectors B8', false);

// Funcion para calcular el area de cada uno de los poligonos cubiertos por agua
var Areas = function(feature) {
  return feature.set({Areas: feature.geometry().area(10)});
};

// Aplicar calculo de areas a los polígonos del umbral de B8
var area_de_poligonos_B8 = ee.Number(vectors_B8.map(Areas));

// Formar arreglo (Tabla/FeatureCollection) de los poligonos con el area calculada
var arr_areas_B8 = ee.FeatureCollection(area_de_poligonos_B8);

// Buscar el area maxima de todas las areas
var max_area_poly_B8 = arr_areas_B8.aggregate_max('Areas');
// print(max_area_poly_B8)

// Encontrar la feature con el area maxima
var max_poly_B8 = arr_areas_B8.filter(ee.Filter.eq("Areas", max_area_poly_B8));
exports.max_poly_B8=max_poly_B8;


// Mostrar el polígono con el area máxima en el mapa
var display_max_poly_B8 = ee.Image(0).updateMask(0).paint(max_poly_B8, '000000', 1);
// Map.addLayer(display_max_poly_B8, {palette: '40d61a'}, 'Poly B8');

Map.centerObject(max_poly_B8, 14);
exports.display_max_poly_B8 = display_max_poly_B8;

// Exportar el poligono a un KML
// Export.table.toDrive({
//   collection: max_poly_B8,
//   description:'ROI',
//   fileFormat: 'KML'
// });
