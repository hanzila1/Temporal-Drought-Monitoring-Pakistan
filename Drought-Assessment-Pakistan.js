// Pakistan Drought Monitoring using Vegetation Condition Index (VCI)
// This script analyzes drought conditions in Pakistan by calculating VCI from MODIS NDVI data
// and overlays it with cropland data to identify agricultural drought.

// Define Pakistan boundary using a rectangle and then get actual boundary from a dataset
var pakistanBounds = ee.Geometry.Rectangle([60.5, 23.5, 77.5, 37.5]);
var countries = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");
var pakistan = countries.filter(ee.Filter.eq('country_na', 'Pakistan'));

// Define data sources
var modis = ee.ImageCollection('MODIS/061/MOD13Q1');
var gcep30 = ee.ImageCollection('projects/sat-io/open-datasets/GFSAD/GCEP30');

// Set time parameters - using a 10-year baseline
var startYear = 2010;
var endYear = 2022; // Extended to 2022 for more recent data
var startDate = ee.Date.fromYMD(startYear, 1, 1);
var endDate = ee.Date.fromYMD(endYear, 12, 31);

// Filter the MODIS collection to the study period
var filtered = modis
  .filter(ee.Filter.date(startDate, endDate))
  .filter(ee.Filter.bounds(pakistan));

// Cloud Masking function
var bitwiseExtract = function(input, fromBit, toBit) {
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
  var mask = ee.Number(1).leftShift(maskSize).subtract(1);
  return input.rightShift(fromBit).bitwiseAnd(mask);
};

var maskSnowAndClouds = function(image) {
  var summaryQa = image.select('SummaryQA');
  // Select pixels which are less than or equals to 1 (0 or 1)
  var qaMask = bitwiseExtract(summaryQa, 0, 1).lte(1);
  var maskedImage = image.updateMask(qaMask);
  return maskedImage.copyProperties(image, ['system:index', 'system:time_start']);
};

// Apply the cloud masking function to all images in the collection
var maskedCol = filtered.map(maskSnowAndClouds);

// Select NDVI band and scale values
var ndviCol = maskedCol.select('NDVI');
// MODIS NDVI values are scaled by 10000 and need to be divided by 10000
var scaleNdvi = function(image) {
  var scaled = image.divide(10000);
  return scaled.copyProperties(image, ['system:index', 'system:time_start']);
};

var ndviScaled = ndviCol.map(scaleNdvi);

// Create NDVI composite for every month
var years = ee.List.sequence(startYear, endYear);
var months = ee.List.sequence(1, 12);

// Map over the years and create a monthly average collection
var monthlyImages = years.map(function(year) {
  return months.map(function(month) {
    var filtered = ndviScaled
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .filter(ee.Filter.calendarRange(month, month, 'month'));
    var monthly = filtered.mean();
    return monthly.set({'month': month, 'year': year});
  });
}).flatten();

// Create a collection of monthly NDVI images
var monthlyCol = ee.ImageCollection.fromImages(monthlyImages);

// Compute Minimum NDVI for each month across all years
var monthlyMinImages = months.map(function(month) {
    var filtered = monthlyCol.filter(ee.Filter.eq('month', month));
    var monthlyMin = filtered.min();
    return monthlyMin.set('month', month);
});
var monthlyMin = ee.ImageCollection.fromImages(monthlyMinImages);

// Compute Maximum NDVI for each month across all years
var monthlyMaxImages = months.map(function(month) {
    var filtered = monthlyCol.filter(ee.Filter.eq('month', month));
    var monthlyMax = filtered.max();
    return monthlyMax.set('month', month);
});
var monthlyMax = ee.ImageCollection.fromImages(monthlyMaxImages);

// Create UI elements for interactive analysis
// Year and month selectors
var yearSelect = ui.Select({
  items: [2018, 2019, 2020, 2021, 2022].map(String),
  value: '2022',
  onChange: updateVCI,
  style: {stretch: 'horizontal'}
});

var monthSelect = ui.Select({
  items: [
    {label: 'January', value: '1'},
    {label: 'February', value: '2'},
    {label: 'March', value: '3'},
    {label: 'April', value: '4'},
    {label: 'May', value: '5'},
    {label: 'June', value: '6'},
    {label: 'July', value: '7'},
    {label: 'August', value: '8'},
    {label: 'September', value: '9'},
    {label: 'October', value: '10'},
    {label: 'November', value: '11'},
    {label: 'December', value: '12'}
  ],
  value: '5',
  onChange: updateVCI,
  style: {stretch: 'horizontal'}
});

// Function to calculate and display VCI for selected year and month
function updateVCI() {
  // Clear previous layers
  Map.layers().reset();
  
  // Add Pakistan boundary
  Map.addLayer(pakistan, {color: 'black'}, 'Pakistan Boundary');
  
  // Get currently selected year and month
  var currentYear = parseInt(yearSelect.getValue());
  var currentMonth = parseInt(monthSelect.getValue());
  
  // Filter monthly collection to the selected year and month
  var filtered = monthlyCol
    .filter(ee.Filter.eq('year', currentYear))
    .filter(ee.Filter.eq('month', currentMonth));
    
  // Get the NDVI image for the selected month
  var currentMonthNdvi = ee.Image(filtered.first());
  
  // Get the historical min and max NDVI for the selected month
  var minNdvi = ee.Image(monthlyMin.filter(ee.Filter.eq('month', currentMonth)).first());
  var maxNdvi = ee.Image(monthlyMax.filter(ee.Filter.eq('month', currentMonth)).first());
  
  // Calculate VCI = 100 * (NDVI - min) / (max - min)
  var image = ee.Image.cat([currentMonthNdvi, minNdvi, maxNdvi]).rename(['ndvi', 'min', 'max']);
  var vci = image.expression('100 * (ndvi - min) / (max - min)',
      {'ndvi': image.select('ndvi'),
       'min': image.select('min'),
       'max': image.select('max')
      }).rename('vci');
  
  // Get cropland mask
  var cropLand = gcep30.mosaic().eq(2);
  
  // Apply cropland mask to VCI
  var vciMasked = vci.updateMask(cropLand);
  
  // Classify VCI into drought categories
  var droughtCondition = vciMasked
    .where(vciMasked.lt(20), 1)    // Extreme Drought
    .where(vciMasked.gte(20).and(vciMasked.lt(40)), 2)  // Severe Drought
    .where(vciMasked.gte(40).and(vciMasked.lt(60)), 3)  // Moderate Drought
    .where(vciMasked.gte(60).and(vciMasked.lt(80)), 4)  // Normal
    .where(vciMasked.gte(80), 5);  // Above Normal
  
  // Define visualization parameters
  var vciPalette = ['#a50026','#d73027','#f46d43','#fdae61',
    '#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837'];
  var vciVisParams = {min: 0, max: 100, palette: vciPalette};
  
  var droughtPalette = ['#730000', '#e60000', '#ffaa00', '#66bd63', '#1a9850'];
  var droughtParams = {min: 1, max: 5, palette: droughtPalette};
  
  // Add layers to map
  Map.addLayer(vci.clip(pakistan), vciVisParams, 'VCI');
  Map.addLayer(vciMasked.clip(pakistan), vciVisParams, 'VCI (Cropland Only)', false);
  Map.addLayer(droughtCondition.clip(pakistan), droughtParams, 'Drought Condition');
  
  // Update the title with current date
  var monthNames = ['January', 'February', 'March', 'April', 'May', 'June',
                   'July', 'August', 'September', 'October', 'November', 'December'];
  titleLabel.setValue('Pakistan Drought Monitoring - ' + 
                      monthNames[currentMonth-1] + ' ' + currentYear);
  
  // Extract extreme drought areas (VCI < 20)
  var extremeDrought = droughtCondition.eq(1);
  Map.addLayer(extremeDrought.clip(pakistan), {palette: ['#FF0000']}, 'Extreme Drought Areas', false);
  
  // Calculate drought statistics
  var droughtStats = droughtCondition.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: pakistan,
    scale: 250,
    maxPixels: 1e9
  });
  
  // Print statistics to console
  var monthNames = ['January', 'February', 'March', 'April', 'May', 'June',
                   'July', 'August', 'September', 'October', 'November', 'December'];
  print('Drought Statistics for ' + monthNames[currentMonth-1] + ' ' + currentYear + ':', 
        droughtStats.get('vci'));
}

// Create control panel
var controlPanel = ui.Panel({
  style: {
    position: 'top-right',
    padding: '8px 15px',
    width: '300px'
  }
});

// Title for control panel
var titleLabel = ui.Label({
  value: 'Pakistan Drought Monitoring',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 10px 0'
  }
});

// Add UI elements to control panel
controlPanel.add(titleLabel);
controlPanel.add(ui.Label('Select Year:'));
controlPanel.add(yearSelect);
controlPanel.add(ui.Label('Select Month:'));
controlPanel.add(monthSelect);

// Create legend
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    backgroundColor: 'rgba(255, 255, 255, 0.8)'
  }
});

// Create legend title
var legendTitle = ui.Label({
  value: 'Drought Classification',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    margin: '0 0 4px 0',
    padding: '0'
  }
});

// Add title to legend
legend.add(legendTitle);

// Create a function to generate a legend item
function makeLegendItem(color, name) {
  var colorBox = ui.Label({
    style: {
      backgroundColor: color,
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });
  
  var description = ui.Label({
    value: name,
    style: {margin: '0 0 4px 6px'}
  });
  
  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
}

// Add legend items
legend.add(makeLegendItem('#730000', 'Extreme Drought (VCI < 20)'));
legend.add(makeLegendItem('#e60000', 'Severe Drought (VCI 20-40)'));
legend.add(makeLegendItem('#ffaa00', 'Moderate/Fair (VCI 40-60)'));
legend.add(makeLegendItem('#66bd63', 'Normal (VCI 60-80)'));
legend.add(makeLegendItem('#1a9850', 'Above Normal (VCI > 80)'));

// Add VCI explanation
var infoPanel = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px',
    backgroundColor: 'rgba(255, 255, 255, 0.8)',
    width: '300px'
  }
});

var infoTitle = ui.Label({
  value: 'About Vegetation Condition Index (VCI)',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    margin: '0 0 4px 0'
  }
});

var infoText = ui.Label({
  value: 'VCI compares current NDVI to the historical range (min and max) to evaluate vegetation health. ' +
         'Low VCI values indicate drought conditions, while high values indicate favorable conditions. ' +
         'This analysis focuses on cropland areas to monitor agricultural drought.',
  style: {fontSize: '12px'}
});

infoPanel.add(infoTitle);
infoPanel.add(infoText);

// Add a note about data sources
infoPanel.add(ui.Label({
  value: 'Data sources: MODIS NDVI (2010-2022) and Global Cropland Extent',
  style: {fontSize: '10px', fontStyle: 'italic', margin: '10px 0 0 0'}
}));

// Add panels to map
Map.add(controlPanel);
Map.add(legend);
Map.add(infoPanel);

// Center map on Pakistan
Map.centerObject(pakistan, 5);

// Initialize the map with the default selection
updateVCI();

// Add export functionality
var exportButton = ui.Button({
  label: 'Export Current View to Drive',
  onClick: function() {
    var currentYear = parseInt(yearSelect.getValue());
    var currentMonth = parseInt(monthSelect.getValue());
    var monthNames = ['January', 'February', 'March', 'April', 'May', 'June',
                     'July', 'August', 'September', 'October', 'November', 'December'];
    var monthName = monthNames[currentMonth-1];
    
    // Get the drought condition layer
    var filtered = monthlyCol
      .filter(ee.Filter.eq('year', currentYear))
      .filter(ee.Filter.eq('month', currentMonth));
    var currentMonthNdvi = ee.Image(filtered.first());
    var minNdvi = ee.Image(monthlyMin.filter(ee.Filter.eq('month', currentMonth)).first());
    var maxNdvi = ee.Image(monthlyMax.filter(ee.Filter.eq('month', currentMonth)).first());
    var image = ee.Image.cat([currentMonthNdvi, minNdvi, maxNdvi]).rename(['ndvi', 'min', 'max']);
    var vci = image.expression('100 * (ndvi - min) / (max - min)',
        {'ndvi': image.select('ndvi'),
         'min': image.select('min'),
         'max': image.select('max')
        }).rename('vci');
    
    var cropLand = gcep30.mosaic().eq(2);
    var vciMasked = vci.updateMask(cropLand);
    
    var droughtCondition = vciMasked
      .where(vciMasked.lt(20), 1)
      .where(vciMasked.gte(20).and(vciMasked.lt(40)), 2)
      .where(vciMasked.gte(40).and(vciMasked.lt(60)), 3)
      .where(vciMasked.gte(60).and(vciMasked.lt(80)), 4)
      .where(vciMasked.gte(80), 5);
    
    // Export the image to Drive
    Export.image.toDrive({
      image: ee.Image.cat([vci, droughtCondition.rename('drought_class')]).clip(pakistan),
      description: 'Pakistan_Drought_' + monthName + '_' + currentYear,
      scale: 250,
      region: pakistan.geometry(),
      maxPixels: 1e9
    });
    
    print('Exporting drought data to Google Drive as "Pakistan_Drought_' + 
          monthName + '_' + currentYear + '"');
  }
});

// Add export button to control panel
controlPanel.add(ui.Panel([exportButton], ui.Panel.Layout.Flow('horizontal')));

// Add a small map for administrative boundaries reference
var adminMap = ui.Map();
var adminPanel = ui.Panel({
  widgets: [
    ui.Label('Administrative Boundaries', {fontWeight: 'bold'}),
    adminMap
  ],
  style: {
    position: 'bottom-right',
    width: '300px',
    height: '200px'
  }
});

// Add provinces to the admin map
var provinces = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Pakistan'));
  
adminMap.addLayer(provinces, {color: 'blue'}, 'Provinces');
adminMap.addLayer(pakistan, {color: 'black'}, 'Country Boundary');
adminMap.setControlVisibility({
  all: false,
  layerList: true
});
adminMap.setCenter(70, 30, 5);

// Add link to get full analysis report
controlPanel.add(ui.Label({
  value: 'Generate Comprehensive Analysis Report',
  style: {
    color: 'blue',
    textDecoration: 'underline',
    margin: '10px 0 0 0'
  }
}));