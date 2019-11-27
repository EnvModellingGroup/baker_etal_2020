import qgis.core as qgis
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
import gdal
import ogr
from PyQt5.QtCore import *
import os
import argparse
import sys
import csv
from subprocess import call
import shutil
import tempfile
from qgis.PyQt.QtCore import QVariant

qgis_prefix="/usr"    
qgis.QgsApplication.setPrefixPath(qgis_prefix, True) 
qgs = qgis.QgsApplication([], False)
qgs.initQgis()
# a few hardcoded things...
rscript = "maxent_severn.R"
# make our temp dir to store some files from R
if not os.path.exists('temp'):
    os.makedirs('temp')
scenarios = ["pd","barrage"]
# unlike the no_web version of this script, we can't run this in parallel as results depend on others

def process_taxon(t,params,output_dir):

    # make temp dir
    temp_dir = tempfile.mkdtemp(prefix=t+"_", dir='temp')

    # call R script
    result_filename = '%s_out.txt' % t
    with open(result_filename, 'w') as result:
        call(['Rscript', rscript, temp_dir, params[1], params[2]],stdout=result,stderr=result);
    
    # parse R output for thresholds
    try:
        with open(os.path.join(temp_dir,"temp.log")) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        found_threshold = False
        found_ssb = False
        threshold = 0
        AUC = 0
        ssb = []
        k_fold_AUCs = []
        for c in content:
            if c.startswith("K-FOLD AUC"):
                tokens = c.split(':')
                k_fold_AUCs = tokens[1].split(',')[0:-1] # ignore last comma
                k_fold_AUCs = [float(i) for i in [float(i) for i in k_fold_AUCs]]
                k_fold_AUCs.append(sum(k_fold_AUCs)/5.0) # add average
            if (c.startswith("AUC")):
                tokens = c.split(':')
                AUC = float(tokens[1])
            if (c.startswith("Threshold")):
                tokens = c.split(":")
                threshold = float(tokens[1])
                found_threshold = False
            if (c.startswith("SSB evaluation")):
                tokens = c.split(':')
                ssb = tokens[1].split(",") # three of these

    except IOError:
        print("\tFailed to find temp.log for "+t+". Carrying on.")
        return [t,False]

    # mv output back
    try:
        shutil.move(os.path.join(temp_dir, "temp_pd.tif"),os.path.join(output_dir,t+"_pd.tif"))
        shutil.move(os.path.join(temp_dir, "temp_barrage.tif"),os.path.join(output_dir,t+"_barrage.tif"))
        shutil.move(os.path.join(temp_dir, "temp_response.pdf"),os.path.join(output_dir,t+"_response.pdf"))
        shutil.move(os.path.join(temp_dir,"temp.log"),os.path.join(output_dir,t+"_R.log"))
    except IOError:
        try:
            shutil.move(os.path.join(temp_dir,"temp.log"),os.path.join(output_dir,t+"_R.log"))
            print("\tError copying raster files for "+t+". Carrying on.")
            return [t,False]
        except IOError:
            print("\tTotally failed to run Maxent for "+t+". Carrying on.")
            return [t,False]


    # need to loop over all rasters
    area_for_t = []
    area_for_t.append(t)
    for s in scenarios:
        shutil.copy(os.path.join(output_dir,t+"_"+s+".tif"),os.path.join(temp_dir, "temp.tif"))
        
        # calulate area
        area = area_from_raster(threshold, temp_dir)

        # mv shapefile back
        try:
            shutil.move(os.path.join(temp_dir, "temp.shp"),os.path.join(output_dir,t+"_"+s+".shp"))
            shutil.move(os.path.join(temp_dir, "temp.shx"),os.path.join(output_dir,t+"_"+s+".shx"))
            shutil.move(os.path.join(temp_dir, "temp.prj"),os.path.join(output_dir,t+"_"+s+".prj"))
            shutil.move(os.path.join(temp_dir, "temp.dbf"),os.path.join(output_dir,t+"_"+s+".dbf"))
        except IOError:
            print("\tError copying shapefiles for "+t+". Carrying on.")
            return [t,False]

        area_for_t.append(area)
        os.remove(os.path.join(temp_dir,"temp.tif"))
        os.remove(os.path.join(temp_dir,"temp2.tif"))


    try:
        area_for_t.append(AUC)
        area_for_t.extend(ssb)
        area_for_t.append(threshold)
        area_for_t.extend(k_fold_AUCs)
    except IndexError:
        return [t,False]

    # save data
    # this  was in the main taxa loop, but in parallel won't work, so, writing out lot's of mini_files
    # write the output to csv - taxa,area,area_RCP2,area_RCP...
    output_file = os.path.join(output_dir,t+"_sdm.csv")
    with open(output_file, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(area_for_t) 

    return [t,True]



def main():

    from random import shuffle

    # do stuff
    parser = argparse.ArgumentParser(
         prog="Run SDM models for Severn Barrage",
         description="Run MaxEnt model on all species listed for without and with barrage, but using previous output as environment layer",
         )
    parser.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            'input_file', 
            metavar='input_file',
            help="Your input csv file containing species to run and which layers to use. species_layers_web.csv"
            )
    parser.add_argument(
            'output_file',
            metavar='output_file',
            help="The output file. Another csv containing taxa and areas calculated."
            )
    parser.add_argument(
            'output_dir',
            metavar='output_dir',
            help="The output dir. "
            )
    parser.add_argument(
            'previous_output',
            metavar='previous_output',
            help="The previous output dir for the standard SDMs"
            )

    args = parser.parse_args()
    verbose = args.verbose
    input_file = args.input_file
    output_file = args.output_file
    output_dir = args.output_dir
    previous_output = args.previous_output


    runs = {}
    with open(input_file, 'r') as f:
        csvreader = csv.reader(f, delimiter=",")
        i = 0
        for row in csvreader:
            if (i == 0):
                headers = row
                i += 1
            else:
                runs[row[0]] = row[1:]
                i += 1

    # taxa have to be run in a certain order to get the food web set up
    # this is the order and we don't have to run some, so we can skip these and use results from the non-web version
    sdm_taxa = ['whelk', 'brown_shrimp', 'velvet_swimming_crab', 'catworm', 'sand_goby', 'bib', 'poor_cod', 'dover_sole']
    taxa = runs.keys()
    # loop over taxa
    for t in taxa:
        # copy the previous from the non-web verion results to the data dir. File name should be taxa_pd.tif and taxa_barrage.tif
        if not os.path.isfile(os.path.join("../data/unaltered",t+".tif")):
            os.remove(filename)
        if not os.path.isfile(os.path.join("../data/with_barrage",t+".tif")):
            os.remove(filename)
        shutil.copy(os.path.join(previous_output, t+"_pd.tif"),os.path.join("../data/unaltered",t+".tif"))
        shutil.copy(os.path.join(previous_output, t+"_barrage.tif"),os.path.join("../data/with_barrage",t+".tif"))
    # we then overwrite some of these as we go, but these grab the taxa that don't have any links (i.e. bottom of the food web)

    success = []
    for t in sdm_taxa:
        success.append(process_taxon(t,runs[t], output_dir))
        # copy the output to the data dir fir use with the next one
        shutil.copy(os.path.join(output_dir, t+"_pd.tif"),os.path.join("../data/unaltered",t+".tif"))
        shutil.copy(os.path.join(output_dir, t+"_barrage.tif"),os.path.join("../data/with_barrage",t+".tif"))
        # these layers are then available for use in the next taxon, e.g. dover_sole will use the layers from poor cod, bib and velvet swimming crab
        # which in turn use layers from the previous run or other from this run.

 
    # join together all the little CSV files
    # set up our output csv
    # write the output to csv - taxa,area,area_RCP2,area_RCP...
    if not os.path.isfile(output_file):
        with open(output_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Taxon', 'Area', "Area Barrage", "AUC", "SSB 1", "SSB 2", "SSB AUC", "Threshold", "k-fold AUC 1", "k-fold AUC 2", "k-fold AUC 3", "k-fold AUC 4", "k-fold AUC 5", "k-fold Ave. AUC" ])
            
    with open(output_file, 'a') as f:
        for s in success:
            t = s[0]
            ok = s[1]
            if ok:
                t_output_file = os.path.join(output_dir,t+"_sdm.csv")
                with open(t_output_file,"r") as infile:
                    f.write(infile.read())

    return

def area_from_raster(threshold, temp_dir):

        # load raster produced by r
        MaxEnt_prediction = qgis.QgsRasterLayer(os.path.join(temp_dir, "temp.tif"))

        # work out the threshold part
        entries = []
        # Define band1
        me_pred = QgsRasterCalculatorEntry()
        me_pred.ref = 'temp@1'
        me_pred.raster = MaxEnt_prediction
        me_pred.bandNumber = 1
        entries.append( me_pred )
        # Process calculation with input extent and resolution
        calc = QgsRasterCalculator( 'temp@1 > '+str(threshold), os.path.join(temp_dir, "temp2.tif"), 'GTiff', MaxEnt_prediction.extent(), MaxEnt_prediction.width(), MaxEnt_prediction.height(), entries )
        calc.processCalculation()

        # here we switch to using gdal so we can polygonize the raster
        sourceRaster = gdal.Open(os.path.join(temp_dir, "temp2.tif"))
        band = sourceRaster.GetRasterBand(1)
        raster_field = ogr.FieldDefn('presence', ogr.OFTInteger)
        outShapefile = "polygonized"
        dest_srs = ogr.osr.SpatialReference()
        dest_srs.ImportFromEPSG(32630)
        driver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(os.path.join(temp_dir,"temp.shp")):
            driver.DeleteDataSource(os.path.join(temp_dir,"temp.shp"))
        outDatasource = driver.CreateDataSource(os.path.join(temp_dir, "temp.shp"))
        outLayer = outDatasource.CreateLayer("polygonized", srs=dest_srs)
        outLayer.CreateField(raster_field)
        gdal.Polygonize( band, None, outLayer, 0, [], callback=None )
        outDatasource.Destroy()
        sourceRaster = None

        # we now have a polygon version of the raster.
        # load this in using Q and then delete the zero-parts (1 is where the critter lives)
        polygon = qgis.QgsVectorLayer(os.path.join(temp_dir, "temp.shp"), "Shapefile", "ogr")
        # check this is valid
        if not polygon.isValid():
            print("Error loading file ", os.path.join(temp_dir, "temp.shp"))
        # build a request to filter the features based on an attribute
        request = qgis.QgsFeatureRequest().setFilterExpression('"presence" != 1')
        polygon.startEditing()
        # loop over the features and delete
        for f in polygon.getFeatures(request):
            polygon.deleteFeature(f.id())

        # now calculate the areas
        lProvider = polygon.dataProvider()
        lProvider.addAttributes( [ qgis.QgsField("Area",QVariant.Double) ] )
        #lProvider.setFields()
        

        sum_area = 0
        # set up the calculator
        elps_crs = qgis.QgsCoordinateReferenceSystem()
        elps_crs.createFromUserInput("WGS84")
        lyr_crs = polygon.crs()
        trans_context = qgis.QgsCoordinateTransformContext()
        trans_context.calculateDatumTransforms(lyr_crs, elps_crs)
        calculator = qgis.QgsDistanceArea()
        calculator.setEllipsoid('WGS84')
        calculator.setSourceCrs(lyr_crs, trans_context)

        for gFeat in polygon.getFeatures():
            geom = gFeat.geometry()
            area = calculator.measureArea(geom)
            gFeat
            sum_area += area
        polygon.commitChanges()
        
        sum_area = sum_area / 1000000. # convert to sq km
        return sum_area



if __name__ == "__main__":
    main()

