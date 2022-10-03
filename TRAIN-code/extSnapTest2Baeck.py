# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:56:42 2021
@author: crisj
"""
import sys
#sys.path.append('C:\Users\vaak\.snap\auxdata\snaphu-v1.4.2_win64\bin')
#sys.path.append('C:\\Users\\vaak\\.snap\\snap-python\snappy')
sys.path.append('C:\\Users\\crisj\\.snap\\snap-python\\snappy')
sys.path.append('C:\\Users\\crisj\\.conda\\envs\\snap_env\\Lib') # anaconda environment created for this script was 'snap_env'
import os
from snappy import GPF
from snappy import ProductIO
from snappy import HashMap
from snappy import PixelPos, GeoPos, Band
#from snappy import jpy
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import numpy as np
import json
import datetime 
from tqdm import tqdm


import subprocess
import oct2py
from oct2py import octave

 # functions
# Hashmap is used to give us access to all JAVA operators
#HashMap = jpy.get_type('java.util.HashMap')

parameters = HashMap()

def read(filename):
    print('Reading...')
    return ProductIO.readProduct(filename)

def topsar_split(product,IW,firstBurstIndex,lastBurstIndex):
    print('Apply TOPSAR Split...')
    parameters = HashMap()
    parameters.put('subswath', IW)
    parameters.put('firstBurstIndex',firstBurstIndex ) # added by me
    parameters.put('lastBurstIndex',lastBurstIndex ) # added by me
    parameters.put('selectedPolarisations', 'VV')
    output=GPF.createProduct("TOPSAR-Split", parameters, product)
    return output

def apply_orbit_file(product):
    print('Applying orbit file ...')
    parameters.put("Orbit State Vectors", "Sentinel Precise (Auto Download)")
    parameters.put("Polynomial Degree", 3)
    parameters.put("Do not fail if new orbit file is not found", True)
    return GPF.createProduct("Apply-Orbit-File", parameters, product)

def back_geocoding(product):
    print('back_geocoding ...')
    parameters.put("Digital Elevation Model", "SRTM 1Sec HGT (Auto Download)")
    parameters.put("DEM Resampling Method", "BILINEAR_INTERPOLATION")
    parameters.put("Resampling Type", "BILINEAR_INTERPOLATION")
    parameters.put("Mask out areas with no elevation", True)
    parameters.put("Output Deramp and Demod Phase", True)
    parameters.put("Disable Reramp", False)
    return GPF.createProduct("Back-Geocoding", parameters, product)

def Enhanced_Spectral_Diversity(product):
    parameters = HashMap()
#    parameters.put("fineWinWidthStr2",512)
#    parameters.put("fineWinHeightStr",512)
#    parameters.put("fineWinAccAzimuth",16)
#    parameters.put("fineWinAccRange",16)
#    parameters.put("fineWinOversampling",128)
#    parameters.put("xCorrThreshold",0.1)
#    parameters.put("cohThreshold",0.3)
#    parameters.put("numBlocksPerOverlap",10)
#    parameters.put("esdEstimator",'Periodogram')
#    parameters.put("weightFunc",'Inv Quadratic')
#    parameters.put("temporalBaselineType",'Number of images')
#    parameters.put("maxTemporalBaseline",4)
#    parameters.put("integrationMethod",'L1 and L2')
#    parameters.put("doNotWriteTargetBands",False)
#    parameters.put("useSuppliedRangeShift",False)
#    parameters.put("overallRangeShift",0)
#    parameters.put("useSuppliedAzimuthShift",False)
#    parameters.put("overallAzimuthShift",0)
    return GPF.createProduct("Enhanced-Spectral-Diversity", parameters, product)

def interferogram(product):
    print('Creating interferogram ...')
    parameters.put("Subtract flat-earth phase", True)
    parameters.put("Degree of \"Flat Earth\" polynomial", 5)
    parameters.put("Number of \"Flat Earth\" estimation points", 501)
    parameters.put("Orbit interpolation degree", 3)
    parameters.put("Include coherence estimation", True)
    parameters.put("Square Pixel", True)
    # Added by mBergeron
    parameters.put("Output Elevation", True)
    #
    parameters.put("Independent Window Sizes", False)
    parameters.put("Coherence Azimuth Window Size", 10)
    parameters.put("Coherence Range Window Size", 2)
    return GPF.createProduct("Interferogram", parameters, product)

def topsar_deburst(source):  
    parameters = HashMap()
    parameters.put("Polarisations", "VV,VH")
    output=GPF.createProduct("TOPSAR-Deburst", parameters, source)
    return output

def topophase_removal(product):
    parameters.put("Orbit Interpolation Degree", 3)
    parameters.put("Digital Elevation Model", "SRTM 1Sec HGT (Auto Download)")
    parameters.put("Tile Extension[%]", 100)
    parameters.put("Output topographic phase band", True)
    parameters.put("Output elevation band", True)
    return GPF.createProduct("TopoPhaseRemoval", parameters, product)

#def Multilook(product, ML_nRgLooks,ML_nAzLooks):  
def Multilook(product, ML_nRgLooks):
    parameters = HashMap()
    parameters.put('grSquarePixel',True)
    parameters.put("nRgLooks", ML_nRgLooks)
    output=GPF.createProduct("Multilook", parameters, product)
    return output

def goldstein_phasefiltering(product):
    parameters.put("Adaptive Filter Exponent in(0,1]:", 1.0)
    parameters.put("FFT Size", 64)
    parameters.put("Window Size", 3)
    parameters.put("Use coherence mask", False)
    parameters.put("Coherence Threshold in[0,1]:", 0.2)
    return GPF.createProduct("GoldsteinPhaseFiltering", parameters, product)

def SNAPHU_export(product,SNAPHU_exp_folder):
    parameters = HashMap()
    parameters.put('targetFolder', SNAPHU_exp_folder) # 
    output = GPF.createProduct('SnaphuExport', parameters, product)
    ProductIO.writeProduct(output, SNAPHU_exp_folder, 'Snaphu')
    return (output)

def snaphu_unwrapping(product,target_Product_File,outFolder,filename):
    parameters = HashMap()
    parameters.put('targetProductFile', target_Product_File) # from SNAPHU_export
    parameters.put('outputFolder', outFolder)
    parameters.put('copyOutputAndDelete', 'Snaphu-unwrapping-after.vm')
    parameters.put('copyFilesTemplate', 'Snaphu-unwrapping-before.vm')
    product = GPF.createProduct('snaphu-unwrapping', parameters, product)
    ProductIO.writeProduct(product, filename+ '.dim', 'BEAM-DIMAP')
    print('Phase unwrapping performed successfully …')

def do_subset_band(source, wkt):
    print('\tSubsetting...')
    parameters = HashMap()
    parameters.put('geoRegion', wkt)
    #parameters.put('outputImageScaleInDb', True)
    output = GPF.createProduct('Subset', parameters, source)
    return output

def do_terrain_correction(source,band):
#def do_terrain_correction(source, proj, downsample):
    print('\tTerrain correction...')
    parameters = HashMap()
    parameters.put('demName', 'SRTM 1Sec HGT') # 'SRTM 3Sec'
#    parameters.put('imgResamplingMethod', 'BILINEAR_INTERPOLATION')
#    #parameters.put('mapProjection', proj)       # comment this line if no need to convert to UTM/WGS84, default is WGS84
    #parameters.put('saveProjectedLocalIncidenceAngle', False)
    parameters.put('sourceBands', band)
    parameters.put('saveSelectedSourceBand', True)
#    parameters.put('nodataValueAtSea', False)
    #parameters.put('pixelSpacingInMeter', 35)
#    while downsample == 1:                      # downsample: 1 -- need downsample to 40m, 0 -- no need to downsample
#        parameters.put('pixelSpacingInMeter', 40.0)
#        break
    output = GPF.createProduct('Terrain-Correction', parameters, source)
    return output

def addElevation(source):
    print('\tAdd elevation band...')
    parameters = HashMap()
    parameters.put('demName', 'SRTM 1Sec HGT') # 'SRTM 3Sec'
    output = GPF.createProduct('AddElevation', parameters, source)
    return output
    
def addLatLon(source):
    print('\tSubsetting...')
    parameters = HashMap()
    #parameters.put('outputImageScaleInDb', True)
    output = GPF.createProduct('Subset', parameters, source)
    return output

def write(product, filename):
    ProductIO.writeProduct(product, filename, "GeoTiff")
    
def write_BEAM_DIMAP_format(product, filename):
    print('Saving BEAM-DIMAP format...')
    ProductIO.writeProduct(product, filename + '.dim', 'BEAM-DIMAP')
#%%
#Input variables
#filename_1 = os.path.join(r'D:\PhD Info\InSAR\Examples\Ecuador_Galapagos','S1A_IW_SLC__1SDV_20170319T002614_20170319T002644_015753_019EFB_FA12.zip')
#filename_2 = os.path.join(r'D:\PhD Info\InSAR\Examples\Ecuador_Galapagos','S1A_IW_SLC__1SDV_20170331T002615_20170331T002645_015928_01A42E_7662.zip')
#out_filename =os.path.join(r'D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos','InSAR_pipeline_I')
#
#IW='IW3'
#firstBurstIndex = 1
#lastBurstIndex = 4

def numpy_snap(array,band_name):
    newband = Band(band_name,type(array[0][0].item()), np.size(array)[1], np.size(array)[1])
    newband.writePixels(0,0,np.size(array)[1],np.size(array)[0],array)
    return newband


def InSAR_pipeline_I(filename_1, filename_2,IW,firstBurstIndex,lastBurstIndex,out_filename):
    product_1 = read(filename_1)
    product_2 = read(filename_2)
    product_TOPSAR_1 = topsar_split(product_1,IW,firstBurstIndex,lastBurstIndex)
    product_TOPSAR_2 = topsar_split(product_2,IW,firstBurstIndex,lastBurstIndex)   
    product_orbitFile_1 = apply_orbit_file(product_TOPSAR_1)
    product_orbitFile_2 = apply_orbit_file(product_TOPSAR_2)
    product = back_geocoding([product_orbitFile_1, product_orbitFile_2])
    product = Enhanced_Spectral_Diversity(product)
    write_BEAM_DIMAP_format(product, out_filename)
    print("InSAR_pipeline_I complete")

def InSAR_pipeline_II(in_filename,ML_nRgLooks,out_filename_II):
    product = read(in_filename) # reads .dim
    product = interferogram(product)
    product = topsar_deburst(product)
    product = topophase_removal(product)
    # Insert subset here
    product = Multilook(product, ML_nRgLooks=ML_nRgLooks)
    product = goldstein_phasefiltering(product)    
    #product = do_terrain_correction(product)
    write_BEAM_DIMAP_format(product, out_filename_II)
    print("InSAR_pipeline_II complete")

"""
def InSAR_pipeline_III(in_filename_III,out_filename_III):
    
    
    product = read(in_filename_III) # reads .dim   
    band_names = product.getBandNames()
    #print("Band names: {}".format(", ".join(band_names)))
    a= format(", ".join(band_names)) # band names as a comma separated string 
    b= a.split(',') # split string into list
    #change band names as to mintpy accepts them
    product.getBand(b[3].strip()).setName('Phase_ifg')
    product.getBand(b[4].strip()).setName('coh')
    interferogram_TC = do_terrain_correction(product,band='Phase_ifg') # interferogram terrain correction
    
    interferogram_TC = addElevation(interferogram_TC)
    
    write_BEAM_DIMAP_format(interferogram_TC, out_filename_III+'_filt_int_sub_tc')
    coherence_TC = do_terrain_correction(product,band='coh') # coherence terrain correction
    write_BEAM_DIMAP_format(coherence_TC, out_filename_III+'_coh_tc')
    print("InSAR_pipeline_III complete")
 """  
    
def InSAR_pipeline_III(in_filename_III,out_filename_III):
    
    product = read(in_filename_III) # reads .dim   
    band_names = product.getBandNames()
    #print(band_names)
    print("Band names: {}".format(", ".join(band_names)))
    a= format(", ".join(band_names)) # band names as a comma separated string 
    b= a.split(',') # split string into list
    #change band names as to mintpy accepts them
    product.getBand(b[3].strip()).setName('Phase_ifg')
    product.getBand(b[4].strip()).setName('coh')
    Interferogram = do_terrain_correction(product,band='Phase_ifg') # coherence terrain correction
    Coherence = do_terrain_correction(product,band='coh') # interferogram terrain correction
    Interferogram = addElevation(Interferogram)
    aband = Interferogram.addBand(Coherence.getBand(Coherence.getBandNames()[0]))
    write_BEAM_DIMAP_format(Interferogram, out_filename_III)
    print("InSAR_pipeline_III complete")

def atmos(a,b,c):
    pass
"""
# extract input parameters
with open(sys.argv[1], 'r') as file:
    data = file.read()
    
mydict_str = data
mydict = json.loads(mydict_str)  # decode json string into dictionnary
pipeline = mydict['pipeline']
# run pipeline    
"""
    
if __name__ == "__main__":
    c = np.float64(299792458)
    #raw1 =  read(filename_2)
    #raw2 =  read(filename_2)
    product2 = ProductIO.readProduct(r"./InSAR_pipeline_II.dim")
    product4 = ProductIO.readProduct(r"C:\Users\max\Downloads\SnappyPython\Pipeline4\Pipeline4\Pipeline_IIII.dim")
    product2 = addElevation(product2)
    data_phase = product4.getBand("Unw_Phase_ifg_20Apr2016_08Apr2016") 
    data_elevation = product2.getBand("elevation") 
    
    width = product2.getSceneRasterWidth()
    height = product2.getSceneRasterHeight()
    print(width)
    print(", ")
    print(height)
    phase = np.zeros((height,width),dtype=np.float64)
    phase_atmos = np.zeros((height,width),dtype=np.float64)
    incident_angle = np.zeros((height,width))
    lat = np.zeros((height,width))
    lon = np.zeros((height,width))
    slantRange = np.zeros((height,width))
    DEM = np.zeros((height,width),dtype=np.float64)
    
    ##### Slave metadata #####
    metadata = product2.getMetadataRoot()
    submeta = metadata.getElement("Slave_Metadata")
    slaveNames =  submeta.getElementNames()
    print(slaveNames)
    slaveMeta=submeta.getElementAt(0)
    slaveData = slaveMeta.getElementAt(0)
    string = slaveData.getAttributeNames()
    orbit1 = slaveData.getElementAt(0)
    time_start = orbit1.getAttributeString("time")
    print(time_start.format())
    time_start = time_start.format()
    
    #### Master metadata #####
    time_end =	product4.getStartTime()
    time_end = time_end.format()
    
    timstamp_format = r"%d-%b-%Y %H:%M:%S.%f"
    time_start = datetime.datetime.strptime(time_start, timstamp_format)
    time_end = datetime.datetime.strptime(time_end, timstamp_format)
    
    #### Get wavelength from metadata ####
    AbsMeta = metadata.getElement("Abstracted_Metadata")
    freq = AbsMeta.getAttributeString("radar_frequency") # In MHz
    freq = np.float64(freq)
    #print(np.float64(freq)+1)
    
    lamb = c/(freq*(10.**6.)) 
    print(lamb)
    
    #latitude = numpy.zeros((width,height))
    #longitude = numpy.zeros((width,height))
    latlon = product2.getSceneGeoCoding()
    var = GeoPos()
    pix = PixelPos(0,0)
    geoPos = latlon.getGeoPos(pix,var)
    print(geoPos.lat)
    print(geoPos.lon)
    for i in tqdm(range(0,height)):
        for j in range(0,width):
            pix = PixelPos(i,j)
            geoPos = latlon.getGeoPos(pix,var)
            lat[i][j] = geoPos.lat
            lon[i][j] = geoPos.lon
    print("FInished loop")
    slantRange = product2.getTiePointGrid("slant_range_time")
    incident_angle = product2.getTiePointGrid("incident_angle")
    #DEM = product2.getBand()
    data_phase.readPixels(0, 0, width, height, phase)
    data_elevation.readPixels(0, 0, width, height, DEM)

    #Latitude: lat
    # Longitude: lon
    # wrapped phase: phase
    # elevation:  DEM

    lat2 =np.reshape(lat,width* height, 1)
    lon2 =np.reshape(lon,width* height, 1)
    lonlat = np.vstack([lon2,lat2]).T
#   phase = np.reshape(phase,width* height, 1)
#   DEM = np.reshape(DEM,width* height, 1)
    phase = np.ravel(phase,order='C');
    DEM = np.ravel(DEM,order='C');
    lat2 = np.ravel(lat,order='C');
    DEM = np.ravel(lon,order='C');
    print(lonlat.shape)
    print(phase.shape)
    print(DEM.shape)

    octave.addpath('C:\\Users\\max\\Downloads\\SnappyPython\\TRAIN-master')
    octave.addpath('C:\\Users\\max\\Downloads\\SnappyPython\\TRAIN-master\matlab\\')
    oc = oct2py.Oct2Py()
    cor = oc.feval('C:\\Users\\max\\Downloads\\SnappyPython\\TRAIN-master\\BekaertLinear',DEM, lon2,lat2, phase,height,width,'C:\\Users\\max\\Downloads\\SnappyPython\\TRAIN-master\\','C:\\Users\\max\\Downloads\\SnappyPython\\TRAIN-master\\matlab')
    corected_phase = phase - cor
    corected_phase2= np.reshape(corected_phase,height, width)

"""
    
   
elif pipeline == "snaphu":
    pass
    
else:
    out_file4 = "FinalData"
    towrite = ProductIO.readProduct(mydict["out_filename_III"] + '_filt_int_sub_tc' + '.dim')
    ProductIO.writeProduct(towrite, out_file4 + '//first_treatment.nc', 'NetCDF4-CF')
    ProductIO.writeProduct(towrite, out_file4 + '//interf_phase.tif', 'GeoTIFF-BigTIFF')

"""


#%% Test one operator indiv.
#product = read(os.path.join(r'D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos','InSAR_pipeline_II.dim')) # reads .dim   
#band_names = product.getBandNames()
##print("Band names: {}".format(", ".join(band_names)))
#a= format(", ".join(band_names)) # band names as a comma separated string 
#b= a.split(',') # split string into list
##%%
#product.getBand(b[3].strip()).setName('Phase_ifg')
#product.getBand(b[4].strip()).setName('coh')
#interferogram_TC = do_terrain_correction(product,band='Phase_ifg') # interferogram
#write_BEAM_DIMAP_format(interferogram_TC, out_filename_III+'_filt_int_sub_tc')
#coherence_TC = do_terrain_correction(product,band='coh') # coherence
#write_BEAM_DIMAP_format(coherence_TC, out_filename_III+'_coh_tc')
#%%
# ML
#parameters = HashMap()
#parameters.put('grSquarePixel',True)
#parameters.put("nRgLooks", 6)
##parameters.put("nAzLooks", 4)
##parameters.put('outputIntensity',False)
#output=GPF.createProduct("Multilook", parameters, product)
#write_BEAM_DIMAP_format(output, os.path.join(r'D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos','Multilook'))
# Terrain correction
#product = do_terrain_correction(product,band=b[3].strip())
#write_BEAM_DIMAP_format(product, os.path.join(r'D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos','TC_ifg'))
#%% Test unwrapping

################################
#product = read(os.path.join(r'D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos','subset_0_of_InSAR_pipeline_II.dim')) # reads .dim       
################################

#SNAPHU_exp_folder = 'D:\\PhD Info\\InSAR\\Examples\\SNAPPY_Ecuador_Galapagos\\exp2\\a' #'D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos\exp_snappy'
#SNAPHU_export(product,SNAPHU_exp_folder)

##################################
#SNAPHU_exp_folder='D:\PhD Info\InSAR\Examples\SNAPPY_Ecuador_Galapagos\exp_subset\subset_0_of_InSAR_pipeline_II'
##################################

"""
target_Product_File=SNAPHU_exp_folder # export
outFolder=SNAPHU_exp_folder#+'\\'#'D:\\PhD Info\\InSAR\\Examples\\SNAPPY_Ecuador_Galapagos\\exp_snappy\\InSAR_pipeline_II'
filename='date_name'+'_unw'
parameters = HashMap()
parameters.put('targetProductFile', target_Product_File) # from SNAPHU_export
parameters.put('outputFolder', outFolder)
#parameters.put('copyOutputAndDelete', 'Snaphu-unwrapping-after.vm')
#parameters.put('copyFilesTemplate', 'Snaphu-unwrapping-before.vm')
product = GPF.createProduct('snaphu-unwrapping', parameters, product)
#ProductIO.writeProduct(product, filename, 'BEAM-DIMAP')
ProductIO.writeProduct(product, filename, 'ENVI')
"""