import fwd as fd
import os, glob
import numpy as np
from scipy.io import readsav
from astropy.io import fits
#------------------------------------------#
#Inputs:

pointing_info_file = None
#pointing_info_file = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/FWD_pointing/FakeData/pointing_info_raw.dat'

DEM_file = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/DEMs/sparse_ObsDEM_AIA_XRT_dem.pkl'
plate_scale = 2.8

MaGIXS_resp = 'magixs_resp/D1Jul2024_magixs2_response_feldman_m_el_with_tables.fits'
resp_unit = 'electron'

Wvl_cal_file = '/Users/bmondal/BM_Works/MaGIXS/ForwardModel/magixs_resp/wavearray_30arcminfa_33disp_distortion.sav'

Vigneting_func_file = '/Users/bmondal/BM_Works/MaGIXS/ForwardModel/magixs_resp/nvigmap_magixs2.sav'

IncludePoisonNoise = True
NoFrames = 138 #Number of frames
IncludeDetecorNoise = True
detector_noise_file = 'magixs_resp/darks_for_biswajit.sav'
phot2elec_conv_file = 'magixs_resp/MaGIXS2_ph2elec.txt'

exposure = 1.3
psf = 16 #in arcsec
cmap = 'hot'

#Outputs
StoreOutputs = True
OutDir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/FWD_pointing/FakeData/MultiFrame'
namestr = 'MaGIXS2_MultiFrame_'
SavePlot = False

StoreRawDetector = True
RawDetectorDim = [1040,2152]
NonActivePix = [50, 8,4]

line_label = ['Mg XI','Fe XVII','Ne IX','Fe XVIII','Fe XVII','Fe XVII','Fe XVII','O VIII','O VII','O VII', 'N VII']
line_wvl = [9.17,12.12,13.45,14.20,15.01,15.26,17.05,18.97,21.60,22.1,24.78]

#------------------------------------------#



Dir = OutDir

dem = fd.util.load_obj(DEM_file[0:-4])
EM_map = dem['EMmap'] # in 1.0e26 cm-5
dem_logt = dem['logt']
ref_map = dem['ref_map']

detector_noise = readsav(detector_noise_file)['data'][4,:,:] #in the unit of electron or photon

m = fd.fwd(EM_map,dem_logt,ref_map)

FullSun_img_data = np.nan_to_num(ref_map.data)  # Replace NaNs with 0


m.magixs_fwd(FullSun_img_data,MaGIXS_resp,Wvl_cal_file,exposure = exposure,
        resp_unit = resp_unit, plate_scale = plate_scale, Vigneting_func_file = Vigneting_func_file,
        psf = psf,line_label = line_label,line_wvl = line_wvl,cmap=ref_map.cmap,
        StoreRawDetector = StoreRawDetector, RawDetectorDim=RawDetectorDim, 
        NonActivePix=NonActivePix, IncludePoisonNoise=IncludePoisonNoise, 
        IncludeDetecorNoise=IncludeDetecorNoise, detector_noise=detector_noise, 
        phot2elec_conv_file=phot2elec_conv_file, StoreOutputs=StoreOutputs, OutDir=OutDir,namestr = namestr,
        pointing_info_file = pointing_info_file,SavePlot=SavePlot, NoFrames = NoFrames)

'''
#Make average data of all frames
Overlappogram_file = os.path.join(Dir,namestr+'FWDmodel.fits')
m.make_average_detector(Overlappogram_file,os.path.join(Dir,'MaGIXS2_Average_FWDmodel.fits'))

##Run Inversion of the simulated overlappogram using 'https://github.com/ECCCO-mission/overlappogram'
#dataDir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/FWD_pointing/FakeData/MultiFrame/'
#overlappogram_file = 'MaGIXS2_Average_FWDmodel.fits'
#resp_file = 'magixs_resp/D2Jul2024_magixs2_response_feldman_m_el_with_tables.fits'
#gnt_file = 'magixs_resp/D16Sep2023_chi_beta_feldman_master_gnt_magixs1_with_tables.fits'
#configfile = 'config1_BM.toml'
#m.run_inversion(dataDir,overlappogram_file, resp_file, gnt_file,configfile, weight_file = None, OutDir=None, prefix = "MaGIXS2_inversion",FOV=1)



#Plot all frame data
EMin_file = os.path.join(Dir,namestr+'EMcubeInput.fits')
Overlappogram_file = os.path.join(Dir,namestr+'FWDmodel.fits')
m.plot_detector(ref_map,Overlappogram_file,EMin_file,cmap,namestr,line_label,line_wvl,FOV = 1, unit= 'photon', FrameRange = [0,5])


#Plot average data of all frames
Overlappogram_file  = os.path.join(Dir,'MaGIXS2_Average_FWDmodel.fits')
EMin_file = os.path.join(Dir,namestr+'EMcubeInput.fits')
SpecPure_file = None #os.path.join(Dir,'inversion','MaGIXS2_inversion_spectral_x2_5.0_5e-05_wpsf.fits')
m.plot_detector_avg(ref_map,Overlappogram_file,EMin_file,SpecPure_file=SpecPure_file,cmap=cmap,line_label=line_label,line_wvl=line_wvl,FOV = 1,Pscale = [2.8*2.8,2.8],unit='photon',specInd = [7,12,5])

'''

'''
## Not implemented yet
#Convert IDL emcube file to python readable format
IDLfile = ''
m = fd.make_em()
m.emcube_idl2pkl(IDLfile,DEM_file[0:-4],AIA_ref_file = None, XRT_ref_file = None)#, removeIDL = False)
'''


