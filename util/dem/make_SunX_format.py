'''
Convert the MC_dem .idl output to sunx .pkl format 
'''

import numpy as np
import matplotlib.pyplot as plt
import glob, os, time
from scipy.io import readsav
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map
import sunx as ar

###
DataDir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/DEMs/'
idl_file = 'ObsDEM_AIA_XRT_dem.sav'
###

OutDir = DataDir

data = readsav(os.path.join(DataDir,idl_file))

IMG_DataDir = data['data_dir'].decode('utf-8') #directory where AIA and XRT input files exist
aia_files = data['aia_files']
xrt_files = data['xrt_files']

emcube = data['emcube'] #[logt, Y, X]
emcube = emcube.T #[Y,X,logt]
emcube_re = np.zeros([len(emcube[0,:,0]),len(emcube[:,1,0]),len(emcube[0,0,:])])

for i in range(len(emcube[0,0,:])):
    emcube_re[:,:,i] = emcube[:,:,i].T #[X,Y,logT]

if len(xrt_files) > 0:
    xrt_ref_map = sunpy.map.Map(os.path.join(IMG_DataDir,'XRT','rebined',xrt_files[0].decode('utf-8')))
aia_ref_map = sunpy.map.Map(os.path.join(IMG_DataDir,'AIA','rebined',aia_files[0].decode('utf-8')))


results = {'DEM_code': 'sparse_dem'}
results['instrument_passbands'] = data['chn_all']
results['EMmap'] = emcube_re
results['logt'] = data['use_lgtaxis']
results['ObsEM_unit'] = '1.0e26 cm-5'
#results['ObsDEMmap_Err'] = DEM_map_err
#results['ObsDEMmap_logt_midbins'] = mlogt
#results['ObsDEMmap_logtErr'] = DEM_map_logTerr
#results['chisq'] = chisq_map
#results['EMmap'] = EMmap
#results['ObsImg'] = Img
results['Tresp']  = data['tresp_all'] #[logt,chn]
results['Tresp_logt'] = data['use_lgtaxis']
results['ref_map'] = aia_ref_map
results['ref_map_xrt'] = xrt_ref_map
results['plate_scale'] = data['plate_scale']
ar.util.save_obj(results, os.path.join(OutDir,'sparse_'+idl_file[0:-4]))
os.remove(os.path.join(DataDir,idl_file))
