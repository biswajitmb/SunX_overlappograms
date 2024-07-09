'''
Prepare AIA level-1 to 1.5. 
Load XRT data. Convert AIA and XRT to a given plate-scale and save the images
'''

# Import some of the stuff we will need
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy import units as u
from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table
from aiapy.calibrate import estimate_error
from aiapy.calibrate import register, update_pointing,fix_observer_location
import os,glob
import matplotlib.pyplot as plt
##########

class make_em(object):
    def __init__(self,date='20240716T00:00:00'):
        self.date = date

    #def download_aia(self,):

    #def download_xrt(self):

    #def prep_xrt(self,XRTdir = '',OutDir = ''):

    #def prepare_data(self,):

    #def run_sparse_em(self):

    def emcube_idl2pkl(self,IDLfile,OutFile,AIA_ref_file = None, XRT_ref_file = None, removeIDL = False, DEM_code = 'sparse_dem'):
        data = readsav(os.path.join(DataDir,idl_file))
        emcube = data['emcube'] #[logt, Y, X]
        emcube = emcube.T #[Y,X,logt]
        emcube_re = np.zeros([len(emcube[0,:,0]),len(emcube[:,1,0]),len(emcube[0,0,:])])
        for i in range(len(emcube[0,0,:])):
            emcube_re[:,:,i] = emcube[:,:,i].T #[X,Y,logT]
        
        if XRT_ref_file is not None: xrt_ref_map = sunpy.map.Map(XRT_ref_file)
        if AIA_ref_file is not None: aia_ref_map = sunpy.map.Map(AIA_ref_file)
        else: raise Exception("%% Error : Provide a valid AIA file")
        
        results = {'DEM_code': DEM_code}
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
        if XRT_ref_file is not None: results['ref_map_xrt'] = xrt_ref_map
        else: results['ref_map_xrt'] = None
        results['plate_scale'] = data['plate_scale']
        save_obj(results, os.path.join(OutDir,'sparse_'+idl_file[0:-4]))
        if removeIDL is True: os.remove(os.path.join(DataDir,idl_file))



        
