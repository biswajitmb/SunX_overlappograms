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
#=======================

DataDir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/'

Store_outputs = True

plate_scale = 2.8 #in arc-sec for which the EM cube would be

left_xy_FOV = np.array([-1030,-1030]) #in arc-sec
right_xy_FOV = np.array([+1030,+1030])
####

AIA_images_dir = os.path.join(DataDir,'AIA')
XRT_images_dir = os.path.join(DataDir,'XRT')

delta_FOV = right_xy_FOV-left_xy_FOV
new_dim = (delta_FOV/plate_scale).astype('int')

#Prepare AIA data
AIA_files = glob.glob(os.path.join(AIA_images_dir,'*image_lev1.fits'))
ffa=sorted(np.array(AIA_files))
amaps=sunpy.map.Map(ffa)
wvn0 = [m.meta['wavelnth'] for m in amaps]
srt_id = sorted(range(len(wvn0)), key=wvn0.__getitem__)
amaps = [amaps[i] for i in srt_id]
print([m.meta['wavelnth'] for m in amaps])
# aiaprep the images, may take a while to run
aprep=[]
for m in amaps:
    m_temp = update_pointing(m)
    m_observer_fixed = fix_observer_location(m_temp)
    aprep.append(register(m_temp))
#  Just save out the prepped submaps to quickly load in later

OutDir = os.path.join(AIA_images_dir,'rebined')
if os.path.isdir(OutDir) is False: os.mkdir(OutDir)
for m in aprep:
    bottom_left = SkyCoord(left_xy_FOV[0]*u.arcsec,left_xy_FOV[1]*u.arcsec, frame=m.coordinate_frame)
    top_right = SkyCoord(right_xy_FOV[0]*u.arcsec,right_xy_FOV[1]*u.arcsec, frame=m.coordinate_frame)
    mm = m.submap(bottom_left=bottom_left, top_right=top_right)
    mm = mm.resample(dimensions = new_dim*u.pix)
    wvn="{0:d}".format(1000+mm.meta['wavelnth'])
    wvn=wvn[1:]
    mm.peek()
    OutFile = 'AIA_rebin_level1.5_wvn'+wvn
    mm.save(os.path.join(OutDir,OutFile+'.fits'),overwrite='True')

#Prepare XRT data
XRT_files = sorted(np.array(glob.glob(os.path.join(XRT_images_dir,'prep','XRTcomp_Be_thin-Open.fits'))))
c_tx = 0
if len(XRT_files) > 0:
    OutDir = os.path.join(XRT_images_dir,'rebined')
    if os.path.isdir(OutDir) is False: os.mkdir(OutDir)
    for i in range(len(XRT_files)):
        xrt_file = XRT_files[i]#os.path.join(XRT_images_dir,'XRTcomp_'+i+'.fits')
        xrt_map = sunpy.map.Map(xrt_file)
        #xrt_map.peek()
        bottom_left = SkyCoord(left_xy_FOV[0]*u.arcsec,left_xy_FOV[1]*u.arcsec, frame=xrt_map.coordinate_frame)
        top_right = SkyCoord(right_xy_FOV[0]*u.arcsec,right_xy_FOV[1]*u.arcsec, frame=xrt_map.coordinate_frame)
        xrt_map = xrt_map.submap(bottom_left=bottom_left, top_right=top_right)
        xrt_map = xrt_map.resample(dimensions = new_dim*u.pix)
        OutFile = 'XRT_rebin_'+os.path.basename(xrt_file)
        xrt_map.save(os.path.join(OutDir,OutFile),overwrite='True')
        xrt_map.peek()

