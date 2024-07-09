'''
Note: Implementing parallel processing requires 'spawn' method, with 'if __name__ == '__main__':'
to avoid running whole scripts in each child process. This requires to pass all the variables in parallel function
'''


# Import some of the stuff we will need
import numpy as np
import math
import matplotlib.pyplot as plt
import glob, os, time
from scipy.io import readsav
from scipy import interpolate
from scipy.ndimage import shift
# So can have single copy of demreg on system, and don't need copy in working directory
from sys import path as sys_path
# Change to your local copy's location...
sys_path.append('/Users/bmondal/BM_Works/softwares/demreg/python')
sys_path.append('/Users/bmondal/BM_Works/softwares/biswajit_python_lib')
from dn2dem_pos import dn2dem_pos

import astropy.time as atime
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map
import sunx as ar

from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table
from aiapy.calibrate import estimate_error

start_time = time.time()

#-------- inputs ---------
DataDir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240323T0330/'

AIA_tresp_file = '/Users/bmondal/BM_Works/MaGIXS/MaGIXS_1/XBP_1/data_old/temperature_response/aia_tresp_30072021.sav'
XRT_tresp_file = '/Users/bmondal/BM_Works/MaGIXS/MaGIXS_1/XBP_1/data_old/temperature_response/xrt_tresp_300720021.sav'

DEM_log = [5.6,7.0,0.1] #Grid [min, max, step] for DEM estimation

Store_outputs = False

OutDir = DataDir+'DEMs/'

AIA_passbands = ['094','131','171','193','211','335']
XRT_passbands = ['Be-thin']

plate_scale = 2.8 #in arc-sec for which the EM cube would be

#nCore = 10 #Nmber of core for which the program will parallel computation. If set 0 then it will not use parallel processing.
#-------------------------

All_passbands = AIA_passbands + XRT_passbands

AIA_images_dir = os.path.join(DataDir,'AIA')
XRT_images_dir = os.path.join(DataDir,'XRT')

#load XRT TStore_outputsresp:
xrt_tresp = readsav(XRT_tresp_file)
filters = np.array(xrt_tresp['filters'])
for i in range(len(filters)):
    filters[i]=filters[i].decode('utf-8')
units = xrt_tresp['units'].decode('utf-8')
xrt_tresp_logT = xrt_tresp['logt']

tresp_logt = np.arange(5.0,8.0,0.05)

nf = len(AIA_passbands)+len(XRT_passbands)
nt = len(tresp_logt)

trmatrix=np.zeros((nt,nf))

#load AIA Tresp:
aia_tresp = readsav(AIA_tresp_file)
for fil_aia in np.arange(len(aia_tresp['channels'])):
    aia_tresp['channels'][fil_aia] = aia_tresp['channels'][fil_aia].decode('utf-8')
chans=np.array(aia_tresp['channels'])
aia_tresp_logT = aia_tresp['logt']

left_xy_FOV = np.array([-1030,-1030]) #in arc-sec
right_xy_FOV = np.array([+1030,+1030])

delta_FOV = right_xy_FOV-left_xy_FOV

new_dim = (delta_FOV/plate_scale).astype('int')
Img = np.zeros((new_dim[0],new_dim[1],nf))

#Preparing data

c_t = 0
Inst = []
durs = [] #exposures 
for i in AIA_passbands:
    aia_file = os.path.join(AIA_images_dir,'AIA_level1.5_wvn'+i+'.fits')
    aia_map = sunpy.map.Map(aia_file)
    bottom_left = SkyCoord(left_xy_FOV[0]*u.arcsec,left_xy_FOV[1]*u.arcsec, frame=aia_map.coordinate_frame)
    top_right = SkyCoord(right_xy_FOV[0]*u.arcsec,right_xy_FOV[1]*u.arcsec, frame=aia_map.coordinate_frame)
    aia_map = aia_map.submap(bottom_left=bottom_left, top_right=top_right)
    aia_map = aia_map.resample(dimensions = new_dim*u.pix)
    exposure = aia_map.exposure_time.value
    #aia_map.peek()
    if c_t == 0: ref_map = aia_map
    Img[:,:,c_t] = aia_map.data[:,:]/exposure #DN/s/px
    ind = np.where(aia_tresp['channels'] == 'A'+format('%d'%int(i)))[0][0]
    tresp_intp_func =  interpolate.interp1d(aia_tresp_logT, aia_tresp['tr'][ind,:])
    tresp = tresp_intp_func(tresp_logt)
    trmatrix[:,c_t] = tresp * (plate_scale/0.6)**2 #corrected for 2.8" new pixcel size
    c_t += 1
    Inst += ['AIA']
    durs += [exposure]
c_tx = 0
if len(XRT_passbands) > 0:
    for i in XRT_passbands:
        xrt_file = os.path.join(XRT_images_dir,'XRTcomp_'+i+'.fits')
        xrt_map = sunpy.map.Map(xrt_file)
        #xrt_map.peek()
        bottom_left = SkyCoord(left_xy_FOV[0]*u.arcsec,left_xy_FOV[1]*u.arcsec, frame=xrt_map.coordinate_frame)
        top_right = SkyCoord(right_xy_FOV[0]*u.arcsec,right_xy_FOV[1]*u.arcsec, frame=xrt_map.coordinate_frame)
        xrt_map = xrt_map.submap(bottom_left=bottom_left, top_right=top_right)
        xrt_map = xrt_map.resample(dimensions = new_dim*u.pix)
        #xrt_map.peek()
        exposure = xrt_map.exposure_time.value
        Img[:,:,c_t] = aia_map.data/exposure #DN/s/px
        ind = np.where(filters == i)[0][0]
        tresp_intp_func =  interpolate.interp1d(xrt_tresp_logT, xrt_tresp['tr'][ind,:])
        tresp = tresp_intp_func(tresp_logt)
        trmatrix[:,c_t] = 2*tresp * (plate_scale/1)**2 ##Apply a cross-calibration factor of 2 and 
                                                       #corrected for 2.8" new pixcel size instead of 1" original XRT responce
        if c_tx == 0: ref_map_xrt = xrt_map
        c_t += 1 
        c_tx += 1
        Inst += ['XRT']
        durs += [exposure]
else: ref_map_xrt = None

temps = np.logspace(DEM_log[0],DEM_log[1],num=int((DEM_log[1]-DEM_log[0])/DEM_log[2]))  #DEM temperatures grids
# Temperature bin mid-points for DEM plotting
mlogt=([np.mean([(np.log10(temps[i])),np.log10((temps[i+1]))]) for i in np.arange(0,len(temps)-1)])   

#plt.close('all')
#fig = plt.figure(figsize=(8,8))
#ax = plt.subplot(projection=data0)
#data0.plot()

def EM2DEM(log10T,EM):
    dT = (shift(log10T, -1, cval=0.0) - shift(log10T, 1, cval=0.0)) * 0.5
    ntemps = len(log10T)
    dT[0] = log10T[1] - log10T[0]
    dT[ntemps-1] = (log10T[ntemps-1] - log10T[ntemps-2])
    DEM = EM / ((10**log10T)*np.log(10.) * dT)
    return DEM



DEM_map = np.zeros([new_dim[0],new_dim[1],len(mlogt)])
DEM_map_logTerr = np.zeros([new_dim[0],new_dim[1],len(mlogt)])
DEM_map_err = np.zeros([new_dim[0],new_dim[1],len(mlogt)])
chisq_map = np.zeros([new_dim[0],new_dim[1]])
for i in range(len(Img[:,0,0])):
    for j in range(len(Img[0,:,0])):
        print(i,j)
        data_dem=[]
        data_dem =  Img[i,j,:] #DN/px/s
        #cor_data=data_dem/degs
        #dn_in=cor_data#/durs  #DN/px/s
        #print('dn_in: ',dn_in)
        if all( ii > 0 for ii in data_dem) is True : #if all the chanels has positive values
            #Estimate errors from sunpy for AIA
            num_pix=1
            edn_in = []
            for k in range(nf):
                if Inst[k] == 'AIA':
                    aerr_temp=estimate_error((data_dem[k]*durs[k])*(u.ct/u.pix),int(AIA_passbands[k])*u.angstrom,num_pix) # Will download error table first time using
                    edn_in.append(aerr_temp.value[0])

                elif Inst[k] == 'XRT': 
                    edn_in.append(np.sqrt((data_dem[k] + (0.2*data_dem[k]))*durs[k])/durs[k])  #consider 20% systematic error for XRT
            edn_in = np.array(edn_in)
 
            # 1. Default - reg runs twice, 1st time to work out weight for constraint matrix, then regs with that
            #         Best option if don't know what doing, hence its the default 
            #         Probably best option of AIA data as well.
            dem0,edem0,elogt0,chisq0,dn_reg0=dn2dem_pos(data_dem,edn_in,trmatrix,tresp_logt,temps,gloci=1,emd_int=1) #gloci=0 is default behaviour
            DEM_map[i,j,:] = dem0
            DEM_map_logTerr[i,j,:] = elogt0
            DEM_map_err[i,j,:] = edem0
            chisq_map[i,j] = chisq0
            
            '''
            ##  Plot it all
            fig = plt.figure(figsize=(8, 4.5))
            plt.errorbar(mlogt,dem0,xerr=elogt0,yerr=edem0,fmt='or',ecolor='lightcoral', \
                         elinewidth=3, capsize=0,label='Def Self LWght, $\chi^2 =$ {:0.2f}'.format(chisq0))
            plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
            plt.ylabel('$\mathrm{DEM\;[cm^{-5}\;K^{-1}]}$')
            plt.ylim([1e16,4e22])

            plot_EM_locii = True
            if plot_EM_locii is True:
                for e in range(len(data_dem)):
                    plt.plot(tresp_logt,EM2DEM(tresp_logt,(data_dem[e]/trmatrix[:,e])))

            plt.xlim([5.6,7.0])
            plt.rcParams.update({'font.size': 16})
            plt.yscale('log')
            plt.legend()
            # plt.savefig('demregpy_aiapxl_slw.png',bbox_inches='tight')
            plt.show()
            '''

name = 'AIA_XRT_dem'

m=ar.fieldalign_model('config.dat')
DEM_logT = mlogt
Non_zero_pixcels = np.where(np.sum(DEM_map,axis=2) > 0)
EMmap = m.DemMap2EMmap(DEM_logT,DEM_map,Non_zero_pixcels)

if Store_outputs is True:
    results = {'Observed DEM Map': 'code: HK_dem'}
    results['instrument_passbands'] = All_passbands
    results['ObsDEMmap'] = DEM_map
    results['ObsDEMmap_Err'] = DEM_map_err
    results['ObsDEMmap_logt_midbins'] = mlogt
    results['ObsDEMmap_logtErr'] = DEM_map_logTerr
    results['chisq'] = chisq_map
    results['EMmap'] = EMmap
    results['ObsImg'] = Img
    results['Tresp']  = trmatrix
    results['Tresp_logt'] = tresp_logt
    results['ref_map'] = ref_map
    results['ref_map_xrt'] = ref_map_xrt
    ar.util.save_obj(results, OutDir+'ObsDEM_'+name)

print("END Calculation--- %s seconds ---" % (time.time() - start_time),'\n')
   
'''
#Estimate AIA degradation
nc = len(AIA_passbands)
degs=np.empty(nc)
for i in range(nc):
    degs[i]=degradation(wvn[i],time_obs,calibration_version=10)
'''



'''
for img in range(len(data_files)):
    data = load_obj(data_files[img][0:-4])
    ind = 0
    #ch_ind = [0,1,2,3,4,5,6] #Channe;s indices to be consider for DEM estimation
    ch_ind = [0,1,2,3,4,5] #for AIA
    channels = np.array(data['instruments'][Inst[ind]]['channels'])[ch_ind]
    data0 = data['instruments'][Inst[ind]]['ObsIMG'][channels[0]]
    name = '2021-07-30T18'
    if Inst[ind] == 'AIA':
        durs = np.array([data['instruments'][Inst[ind]]['ObsIMG'][m].meta['exptime'] for m in channels])
        #Estimate the degradation; If will take some time thus for the 1st time it will save the 
        # degradation in a .pkl file and for the nexttime it will restore from that file.
        time_obs = data0.date 
        name = time_obs.value.split(':')[0]
        fname = DataFileDir + name+'_AIA_degradations.pkl'
        wvn = [ int(m[1::]) for m in channels] * u.angstrom
        nc = len(channels)
        if os.path.isfile(fname) is False:
            degs=np.empty(nc)
            for i in range(nc):
                degs[i]=degradation(wvn[i],time_obs,calibration_version=10)
            caldata = {'AIA calibration data': 'method: Sunpy'}
            caldata['date'] = time_obs
            caldata['channels'] = channels
            caldata['degradation'] = degs
            save_obj(caldata,fname[0:-4])
        else:
            AIA_CalData = load_obj(fname[0:-4])
            degs = AIA_CalData['degradation']
        print('degradation: ',degs)
        Img_dim = np.array(data0.data.shape)
    elif Inst[ind] == 'MaGIXS':
        durs = np.array(data['instruments'][Inst[ind]]['Exposure'])[ch_ind]
        degs = 1
        Img_dim = np.array(data0.shape)
    # Now load in the response functions
    tresp_logt = data['instruments'][Inst[ind]]['Tresp_log10T']
    tresp = data['instruments'][Inst[ind]]['Tresp']
    nt=len(tresp_logt)
    nf=len(channels)    
    trmatrix=np.zeros((nt,nf))
    for i in range(0,nf):
        trmatrix[:,i] = tresp[i]  


    temps = np.logspace(DEM_log[0],DEM_log[1],num=int((DEM_log[1]-DEM_log[0])/DEM_log[2]))  #DEM temperatures grids
    # Temperature bin mid-points for DEM plotting
    mlogt=([np.mean([(np.log10(temps[i])),np.log10((temps[i+1]))]) for i in np.arange(0,len(temps)-1)])   

    plt.close('all')
    fig = plt.figure(figsize=(8,8))
    if Inst[ind] == "AIA":
        ax = plt.subplot(projection=data0)
        data0.plot()
    else:
        ax = plt.subplot()#projection=data0)
        ax.imshow(data0,origin='lower',cmap='hot')
    plt.show()

    DEM_map = np.zeros([Img_dim[0],Img_dim[1],len(mlogt)])
    DEM_map_logTerr = np.zeros([Img_dim[0],Img_dim[1],len(mlogt)])
    DEM_map_err = np.zeros([Img_dim[0],Img_dim[1],len(mlogt)])
    chisq_map = np.zeros([Img_dim[0],Img_dim[1]])
    for i in range(Img_dim[0]):
        for j in range(Img_dim[1]):
            print(i,j)
            data_dem=[]
            for k in range(len(channels)):
                m = data['instruments'][Inst[ind]]['ObsIMG'][channels[k]]
                if Inst[ind] == 'AIA':data_dem.append(m.data[i,j])
                elif Inst[ind] == 'MaGIXS': data_dem.append(m[i,j])
            data_dem=np.array(data_dem) #in The saved images by "EMmap2Image.py" are in the unit of DN/s/px for AIA and XRT
            cor_data=data_dem/degs
            dn_in=cor_data#/durs  #DN/px/s
            #print('dn_in: ',dn_in)
            if all( ii > 0 for ii in dn_in) is True : #if all the chanels has positive values
                if Inst[ind] == 'AIA':
                    #Estimate errors from sunpy
                    num_pix=1
                    edn_in = []
                    for k in range(len(wvn)):
                        aerr_temp=estimate_error((data_dem[k]*durs[k])*(u.ct/u.pix),wvn[k],num_pix) # Will download error table first time using
                        edn_in.append(aerr_temp.value[0])
                    edn_in = np.array(edn_in)
                elif Inst[ind] == 'MaGIXS': 
                    edn_in = np.sqrt((cor_data + (0.2*cor_data))*durs)/durs  
                # 1. Default - reg runs twice, 1st time to work out weight for constraint matrix, then regs with that
                #         Best option if don't know what doing, hence its the default 
                #         Probably best option of AIA data as well.
                dem0,edem0,elogt0,chisq0,dn_reg0=dn2dem_pos(dn_in,edn_in,trmatrix,tresp_logt,temps,gloci=1,emd_int=1) #gloci=0 is default behaviour
                DEM_map[i,j,:] = dem0
                DEM_map_logTerr[i,j,:] = elogt0
                DEM_map_err[i,j,:] = edem0
                chisq_map[i,j] = chisq0
                

                ##  Plot it all
                #fig = plt.figure(figsize=(8, 4.5))
                #plt.errorbar(mlogt,dem0,xerr=elogt0,yerr=edem0,fmt='or',ecolor='lightcoral', \
                #             elinewidth=3, capsize=0,label='Def Self LWght, $\chi^2 =$ {:0.2f}'.format(chisq0))
                #plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
                #plt.ylabel('$\mathrm{DEM\;[cm^{-5}\;K^{-1}]}$')
                #plt.ylim([1e16,4e22])

                #plot_EM_locii = True
                #if plot_EM_locii is True:
                #    for e in range(len(dn_in)):
                #        plt.plot(tresp_logt,EM2DEM(tresp_logt,(dn_in[e]/trmatrix[:,e])))

                #plt.xlim([5.6,7.0])
                #plt.rcParams.update({'font.size': 16})
                #plt.yscale('log')
                #plt.legend()
                ## plt.savefig('demregpy_aiapxl_slw.png',bbox_inches='tight')
                #plt.show()

'''
