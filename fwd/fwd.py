from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob,os
import sunpy
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import rotate,shift
from scipy.io import readsav
from fwd.util import *
from scipy import interpolate
##########

class fwd(object):
    def __init__(self,emcube,emlogt,ref_map):
        self.emcube = emcube
        self.emlogt = emlogt    
        self.ref_map = ref_map

    def create_ovlp(self,img_data,EM_map,nFA,mag_resp,wvl,pix_num,cmap, ref_map,delta_pix_no,vig_func_data,vig_nyFA,mid_vig_fa_ind,mid_vig_fa_ind_y,detector_noise, ph2el,resp_unit=0,exposure=0,magixs_sigma = None, IncludePoisonNoise = False, StoreRawDetector = False, RawDetectorDim = [1040,2152],NonActivePix = [50, 8,4],OutDir = None,hdul_all=[],pointing_info_file = None,SavePlot=False,line_label=[],line_wvl=[],namestr='',NoFrames=1,pix2wvl=None):

        print('====',NoFrames)

        if pointing_info_file is None: raise Exception("%% Error : Provide a valid 'pointing_info_file'")
        FOV,roll,Tx,Ty,botx,boty,topx,topy = np.loadtxt(pointing_info_file,unpack=True)
        #try: 
        #    n__ = len(botx)
        #except:
        #    botx = [botx] ; boty = [boty]; topx = [topx]; topy = [topy]; FOV = [FOV]; roll = [roll] ; Tx = [Tx]; Ty = [Ty]
        for k in range(len(botx)): 
            center_point = ((topx[k] + botx[k]) / 2, (topy[k] + boty[k]) / 2)
            center_point_arcsec = ref_map.pixel_to_world(center_point[0]*u.pix, center_point[1]*u.pix)
            shift_x = img_data.shape[1] / 2 - center_point[0]
            shift_y = img_data.shape[0] / 2 - center_point[1]
   
            # Calculate angle to rotate the image
            delta_x = botx[k] - topx[k]
            delta_y = boty[k] - topy[k] 
            angle = np.degrees(np.arctan2(delta_y, delta_x)) - 90

            ny = len(EM_map[0,:,0])
            EM_new = np.zeros([ny,nFA,len(EM_map[0,0,:])])
            mid_y_ind = ny//2
            for i_ in range(len(EM_map[0,0,:])):
                shifted_em_data = shift(EM_map[:,:,i_], (shift_y, shift_x), mode='constant')
                # Rotate the image
                rotated_em_data = rotate(shifted_em_data, angle, reshape=False,mode='constant')
                center_point = np.array(rotated_em_data.shape)//2
                EM_new[:,:,i_] = rotated_em_data[:,center_point[1]-delta_pix_no : center_point[1]+delta_pix_no+1]
                if vig_func_data != 1: #
                    EM_new[mid_y_ind-vig_nyFA:mid_y_ind+vig_nyFA,:,i_] = EM_new[mid_y_ind-vig_nyFA:mid_y_ind+vig_nyFA,:,i_]*vig_func_data[:,:] #include vigneting
                    EM_new[0:mid_y_ind-vig_nyFA,:,i_] = EM_new[0:mid_y_ind-vig_nyFA,:,i_]*vig_func_data[0,:] #multiply the vigneting edge to the rest data for which vigneting function not-available.
                    EM_new[mid_y_ind+vig_nyFA::,:,i_] = EM_new[mid_y_ind+vig_nyFA::,:,i_]*vig_func_data[-1,:] #multiply the vigneting edge to the rest data
                if magixs_sigma is not None: EM_new[:,:,i_] = ndimage.gaussian_filter(EM_new[:,:,i_], magixs_sigma)#Include` PSF
            shifted_img_data = shift(img_data, (shift_y, shift_x), mode='constant')
            rotated_img_data = rotate(shifted_img_data, angle, reshape=False,mode='constant')
    
            # Display the rotated image
            Orig_IMG = rotated_img_data[:,center_point[1]-delta_pix_no : center_point[1]+delta_pix_no+1]
            '''
            if magixs_sigma is not None: Orig_IMG_vig = ndimage.gaussian_filter(Orig_IMG, magixs_sigma)
            Orig_IMG_vig[mid_y_ind-vig_nyFA:mid_y_ind+vig_nyFA,:] = Orig_IMG[mid_y_ind-vig_nyFA:mid_y_ind+vig_nyFA,:]*vig_func_data
            Orig_IMG_vig[0:mid_y_ind-vig_nyFA,:] = Orig_IMG_vig[0:mid_y_ind-vig_nyFA,:]*vig_func_data[0,:]
            Orig_IMG_vig[mid_y_ind+vig_nyFA::,:] = Orig_IMG_vig[mid_y_ind+vig_nyFA::,:]*vig_func_data[-1,:]
            '''
            #mag_resp #[T, FA, X-pix]
            IMG = np.zeros([ny,len(mag_resp[0,0,:])])
            for iii in range(nFA): #Field angle 
                IMG[:,:] += np.dot(EM_new[:,iii,:],mag_resp[:,iii,:])
            IMG = IMG*exposure
            IMG[IMG < 0] = 0
            ph2el_allframe = np.zeros([ny,len(mag_resp[0,0,:]),NoFrames])
            for d_ in range(ny):
                for d__ in range(NoFrames):
                    ph2el_allframe[d_,:,d__] = ph2el 
            #IMG = np.nan_to_num(IMG)## Replace NaNs with 0 
            IMG_allframe = np.repeat(IMG[:, :, np.newaxis], NoFrames, axis=2)
            if IncludePoisonNoise is True:
                if resp_unit == 'electron':
                    IMG_allframe = IMG_allframe/ph2el_allframe
                    IMG_allframe = np.nan_to_num(IMG_allframe)## Replace NaNs with 0
                    IMG_allframe = np.random.poisson(lam=IMG_allframe) #Add poissin statistics
                    IMG_allframe = IMG_allframe*ph2el_allframe
                else:
                    IMG_allframe = np.random.poisson(lam=IMG_allframe)
            IMG_x = IMG.shape[1]//2
            detector_noise_allframe = np.repeat(detector_noise[:, :, np.newaxis], NoFrames, axis=2)
            if StoreRawDetector is True:
                IMG_Magixs2_raw = np.zeros(RawDetectorDim+[NoFrames]) #[1040,2152]
                extra_y = IMG_Magixs2_raw.shape[0]-IMG.shape[0]
                extra_y = extra_y//2
                #IMG_Magixs2_raw[8:1033,50:1075] = IMG[:,0:1025]
                #IMG_Magixs2_raw[8:1033,1079:2104] = IMG[:,1025::]  #NonActivePix = [50, 8,4,0]
                IMG_Magixs2_raw[extra_y:extra_y+IMG.shape[0],int(NonActivePix[0]):int(IMG_x+NonActivePix[0]),:] = IMG_allframe[:,0:IMG_x,:]
                IMG_Magixs2_raw[extra_y:extra_y+IMG.shape[0],int(IMG_x+NonActivePix[0]+NonActivePix[2]):int(IMG_Magixs2_raw.shape[1]-NonActivePix[0]),:] = IMG_allframe[:,IMG_x::,:]
                if len(detector_noise) > 0: IMG_Magixs2_raw = IMG_Magixs2_raw+detector_noise_allframe #include detector noise to the raw image
                #display_rotated_data(IMG_Magixs2_raw,fig, ax3,cmap=cmap)
            if len(detector_noise) > 0:
                
                detector_noise_IMG_fov = np.concatenate((detector_noise_allframe[int(NonActivePix[1]):-int(NonActivePix[1]),int(NonActivePix[0]):int(IMG_x+NonActivePix[0]),:], detector_noise_allframe[int(NonActivePix[1]):-int(NonActivePix[1]),int(IMG_x+NonActivePix[0]+NonActivePix[2]):int(RawDetectorDim[1]-NonActivePix[0]),:]), axis = 1) #remove non-active pixcels
                extra_y = detector_noise_IMG_fov.shape[0]-IMG.shape[0]
                extra_y = extra_y//2
                detector_noise_IMG_fov = detector_noise_IMG_fov[extra_y:extra_y+IMG.shape[0],:,:]#Adjust the FOV in cross-dispersion director with FOV of IMG 
                if resp_unit == 'electron':
                    IMG_allframe = np.add(IMG_allframe,detector_noise_IMG_fov,out=IMG_allframe, casting="unsafe")
                else:
                    IMG_allframe = np.add(IMG_allframe,detector_noise_IMG_fov/ph2el_allframe,out=IMG_allframe, casting="unsafe")

            if SavePlot is True:
                fig,ax,ax2,ax3,ax4 = GUI(img_data,cmap,ref_map)

                line = plt.Line2D((botx[k], topx[k]),
                                  (boty[k], topy[k]),
                                  color='white')
                ax.add_line(line)
                fig.canvas.draw()
                display_rotated_data(Orig_IMG, fig, ax2,cmap=cmap)
                display_rotated_data(IMG_allframe[:,:,0],fig, ax3,cmap=cmap)
                fname = namestr+'FOV_'+format('%d'%FOV[k])
                fig.canvas.mpl_connect('motion_notify_event',lambda event: on_motion_cursor_2(event,IMG_allframe[:,:,0],fig, wvl,pix_num,line_label,line_wvl,resp_unit,OutDir=OutDir,namestr=fname))
                plt.show()

            if OutDir is not None:
                if FOV[k] == 1:
                    hdul_all[0].append(fits.ImageHDU(pix2wvl,name='pix2wvl'))
                    hdul_all[0].append(fits.ImageHDU(ph2el,name='ph2el'))
                    hdul_all[0][0].header.append(('exposere',exposure))
                    hdul_all[0][0].header.append(('resp_unit',resp_unit))
                hdul_all[0].append(fits.ImageHDU(IMG_allframe.transpose(2, 0, 1),name='FOV-'+format('%d'%FOV[k]))) #Data is located from 3rd extension
                hdul_all[0][int(FOV[k])+1].header.append(('Roll', angle, 'degrees'))  
                hdul_all[0][int(FOV[k])+1].header.append(('Tx', center_point_arcsec.Tx.value,'arcsec'))
                hdul_all[0][int(FOV[k])+1].header.append(('Ty', center_point_arcsec.Ty.value,'arcsec'))
                hdul_all[0][int(FOV[k])+1].header.append(('exposere',exposure,'second'))
                hdul_all[0][int(FOV[k])+1].header.append(('resp_unit',resp_unit,))
                hdul_all[0][int(FOV[k])+1].header.append(('LOS_p1',str([topx[k],topy[k]]),'pix'))
                hdul_all[0][int(FOV[k])+1].header.append(('LOS_p2',str([botx[k],boty[k]]),'pix'))   
 
                hdul_all[1].append(fits.ImageHDU(Orig_IMG,name='FOV-'+format('%d'%FOV[k])))
                hdul_all[1][int(FOV[k])-1].header.append(('Roll', angle, 'degrees'))
                hdul_all[1][int(FOV[k])-1].header.append(('Tx', center_point_arcsec.Tx.value,'arcsec'))
                hdul_all[1][int(FOV[k])-1].header.append(('Ty', center_point_arcsec.Ty.value,'arcsec'))
    
                hdul_all[3].append(fits.ImageHDU(EM_new.transpose(2, 0, 1),name='FOV-'+format('%d'%FOV[k]))) #Data is located from 2nd extension
                hdul_all[3][int(FOV[k])].header.append(('Roll', angle, 'degrees'))
                hdul_all[3][int(FOV[k])].header.append(('Tx', center_point_arcsec.Tx.value,'arcsec'))
                hdul_all[3][int(FOV[k])].header.append(('Ty', center_point_arcsec.Ty.value,'arcsec'))
    
                hdul_all[2].write('%d\t'%(FOV[k]))
                hdul_all[2].write('%0.3f\t'%angle)
                hdul_all[2].write('%0.3f\t'%center_point_arcsec.Tx.value)
                hdul_all[2].write('%0.3f\t'%center_point_arcsec.Ty.value)
                hdul_all[2].write('%0.3f\t'%botx[k])
                hdul_all[2].write('%0.3f\t'%boty[k])
                hdul_all[2].write('%0.3f\t'%topx[k])
                hdul_all[2].write('%0.3f\n'%topy[k])
    
                if StoreRawDetector is True:
                    hdul_all[4].append(fits.ImageHDU(IMG_Magixs2_raw.transpose(2, 0, 1),name='FOV-'+format('%d'%FOV[k])))
                    hdul_all[4][int(FOV[k])-1].header.append(('Roll', angle, 'degrees'))
                    hdul_all[4][int(FOV[k])-1].header.append(('Tx', center_point_arcsec.Tx.value,'arcsec'))
                    hdul_all[4][int(FOV[k])-1].header.append(('Ty', center_point_arcsec.Ty.value,'arcsec'))
                    hdul_all[4][int(FOV[k])-1].header.append(('exposere',exposure,'second'))
                    hdul_all[4][int(FOV[k])-1].header.append(('resp_unit',resp_unit))
                    hdul_all[4][int(FOV[k])-1].header.append(('LOS_p1',str([topx[k],topy[k]]),'pix'))
                    hdul_all[4][int(FOV[k])-1].header.append(('LOS_p2',str([botx[k],boty[k]]),'pix'))
            print('Completed')



    def magixs_fwd(self, FullSun_img_data, resp_file, pix2wvl = None,resp_unit = None, plate_scale = None,Vigneting_func_file=None,psf=None,exposure = 1, line_label = [],line_wvl = [],cmap='hot',StoreRawDetector = False,RawDetectorDim = [1040,2152],NonActivePix = [50, 8,4], IncludePoisonNoise=False,IncludeDetecorNoise=False,detector_noise = [],phot2elec_conv_file = None,StoreOutputs = False,OutDir=None,namestr='MaGIXS2_',pointing_info_file=None,SavePlot = False,NoFrames=1):
        '''

        Purpose:
            calculate forword modeled detector using response file and an EM cube.

        Inputs-
            FullSun_img_data : Sunpy map object. For ploting purpose.

            resp_file : fits file.
                primaty hdu - image_array [wvl,FA,logT]
                secondary hdu - table_array [INDEX,LogT]
                third hdu - table_array [INDEX, FIELD_ANGLE]

            resp_unit : = 'electron' or 'photon' -> unit of the response file.

            pix2wvl : 1D array of onaxis pixcel to wvl relation, for ploting purpose.

            plate_scale :
 
        Optional Inputs-
            Vigneting_func_file : 
            psf :
            exposure : 
            line_label : 
            line_wvl : 
            plate_scale : 
            cmap : 
            StoreRawDetector : 
            RawDetectorDim : [X,Y]
            NonActivePix : [X,Y,midX] #Non-active pixcels is raw detector: X,Y: begining and end of X and Y directions. midX: middle of X directions.  
            IncludePoisonNoise : 
            IncludeDetecorNoise : 
            detector_noise : 
            phot2elec_conv_file : 
            StoreOutputs : 
            OutDir : 
            namestr : 
            NoFrames : 
          
        '''
        ref_map = self.ref_map
        if resp_unit is None: raise Exception("%% Error : 'resp_unit' should be in 'electron' or 'photon'")
        if len(line_label) != len(line_wvl): raise Exception("%% Error : dimensions of 'line_label' and 'line_wvl' should be same.")
        if plate_scale is None: raise Exception("%% Error : Provide the plate scale of the input data.")

        ph2el = 1
        if IncludeDetecorNoise is True: 
            #if detector_noise_file is None: raise Exception("%% Error : Provide a valid detector noise data file")
            #detector_noise = readsav(detector_noise_file)['data'][2,:,:] #in the unit of electron or photon
            if phot2elec_conv_file is None: raise Exception("%% Error : Provide a valid 'phot2elec_conv_file'")
            pix_noise_ind,ph2el = np.loadtxt(phot2elec_conv_file,unpack=True)
        hdul_all = []
        if StoreOutputs is True:
            if OutDir is None: 
                os.mkdir('outputs')
                OutDir = 'outputs'
            hdul = fits.HDUList()
            hdul_ref = fits.HDUList()
            pointing_info = open(os.path.join(OutDir,'pointing_info.dat'),'w')
            pointing_info.write('# FOV, Roll (deg), Solar-X (arc-sec), Solar-Y (arc-sec), bot-X, bot-Y, top-X, top-Y\n')
            hdul_emIn = fits.HDUList()
            hdul_emIn.append(fits.ImageHDU(self.emlogt,name='log10T'))
            hdul_all += [hdul,hdul_ref,pointing_info,hdul_emIn]
            if StoreRawDetector is True:
                hdul_raw = fits.HDUList()
                hdul_all += [hdul_raw]
            '''
            if pointing_info_file is None:
                pointing_info_raw = open(os.path.join(OutDir,'pointing_info_raw.dat'),'w')
                pointing_info_raw.write('# FOV, bot-X, bot-Y, top-X, top-Y \n')
                hdul_all += [pointing_info_raw]
            '''
        else: OutDir = None


        magixs_sigma = (psf/plate_scale)/2.355 #in pixcel unit. As psf is in arc-sec

        with fits.open(resp_file) as hdul:
            nloops = len(hdul)
            mag_resp = hdul[0].data #[T, FA, X-pix]
            mag_FA = hdul[2].data
            mag_resp_logtemp = hdul[1].data
            resp_header = hdul[0].header
        mag_FA = mag_FA.FIELD_ANGLE
        mag_resp_logtemp = mag_resp_logtemp.LOGT
        nFA = len(mag_FA)
        delta_pix_no = (nFA-1)//2
        if len(mag_resp_logtemp) != len(self.emlogt): raise Exception("%% Error : Inconsistent dimmension logT in respose and emcube")
        elif max(abs(mag_resp_logtemp - self.emlogt)) > 0.01 : raise Exception("%% Error : Inconsistent logT in respose and emcube")

        ''' 
        wvl_resp_data = readsav(Wvl_cal_file)
        ind_onaxis = np.where(wvl_resp_data['fieldangles'] == 0)
        wvl = wvl_resp_data['wavearray'][ind_onaxis,:][0,0,:] #Using On-axis wavelength response
        '''
        if pix2wvl is not None:
            wvl = pix2wvl
            pix_num = np.arange(len(wvl))
        else:
            pix_num = np.arange(len(mag_resp[0,0,:])) 
            wvl = pix_num
       

        if Vigneting_func_file is not None:
            vig_func = readsav(Vigneting_func_file)
            vig_func_data = vig_func['nvigmap']
            vig_FA_x = vig_func['xpix_fa']
            vig_FA_y = vig_func['ypix_fa']
            vig_nyFA = int(len(vig_FA_y)/2)
            mid_vig_fa_ind = np.where(vig_func['xpix_fa'] == 0)[0][0]
            vig_func_data = vig_func_data[:,mid_vig_fa_ind-delta_pix_no:mid_vig_fa_ind+delta_pix_no+1]
            vig_FA_x = vig_FA_x[mid_vig_fa_ind-delta_pix_no:mid_vig_fa_ind+delta_pix_no+1]
            mid_vig_fa_ind_y = np.where(vig_func['ypix_fa'] == 0)[0][0]
            
            vig_func_data = np.array(vig_func_data)
            #vig_func_data[100:300,:] = 0 #check if its working
        else: 
            vig_func_data = 1
            vig_nyFA = 0
            mid_vig_fa_ind = 0
            mid_vig_fa_ind_y = 0
              

        if pointing_info_file is None:

            fig,ax,ax2,ax3,ax4 = GUI(FullSun_img_data,cmap,ref_map)
 
            # Connect mouse events to the callback functions
            fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event,FullSun_img_data,self.emcube,nFA,mag_resp,cmap, ref_map, fig,delta_pix_no, vig_func_data,vig_nyFA,mid_vig_fa_ind,mid_vig_fa_ind_y,detector_noise, ph2el,resp_unit=resp_unit,exposure=exposure,magixs_sigma = magixs_sigma,IncludePoisonNoise=IncludePoisonNoise,StoreRawDetector=StoreRawDetector,RawDetectorDim=RawDetectorDim,NonActivePix = NonActivePix,OutDir=OutDir,hdul_all=hdul_all,pix2wvl = wvl,NoFrames = NoFrames))
            fig.canvas.mpl_connect('motion_notify_event', lambda event: on_motion(event,fig))
            fig.canvas.mpl_connect('motion_notify_event',lambda event: on_motion_cursor(event,fig, wvl,pix_num,line_label,line_wvl,resp_unit,OutDir=OutDir,namestr=namestr))
            
            plt.show()

        else:
            if SavePlot is False: plt.close('all') 
            self.create_ovlp(FullSun_img_data,self.emcube,nFA,mag_resp,wvl,pix_num,cmap, ref_map,delta_pix_no,vig_func_data,vig_nyFA,mid_vig_fa_ind,mid_vig_fa_ind_y,detector_noise, ph2el,resp_unit=resp_unit,exposure=exposure,magixs_sigma = magixs_sigma, IncludePoisonNoise = IncludePoisonNoise, StoreRawDetector = StoreRawDetector, RawDetectorDim = RawDetectorDim,NonActivePix = NonActivePix,OutDir = OutDir,hdul_all=hdul_all,pointing_info_file = pointing_info_file,SavePlot=SavePlot,line_label=line_label,line_wvl=line_wvl,namestr=namestr,NoFrames=NoFrames,pix2wvl=wvl)

        if StoreOutputs is True:
            hdul_all[0].writeto(os.path.join(OutDir,namestr+'FWDmodel.fits'),overwrite=True)
        
            hdul_all[1].writeto(os.path.join(OutDir,namestr+'FOV.fits'),overwrite=True)
        
            hdul_all[3].writeto(os.path.join(OutDir,namestr+'EMcubeInput.fits'),overwrite=True)
        
            hdul_all[2].close()
        
            if StoreRawDetector is True:
                hdul_all[4].writeto(os.path.join(OutDir,namestr+'FWDmodel_rawDetector.fits'),overwrite=True)
            hdul_all[-1].close()
    

    def make_average_detector(self,Overlappogram_file,OutFile):
        hdul = fits.open(Overlappogram_file)
        nframe = len(hdul) - 2
        for i in range(nframe):
            data = hdul[i+2].data
            average_data = np.average(data,axis=0)
            hdul[i+2].data = average_data
        hdul.writeto(OutFile, overwrite=True)    

    def run_inversion(self,dataDir,overlappogram_file, resp_file, gnt_file,configfile, weight_file = None, OutDir=None, prefix = "MaGIXS2_inversion",FOV=1):
        #search_and_replace_line(configfile, 'weights','weights='+os.path.join(dataDir,weight_file))
        search_and_replace_line(configfile, 'response','response="'+os.path.join(resp_file)+'"')
        search_and_replace_line(configfile, 'gnt','gnt="'+os.path.join(gnt_file)+'"')

        if OutDir is None:
            OutDir = os.path.join(dataDir,'inversion')
            if os.path.isdir(OutDir) is False: os.mkdir(OutDir)
        search_and_replace_line(configfile, 'directory', 'directory="'+OutDir+'"')
        
        search_and_replace_line(configfile, 'prefix','prefix="'+prefix+'"')
        
        spec_pure_hdul = fits.HDUList()

        hdul = fits.open(os.path.join(dataDir,overlappogram_file))
        ph2el = hdul[1].data
        if FOV is None: #Consider all FOV
            hdul = fits.open(overlappogram_file)
            nFOV = len(hdul)
            new_hdul = fits.HDUList([primary_hdu] + [hdul[i] for i in range(len(hdul)) if i != 0 and i != 1])
        else:
            overlappogram = hdul[FOV+1]
            primary_hdu = fits.PrimaryHDU(data=(overlappogram.data / float(overlappogram.header['EXPOSERE'])), header=overlappogram.header)
            new_hdul = fits.HDUList(primary_hdu) 
            temp_file_path = os.path.join(OutDir,'InputOverlapogram.fits')
            new_hdul.writeto(temp_file_path, overwrite=True)
            search_and_replace_line(configfile, 'overlappogram', 'overlappogram="'+os.path.join(temp_file_path)+'"')

            if weight_file is None:
                '''
                IMG_x = overlappogram.data.shape[1]//2
                detector_noise_IMG_fov = np.concatenate((detector_noise[int(NonActivePix[1]):-int(NonActivePix[1]),
                              int(NonActivePix[0]):int(IMG_x+NonActivePix[0]),:], detector_noise[int(NonActivePix[1]):-int(NonActivePix[1]),
                              int(IMG_x+NonActivePix[0]+NonActivePix[2]):int(detector_noise.shape[1]-NonActivePix[0]),:]), 
                              axis = 1) #remove non-active pixcels
                extra_y = detector_noise_IMG_fov.shape[0]-overlappogram.data.shape[0]
                extra_y = extra_y//2
                detector_noise_IMG_fov = detector_noise_IMG_fov[extra_y:extra_y+overlappogram.data.shape[0],:,:]#Adjust the FOV in cross-dispersion director with FOV of IMG 
                '''
                error = np.sqrt(overlappogram.data/ph2el) / float(overlappogram.header['EXPOSERE'])#ph
                sys_error = (overlappogram.data/ph2el)*0.05 / float(overlappogram.header['EXPOSERE'])
                error = np.sqrt(error*2 + sys_error**2)*ph2el
                error = 1/error
                error = np.nan_to_num(error)
                primary_hdu2 = fits.PrimaryHDU(data=error, header=overlappogram.header)
                new_hdul2 = fits.HDUList(primary_hdu2)
                weight_file = os.path.join(OutDir,'Weights.fits')
                new_hdul2.writeto(weight_file, overwrite=True)
            
            search_and_replace_line(configfile, 'weights','weights="'+weight_file+'"')
            os.system('unfold '+configfile)

            hdul_spmap = fits.open(os.path.join(OutDir,prefix+'_spectral_x2_5.0_5e-05_wpsf.fits'))
            spec_pure_hdul.append(fits.ImageHDU(hdul_spmap[1].data,name='IonList'))
            spec_pure_hdul.append(fits.ImageHDU(hdul_spmap[2].data,name='logT'))
            spec_pure_hdul.append(fits.ImageHDU(hdul_spmap[0].data,name='FOV'+format('%d'%FOV)))
            specpure_outfile = os.path.join(OutDir,prefix+'SpecPure.fits')
            spec_pure_hdul.writeto(temp_file_path3, overwrite=True)
            #os.remove(temp_file_path)
            #os.remove(weight_file)

    def plot_detector(self,ref_map,Overlappogram_file,EMin_file,cmap='hot',line_label=[],line_wvl=[],FOV = 1, LOS_line_color = 'w',FrameRange = [0,1],unit='electron'):

        FullSun_img_data = np.nan_to_num(ref_map.data)
        if FOV < 1: FOV = 1
        with fits.open(EMin_file) as hdul:
            logt = hdul[FOV-1].data
            EM = hdul[FOV].data
            summedEM = np.sum(EM[:,:,:],axis=0)
        with fits.open(Overlappogram_file) as hdul:
            pix2wvl = hdul[0].data
            ph2el = hdul[1].data
            data = hdul[FOV+1].data
            resp_unit = hdul[0].header['RESP_UNIT']
            exposure = hdul[FOV+1].header['EXPOSERE']
            roll = hdul[FOV+1].header['ROLL']
            Tx = hdul[FOV+1].header['Tx']
            Ty = hdul[FOV+1].header['Ty']
            start_point = hdul[FOV+1].header['LOS_P1'].split(',')
            end_point = hdul[FOV+1].header['LOS_P2'].split(',')
        start_point = [float(start_point[0][1::]),float(start_point[1][0:-1])]
        end_point = [float(end_point[0][1::]),float(end_point[1][0:-1])]
        #center_point = ref_map.world_to_pixel(SkyCoord(Tx*u.arcsec,Ty*u.arcsec,frame=ref_map.coordinate_frame))
        #start_point,end_point = get_line_start_end(-roll, center_point.x.value,center_point.y.value,length = linelength)
        if FrameRange is None:
            FrameRange = range(0,len(data[:,0,0]))
        else: FrameRange = range(FrameRange[0],FrameRange[1])
        for l in FrameRange:
            fig,ax,ax2,ax3,ax4 = GUI(FullSun_img_data,cmap,ref_map)
            line = plt.Line2D((start_point[0], end_point[0]),
                          (start_point[1], end_point[1]),
                          color=LOS_line_color,lw=2,ls='--')
            ax.add_line(line)
            fig.canvas.draw()
            display_rotated_data(summedEM, fig, ax2,cmap=cmap)
            pix_num = np.arange(len(pix2wvl))
            if resp_unit == 'electron':
                if unit != 'electron':
                    data2 = data[l,:,:] / (ph2el*exposure)
                    resp_unit = 'photon'
            elif resp_unit == 'photon':
                if unit != 'photon':
                    data2 = data[l,:,:]*ph2el/exposure
                    resp_unit = 'electron'
            display_rotated_data(data2,fig, ax3,cmap=cmap) 
            fig.canvas.mpl_connect('motion_notify_event',lambda event: on_motion_cursor_2(event,data2,fig, pix2wvl,pix_num,line_label,line_wvl,resp_unit,OutDir='',namestr=''))        
            plt.show()
 

    def plot_detector_avg(self,ref_map,Overlappogram_file,EMin_file,SpecPure_file=None,EffAreaFile = None,cmap='hot',line_label=[],line_wvl=[],FOV = 1,specInd = [7,12,5],Pscale = [2.8*2.8,2.8],unit='electron',LOS_line_color='w'):
        
        FullSun_img_data = np.nan_to_num(ref_map.data)
        if FOV < 1: FOV = 1 
        with fits.open(EMin_file) as hdul:
            logt = hdul[FOV-1].data
            EM = hdul[FOV].data
            summedEM = np.sum(EM[:,:,:],axis=0)
        with fits.open(Overlappogram_file) as hdul:
            pix2wvl = hdul[0].data
            ph2el = hdul[1].data
            data = hdul[FOV+1].data
            resp_unit = hdul[0].header['RESP_UNIT']
            exposure = hdul[FOV+1].header['EXPOSERE']
            roll = hdul[FOV+1].header['ROLL']
            Tx = hdul[FOV+1].header['Tx']
            Ty = hdul[FOV+1].header['Ty']
            start_point = hdul[FOV+1].header['LOS_P1'].split(',')
            end_point = hdul[FOV+1].header['LOS_P2'].split(',')
        start_point = [float(start_point[0][1::]),float(start_point[1][0:-1])]
        end_point = [float(end_point[0][1::]),float(end_point[1][0:-1])]
        #center_point = ref_map.world_to_pixel(SkyCoord(Tx*u.arcsec,Ty*u.arcsec,frame=ref_map.coordinate_frame))
        #start_point,end_point = get_line_start_end(-roll, center_point.x.value,center_point.y.value,length = linelength)
        if SpecPure_file is not None:
            with fits.open(SpecPure_file) as hdul:
                dataSpec = hdul[0].data
                ions = hdul[2].data
            Au2Km = 1.496e8 #Au to Km
            f = (Pscale[0]*Pscale[1]*725*725)/(Au2Km**2)
            if EffAreaFile is not None:
                wvl_,effarea = np.loadtxt(EffAreaFile,unpack=True)
                intpf = interpolate.interp1d(wvl_, effarea)

                for w in range(len(specInd)):
                    ion_wvl = float(ions[specInd[w]][1].split('@')[1])
                    EA = intpf(ion_wvl)
                    unit = 'Ph s$^{-1}$'
            else: EA = 1 ; unit = 'Ph cm$^{-2}$ s$^{-1}$' 
            dataSpec = dataSpec*f*EA
        fig,ax,ax2,ax3,ax4,ax5,ax6,ax7 = GUI_plot(FullSun_img_data,cmap,ref_map)
        line = plt.Line2D((start_point[0], end_point[0]),
                          (start_point[1], end_point[1]),
                          color=LOS_line_color,lw=2,ls='--')
        ax.add_line(line)
        fig.canvas.draw()

        if SpecPure_file is not None:
            ax5.imshow(dataSpec[specInd[0],:,:],origin='lower',cmap=cmap)
            add_colorbar(ax5, dataSpec[specInd[0],:,:], cmap, unit)

            ax6.imshow(dataSpec[specInd[1],:,:],origin='lower',cmap=cmap)
            add_colorbar(ax6, dataSpec[specInd[1],:,:], cmap, unit)

            ax7.imshow(dataSpec[specInd[2],:,:],origin='lower',cmap=cmap)
            add_colorbar(ax7, dataSpec[specInd[2],:,:], cmap, unit)

            ax5.imshow(dataSpec[specInd[0],:,:],origin='lower',cmap=cmap)
            ax5.text(0.18,1.01,ions[specInd[0]][1],transform = ax5.transAxes)
            ax6.imshow(dataSpec[specInd[1],:,:],origin='lower',cmap=cmap)
            ax6.text(0.18,1.01,ions[specInd[1]][1],transform = ax6.transAxes)
            ax7.imshow(dataSpec[specInd[2],:,:],origin='lower',cmap=cmap)
            ax7.text(0.18,1.01,ions[specInd[2]][1],transform = ax7.transAxes)
        display_rotated_data(summedEM, fig, ax2,cmap=cmap)
        pix_num = np.arange(len(pix2wvl))
        if resp_unit == 'electron':
            if unit != 'electron':
                data = data / (ph2el*exposure)
                resp_unit = 'photon'
        elif resp_unit == 'photon':
            if unit != 'photon':
                data = data*ph2el/exposure
                resp_unit = 'electron'
        display_rotated_data(data,fig, ax3,cmap=cmap) 
        fig.canvas.mpl_connect('motion_notify_event',lambda event: on_motion_cursor_2(event,data,fig, pix2wvl,pix_num,line_label,line_wvl,resp_unit,OutDir=OutDir,namestr=''))                    
        plt.show()


