from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob,os
import sunpy
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import rotate,shift
from scipy.io import readsav
from scipy import ndimage
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable

#=== Function to save/load python ductionary ===
def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

# Initialize global variables
start_point = None
end_point = None
line = None
img_data = None
rotated_em_data = None
rotated_img_data = None
FOV=0
# Mouse motion event callback function
def on_motion(event,fig):
    global start_point, end_point, line
    axes = fig.get_axes()
    ax = axes[0]
    if event.inaxes == ax and start_point is not None and end_point is None:
        if line is not None:
            line.remove()
        line = plt.Line2D((start_point[0], event.xdata),
                          (start_point[1], event.ydata),
                          color='white')
        ax.add_line(line)
        fig.canvas.draw()
'''
def pix2wvl(x):
    return np.interp(x,pix_num,wvl)
def wvl2pix(x):
    return np.interp(x,wvl,pix_num)
'''
IMG_frame0 = []

# Mouse motion event callback function
def on_motion_cursor(event,fig, wvl,pix_num,line_label,line_wvl,resp_unit,OutDir=None,namestr='',):
    global cursor_position, IMG_frame0, FOV#,ax4, wvl

    def pix2wvl(x):
        return np.interp(x,pix_num,wvl)
    def wvl2pix(x):
        return np.interp(x,wvl,pix_num)

    axes = fig.get_axes()
    ax3 = axes[2]
    ax4 = axes[3]
    #if event.inaxes is not None:
    if event.inaxes == ax3:
        x = event.xdata ; y = event.ydata
        cursor_position = (x, y)
        #ax3.hline(y=int(event.ydata),xmin=0,xmax=1,color='w')
        line = plt.Line2D((0, len(IMG_frame0[0,:])),
                          (int(y), int(event.ydata)),
                          color='white')
        ax3.add_line(line)
        fig.canvas.draw()
        #print(f"Cursor position: {cursor_position}")
        ax4.clear()
        #ax4_second = ax4.secondary_xaxis('top',function(pix2wvl,wvl2pix))
        #spec = np.average(IMG[int(y)-1:int(y)+2,:],axis=0)
        spec = IMG_frame0[int(y),:]
        #specp = ax4.plot(np.arange(0,len(IMG[0,:])),spec,'b')
        specp = ax4.plot(wvl,spec,'b')
        if resp_unit == 'electron':ax4.set_ylabel('el s$^{-1}$ pix$^{-1}$')
        if resp_unit == 'photon': ax4.set_ylabel('Ph s$^{-1}$ pix$^{-1}$')
        ax4_second = ax4.secondary_xaxis('top',functions=(wvl2pix,pix2wvl))
        ax4_second.set_xlabel('Pixcel number')
        ax4.set_xlabel('Wavelength (${\AA}$)')
        ax4.set_xlim([wvl[0],wvl[-1]])
        #ax4.set_ylim([0,ax4.ylim()[1]])
        #primary_ticks = wvl2pix(ax4.get_xticks()[0:-1])
        #ax4_second.set_xticks(primary_ticks)
        #leb = [f'{wvl2pix(tick):.1f}' for tick in ax4.get_xticks()[0:-1]]
        #ax4_second.set_xticklabels(leb)
        ylim = ax4.get_ylim()
        for w in range(len(line_wvl)):
            x_ = line_wvl[w]
            ax4.text(x_,ylim[1]*0.78,line_label[w],fontsize=8,rotation=90,verticalalignment = 'bottom',horizontalalignment = 'center')
            ax4.plot([x_,x_],[0,ylim[1]*0.78],'--r',lw=1,alpha=0.6)
        fname = namestr+'FOV_'+format('%d'%FOV)
        #extent = ax4.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if OutDir is not None: fig.savefig(os.path.join(OutDir,fname+'.png'))#,bbox_inches=extent)
        fig.canvas.draw()
        if line is not None:
            line.remove()


def finalize_line(img_data,EM_map,nFA,mag_resp,cmap, ref_map,fig,delta_pix_no,vig_func_data,vig_nyFA,mid_vig_fa_ind,mid_vig_fa_ind_y,detector_noise, ph2el,resp_unit=0,exposure=0,magixs_sigma=None,IncludePoisonNoise=False,StoreRawDetector=False,RawDetectorDim = [1040,2152],NonActivePix = [50, 8,4],OutDir=None,hdul_all=[],NoFrames=1,pix2wvl=None):
    global start_point, end_point, line,IMG_frame0
    axes = fig.get_axes()
    ax = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]
    if start_point is not None and end_point is not None:
        if line is not None:
            line.remove()
        line = plt.Line2D((start_point[0], end_point[0]),
                          (start_point[1], end_point[1]),
                          color='white')
        ax.add_line(line)
        fig.canvas.draw()

        # Calculate center point of the line
        center_point = ((start_point[0] + end_point[0]) / 2, (start_point[1] + end_point[1]) / 2)
        #print(f"Center point: {center_point}")

        shift_x = img_data.shape[1] / 2 - center_point[0]
        shift_y = img_data.shape[0] / 2 - center_point[1]

        # Calculate angle to rotate the image
        delta_x = end_point[0] - start_point[0]
        delta_y = end_point[1] - start_point[1]
        angle = np.degrees(np.arctan2(delta_y, delta_x)) - 90
        center_point_arcsec = ref_map.pixel_to_world(center_point[0]*u.pix, center_point[1]*u.pix)

        ny = len(EM_map[0,:,0])
        EM_new = np.zeros([ny,nFA,len(EM_map[0,0,:])])
        mid_y_ind = ny//2
        for i_ in range(len(EM_map[0,0,:])):
            shifted_em_data = shift(EM_map[:,:,i_], (shift_y, shift_x), mode='constant')
            # Rotate the image
            rotated_em_data = rotate(shifted_em_data, angle, reshape=False,mode='constant')
            center_point = np.array(rotated_em_data.shape)//2
            EM_new[:,:,i_] = rotated_em_data[:,center_point[1]-delta_pix_no : center_point[1]+delta_pix_no+1]
            EM_new[mid_y_ind-vig_nyFA:mid_y_ind+vig_nyFA,:,i_] = EM_new[mid_y_ind-vig_nyFA:mid_y_ind+vig_nyFA,:,i_]*vig_func_data[:,:] #include vigneting
            EM_new[0:mid_y_ind-vig_nyFA,:,i_] = EM_new[0:mid_y_ind-vig_nyFA,:,i_]*vig_func_data[0,:] #multiply the vigneting edge to the rest data for which vigneting function not-available.
            EM_new[mid_y_ind+vig_nyFA::,:,i_] = EM_new[mid_y_ind+vig_nyFA::,:,i_]*vig_func_data[-1,:] #multiply the vigneting edge to the rest data
            if magixs_sigma is not None: EM_new[:,:,i_] = ndimage.gaussian_filter(EM_new[:,:,i_], magixs_sigma)#Include` PSF

        shifted_img_data = shift(img_data, (shift_y, shift_x), mode='constant')
        rotated_img_data = rotate(shifted_img_data, angle, reshape=False,mode='constant')

        # Display the rotated image
        Orig_IMG = rotated_img_data[:,center_point[1]-delta_pix_no : center_point[1]+delta_pix_no+1]
        display_rotated_data(Orig_IMG, fig, ax2,cmap=cmap)
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
            IMG_Magixs2_raw = np.zeros(RawDetectorDim+[NoFrames]) #[1040,2152,NoFrame]
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
        IMG_frame0 = IMG_allframe[:,:,0]
        display_rotated_data(IMG_allframe[:,:,0],fig, ax3,cmap=cmap)
        if OutDir is not None:
            if FOV == 0:
                '''
                c1 = fits.Column(name='pix2wvl',array=pix2wvl,format='D')
                c2 = fits.Column(name='ph2el',array=ph2el,format='D')
                hdul_all[0].append(fits.TableHDU.from_columns([c1, c2],name='par'))
                '''
                hdul_all[0].append(fits.ImageHDU(pix2wvl,name='pix2wvl'))
                hdul_all[0].append(fits.ImageHDU(ph2el,name='ph2el'))
                #hdul_all[0][int(FOV)-1].header.append(('Roll', angle))
                #hdul_all[0][int(FOV)-1].header.append(('Tx', center_point_arcsec.Tx.value))
                #hdul_all[0][int(FOV)-1].header.append(('Ty', center_point_arcsec.Ty.value))
                #hdul_all[4][int(FOV)-1].header.append(('resp_unit',resp_unit))


            hdul_all[0].append(fits.ImageHDU(IMG_allframe.transpose(2, 0, 1),name='FOV-'+format('%d'%FOV)))
            hdul_all[0][int(FOV)-1].header.append(('Roll', angle))
            hdul_all[0][int(FOV)-1].header.append(('Tx', center_point_arcsec.Tx.value))
            hdul_all[0][int(FOV)-1].header.append(('Ty', center_point_arcsec.Ty.value))
            hdul_all[0][int(FOV)-1].header.append(('resp_unit',resp_unit))
            
            hdul_all[1].append(fits.ImageHDU(Orig_IMG,name='FOV-'+format('%d'%FOV)))
            hdul_all[1][int(FOV)-1].header.append(('Roll', angle))
            hdul_all[1][int(FOV)-1].header.append(('Tx', center_point_arcsec.Tx.value))
            hdul_all[1][int(FOV)-1].header.append(('Ty', center_point_arcsec.Ty.value))

            hdul_all[3].append(fits.ImageHDU(EM_new.transpose(2, 0, 1),name='FOV-'+format('%d'%FOV)))
            hdul_all[3][int(FOV)-1].header.append(('Roll', angle))
            hdul_all[3][int(FOV)-1].header.append(('Tx', center_point_arcsec.Tx.value))
            hdul_all[3][int(FOV)-1].header.append(('Ty', center_point_arcsec.Ty.value))

            hdul_all[2].write('%d\t'%(FOV+1))
            hdul_all[2].write('%0.3f\t'%angle)
            hdul_all[2].write('%0.3f\t'%center_point_arcsec.Tx.value)
            hdul_all[2].write('%0.3f\n'%center_point_arcsec.Ty.value)

            if StoreRawDetector is True:
                hdul_all[4].append(fits.ImageHDU(IMG_Magixs2_raw.transpose(2, 0, 1),name='FOV-'+format('%d'%FOV)))
                hdul_all[4][int(FOV)-1].header.append(('Roll', angle))
                hdul_all[4][int(FOV)-1].header.append(('Tx', center_point_arcsec.Tx.value))
                hdul_all[4][int(FOV)-1].header.append(('Ty', center_point_arcsec.Ty.value))
                hdul_all[4][int(FOV)-1].header.append(('resp_unit',resp_unit))
            hdul_all[-1].write('%d\t'%(FOV+1))
            hdul_all[-1].write('%0.3f\t'%end_point[0]) 
            hdul_all[-1].write('%0.3f\t'%end_point[1])
            hdul_all[-1].write('%0.3f\t'%start_point[0])
            hdul_all[-1].write('%0.3f\n'%start_point[1])
            
        '''
        if event.inaxes is not None:
            cursor_position = (event.xdata, event.ydata)
            print(f"Cursor position: {cursor_position}")
        '''
        # Reset points for new selection
        start_point = None
        end_point = None
        print('Completed') 
def display_rotated_data(rotated_data, fig,ax,cmap='hot'):
    ax.clear()
    ax.imshow(rotated_data, cmap=cmap,origin='lower',vmax=rotated_data.max()/2)
    ax.set_axis_off()
    fig.canvas.draw()

resp_unit=0
exposure=0
RawDetectorDim = [1040,2152]
NonActivePix = [50, 8,4]
OutDir=None
hdul_all=[]
on_motion_cursor_2magixs_sigma = None
IncludePoisonNoise = False
StoreRawDetector = False
magixs_sigma = 0
# Mouse click event callback function
def on_click(event,img_data,EM_map,nFA,mag_resp,cmap, ref_map,fig,delta_pix_no, vig_func_data,vig_nyFA,mid_vig_fa_ind,mid_vig_fa_ind_y,detector_noise, ph2el,resp_unit=resp_unit,exposure=exposure,magixs_sigma=magixs_sigma,IncludePoisonNoise=IncludePoisonNoise,StoreRawDetector=StoreRawDetector,RawDetectorDim = RawDetectorDim,NonActivePix = NonActivePix,OutDir=OutDir,hdul_all=hdul_all,NoFrames=1,pix2wvl=None):
    global start_point, end_point, line, FOV
    axes = fig.get_axes()
    ax = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]
    #if event.inaxes is not None:
    if event.inaxes == ax:
        if event.button == 1:  # Left mouse button
            if start_point is None:
                start_point = (event.xdata, event.ydata)
                print(f"Start point: {start_point}")
            else:
                end_point = (event.xdata, event.ydata)
                print(f"End point: {end_point}")
                finalize_line(img_data,EM_map,nFA,mag_resp,cmap, ref_map,fig,delta_pix_no, vig_func_data,vig_nyFA,mid_vig_fa_ind,mid_vig_fa_ind_y,detector_noise, ph2el,resp_unit=resp_unit,exposure=exposure,magixs_sigma=magixs_sigma,IncludePoisonNoise=IncludePoisonNoise,StoreRawDetector=StoreRawDetector,RawDetectorDim = RawDetectorDim,NonActivePix = NonActivePix,OutDir=OutDir,hdul_all=hdul_all,NoFrames=NoFrames,pix2wvl=pix2wvl)
                FOV += 1


# Mouse motion event callback function
def on_motion_cursor_2(event,IMG,fig, wvl,pix_num,line_label,line_wvl,resp_unit,OutDir=None,namestr='',):
    global line,cursor_position
    def pix2wvl(x):
        return np.interp(x,pix_num,wvl)
    def wvl2pix(x):
        return np.interp(x,wvl,pix_num)

    axes = fig.get_axes()
    ax3 = axes[2]
    ax4 = axes[3]
    #if event.inaxes is not None:
    if event.inaxes == ax3:
        x = event.xdata ; y = event.ydata
        cursor_position = (x, y)
        line = plt.Line2D((0, len(IMG[0,:])),(int(y), int(y)),color='white')
        ax3.add_line(line)
        fig.canvas.draw()
        #print(f"Cursor position: {cursor_position}")
        ax4.clear()
        #ax4_second = ax4.secondary_xaxis('top',function(pix2wvl,wvl2pix))
        #spec = np.average(IMG[int(y)-1:int(y)+2,:],axis=0)
        spec = IMG[int(y),:]
        #specp = ax4.plot(np.arange(0,len(IMG[0,:])),spec,'b')
        specp = ax4.plot(wvl,spec,'b')
        if resp_unit == 'electron':ax4.set_ylabel('el s$^{-1}$ pix$^{-1}$')
        if resp_unit == 'photon': ax4.set_ylabel('Ph s$^{-1}$ pix$^{-1}$')
        ax4_second = ax4.secondary_xaxis('top',functions=(wvl2pix,pix2wvl))
        ax4_second.set_xlabel('Pixcel number')
        ax4.set_xlabel('Wavelength (${\AA}$)')
        ax4.set_xlim([wvl[0],wvl[-1]])
        #ax4.set_ylim([0,ax4.ylim()[1]])
        #primary_ticks = wvl2pix(ax4.get_xticks()[0:-1])
        #ax4_second.set_xticks(primary_ticks)
        #leb = [f'{wvl2pix(tick):.1f}' for tick in ax4.get_xticks()[0:-1]]
        #ax4_second.set_xticklabels(leb)
        ylim = ax4.get_ylim()
        for w in range(len(line_wvl)):
            x_ = line_wvl[w]
            ax4.text(x_,ylim[1]*0.78,line_label[w],fontsize=8,rotation=90,verticalalignment = 'bottom',horizontalalignment = 'center')
            ax4.plot([x_,x_],[0,ylim[1]*0.78],'--r',lw=1,alpha=0.6)
        #extent = ax4.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if OutDir is not None: fig.savefig(os.path.join(OutDir,namestr+'.png'))#,bbox_inches=extent)
        fig.canvas.draw()
        if line is not None:
            line.remove()

def GUI(FullSun_img_data,cmap,ref_map):
    fig = plt.figure(figsize=(14,8))
    y_loc = 0.42
    ax = fig.add_axes([0.0003,y_loc,0.21,0.7])
    ax2 = fig.add_axes([0.215,y_loc,0.185,0.7])
    ax3 = fig.add_axes([0.405,y_loc,0.59,0.7])
    ax4 = fig.add_axes([0.405,0.08,0.59,0.38])

    fig.text(0.01,0.2,'MaGIXS-2 Forward model\n\n  1. Select FOV center-line (left clicks)\n\n  2. Point cursor on overlapogram.\n\n  3. Close the widow to exit.',fontsize=12)
    ax.text(0.1,1.01,'Full Sun ('+ref_map.date.value[0:-7]+')',fontsize=12,transform = ax.transAxes)
    fig.text(0.26,0.96,'MaGIXS-2 FOV',fontsize=12,color='k')
    fig.text(0.65,0.96,'MaGIXS-2 Detector',fontsize=12,color='k')

    ax.imshow(FullSun_img_data, cmap=cmap,origin='lower',vmin=0,vmax=FullSun_img_data.max()/10)
    fig.text(0.01,0.2,'MaGIXS-2 Forward model\n\n  1. Select FOV center-line (left clicks)\n\n  2. Point cursor on overlapogram.\n\n  3. Close the widow to exit.',fontsize=12)
    ax.set_axis_off()
    ax2.set_axis_off()
    ax3.set_axis_off()
    ax4.set_axis_off()

    return fig,ax,ax2,ax3,ax4

def GUI_plot(FullSun_img_data,cmap,ref_map):
    fig = plt.figure(figsize=(14,8))
    y_loc = 0.42
    ax = fig.add_axes([0.0003,y_loc,0.21,0.7])
    ax2 = fig.add_axes([0.215,y_loc,0.185,0.7])
    ax3 = fig.add_axes([0.405,y_loc,0.59,0.7])
    ax4 = fig.add_axes([0.405,0.12,0.59,0.38])

    #fig.text(0.01,0.2,'MaGIXS-2 Forward model\n\n  1. Select FOV center-line (left clicks)\n\n  2. Point cursor on overlapogram.\n\n  3. Close the widow to exit.',fontsize=12)
    ax.text(0.1,1.01,'Full Sun ('+ref_map.date.value[0:-7]+')',fontsize=12,transform = ax.transAxes)
    fig.text(0.26,0.96,'MaGIXS-2 FOV',fontsize=12,color='k')
    fig.text(0.65,0.96,'MaGIXS-2 Detector',fontsize=12,color='k')

    ax.imshow(FullSun_img_data, cmap=cmap,origin='lower',vmin=0,vmax=FullSun_img_data.max()/10)
    #fig.text(0.01,0.2,'MaGIXS-2 Forward model\n\n  1. Select FOV center-line (left clicks)\n\n  2. Point cursor on overlapogram.\n\n  3. Close the widow to exit.',fontsize=12)

    x_w = 0.12
    y_loc = 0.12
    ax5 = fig.add_axes([0.001,y_loc,x_w,0.4])
    ax6 = fig.add_axes([0.001+x_w+0.002,y_loc,x_w,0.4])
    ax7 = fig.add_axes([0.001+2*(x_w+0.002),y_loc,x_w,0.4])

    ax.set_axis_off()
    ax2.set_axis_off()
    ax3.set_axis_off()
    ax4.set_axis_off()
    ax5.set_axis_off()
    ax6.set_axis_off()
    ax7.set_axis_off()

    return fig,ax,ax2,ax3,ax4,ax5,ax6,ax7

def search_and_replace_line(file_path, search_text_start, replace_text):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    found = False
    for i, line in enumerate(lines):
        if line.startswith(search_text_start):
            lines[i] = replace_text + '\n'  # Add a newline character to maintain file formatting
            found = True
            break
    
    if not found:
        print("Search text not found in the file.")
    
    with open(file_path, 'w') as file:
        file.writelines(lines)

def add_colorbar(ax, array, cmap, label):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    norm = plt.Normalize(vmin=array.min(), vmax=array.max())
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,orientation='horizontal')
    cbar.set_label(label, rotation=0, labelpad=15)

def get_line_start_end(roll, centerX,centerY,length = 50):
    angle_rad = np.deg2rad(roll)
    x_offset = length / 2 * np.cos(angle_rad)
    y_offset = length / 2 * np.sin(angle_rad)
    start_point = (int(centerX - x_offset), int(centerY - y_offset))
    end_point = (int(centerX + x_offset), int(centerY + y_offset))
    return start_point,end_point
