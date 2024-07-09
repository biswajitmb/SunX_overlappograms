'''
Purpose: It will use the outputs of "Get_DEM_map_from_ObsImage.py" or "Get_DEM_map_from_PredImage.py"


Note: Implementing parallel processing requires 'spawn' method, with 'if __name__ == '__main__':'
to avoid running whole scripts in each child process. This requires to pass all the variables in parallel function
'''


# Import some of the stuff we will need
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
import scipy.io as io
import glob
from aiapy.calibrate import estimate_error
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.colors as colors

# So can have single copy of demreg on system, and don't need copy in working directory
from sys import path as sys_path
# Change to your local copy's location...
sys_path.append('/Users/bmondal/BM_Works/softwares/biswajit_python_lib')

import astropy.time as atime
from astropy.coordinates import SkyCoord
from astropy import units as u
import sunpy.map

from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table
from aiapy.calibrate import register, update_pointing

import warnings
warnings.simplefilter('ignore')
matplotlib.rcParams['font.size'] = 16
from biswajit_utils_functions import*
#from scipy.odr import ODR, Model, Data, RealData
from sklearn.linear_model import LinearRegression
from sklearn.impute import SimpleImputer

start_time = time.time()

#-------- inputs ---------
DEM_dataDir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/DEM_data/20240323T0330/DEMs/'
dem_files = glob.glob(DEM_dataDir+'*.pkl')

Store_outputs = False

OutDir = DEM_dataDir+'slopes/'

FitDEM = False #if True then the DEM will be fitted by two power-law (clod and hot) corresponding to below and above peak-T.
FitDEM_MinMaxlogT = [5.8,6.7] #logT limits for [min,max] for the cold and hot slope

#store_outputs = True

plot_EM_locii = True

#-------------------------

xr = [5.6,7.0]
yr= [1.0e16,1.0e24]

# based off of my previous get_error_area.pro
# Mostly works, but interp can go a bit funny at times
def get_error_bound(d,ed,l,el,xlim=[5.7,7.2],ylim=[2e20,4e23]):
    n=len(l)

    # interpolate horizontal error limits the temperature binning
    d0=d
    d0[d0<= ylim[0]]=ylim[0]

    m1=np.interp(l,l-el,d0)
    m2=np.interp(l,l+el,d0)
    # m1=np.interp(l-el,l,d0)
    # m2=np.interp(l+el,l,d0)

    mind=d-ed
    # some checks to make sure within chosen limts, if not set to them
    mind[mind<=ylim[0]]=ylim[0]
    mind[mind>=ylim[1]]=ylim[1]

    bot=np.zeros(n)
    for i in np.arange(0,n): 
        bot[i]=np.min([m1[i],m2[i],mind[i]])

    # do the top most points
    maxd=d+ed
    maxd[maxd<=ylim[0]]=ylim[0]
    maxd[maxd>=ylim[1]]=ylim[1]

    top=np.zeros(n)
    for i in np.arange(0,n): 
        top[i]=np.max([m1[i],m2[i],maxd[i]])

    # create the arrays ready for the polyfill plotting
    xx=np.concatenate(([l[0]-el[0]],[l[0]],l,[l[-1]+el[-1]],[l[-1]+el[-1]],np.flip(l),[l[0]],[l[0]-el[0]]))
    yy=np.concatenate(([top[0]],[top[0]],top,[top[-1]],[bot[-1]],np.flip(bot),[l[0]],[l[0]]))

    yy[yy < ylim[0]]=ylim[0]
    yy[yy > ylim[1]]=ylim[1]
    xx[xx < xlim[0]]=xlim[0]
    xx[xx > xlim[1]]=xlim[1]

    return xx,yy


def DEM2EM(log10T,DEM):
    dT = (shift(log10T, -1, cval=0.0) - shift(log10T, 1, cval=0.0)) * 0.5
    ntemps = len(log10T)
    dT[0] = log10T[1] - log10T[0]
    dT[ntemps-1] = (log10T[ntemps-1] - log10T[ntemps-2])
    EM = DEM * (10**log10T) *np.log(10.) * dT
    return EM

def EM2DEM(log10T,EM):
    dT = (shift(log10T, -1, cval=0.0) - shift(log10T, 1, cval=0.0)) * 0.5
    ntemps = len(log10T)
    dT[0] = log10T[1] - log10T[0]
    dT[ntemps-1] = (log10T[ntemps-1] - log10T[ntemps-2])
    DEM = EM / ((10**log10T)*np.log(10.) * dT)
    return DEM

'''
def EM_polow_cool(par,log10T):
    return (par[0] + par[1]*log10T) #return log10(EM)
'''
def EM_polow_cool(log10T,alpha,norm):
    return (norm + alpha*log10T) #return log10(EM)

def EM_polow_hot(log10T,beta,norm):
    return (norm + beta*log10T) #return log10(EM)

def EM2DEM(log10T,EM):
    dT = (shift(log10T, -1, cval=0.0) - shift(log10T, 1, cval=0.0)) * 0.5
    ntemps = len(log10T)
    dT[0] = log10T[1] - log10T[0]
    dT[ntemps-1] = (log10T[ntemps-1] - log10T[ntemps-2])
    DEM = EM / ((10**log10T)*np.log(10.) * dT)
    return DEM


def mouse_event(event):
    global fig2,ax2,ax22,DEM_map,DEM_map_err,DEM_logT_mid,DEM_map_logT_err,DEM_map_chisq,tresp_logt,tresp,Obs_Img, xr,yr
    if event.button == 3: 
        ix,iy  = event.xdata,event.ydata
        #mouse_event.ix = ix #to access the value outside the function
        #mouse_event.iy = iy
        
        pix_loc = [int(ix),int(iy)]
        ax2.cla() #clear the previous datapoints
        ax2.set_yscale('log')
        ax2.set_ylabel('DEM (cm$^{-5}$ K$^{-1}$)')
        ax2.set_xlabel('log$_{10}$(T)')
        ax2.set_xlim(xr)
        ax2.set_ylim(yr)

        #Plot EM locii
        plot_EM_locii = True
        if plot_EM_locii is True:
            for e in range(len(Obs_Img[0,0,:])):
                ax2.plot(tresp_logt,EM2DEM(tresp_logt,(Obs_Img[int(ix),int(iy),e]/tresp[:,e])))

        ax2.errorbar(DEM_logT_mid,DEM_map[pix_loc[0],pix_loc[1],:],xerr=DEM_map_logT_err[pix_loc[0],pix_loc[1],:],yerr=DEM_map_err[pix_loc[0],pix_loc[1],:],fmt='.r',\
             ecolor='lightcoral', elinewidth=3, capsize=0)
        ax2.step(DEM_logT_mid,DEM_map[pix_loc[0],pix_loc[1],:],where='mid',color='red',label='$\Delta\chi$ = '+format('%0.2f'%DEM_map_chisq[pix_loc[0],pix_loc[1]]))
        xx0,yy0=get_error_bound(DEM_map[pix_loc[0],pix_loc[1],:],DEM_map_err[pix_loc[0],pix_loc[1],:],DEM_logT_mid,DEM_map_logT_err[pix_loc[0],pix_loc[1],:],xlim=xr,ylim = yr)
        ax2.fill(xx0,yy0,"lightcoral",alpha=0.5)
        #ax2.legend(loc='upper right')
        
        #Fit the cool and hot EM curve
        DEM_logT_mid = np.array(DEM_logT_mid)
        x = DEM_logT_mid#-(0.05) #temperature edge
        y = np.log10(DEM2EM(x,DEM_map[pix_loc[0],pix_loc[1],:]))
        yerr = np.array(np.log10(DEM2EM(x,DEM_map[pix_loc[0],pix_loc[1],:]+DEM_map_err[pix_loc[0],pix_loc[1],:])) - y)
        ind = np.where(((x >= FitDEM_MinMaxlogT[0])&(x <= FitDEM_MinMaxlogT[1])))
        x=x[ind]
        y = y[ind]
        yerr = yerr[ind]
        xerr = DEM_map_logT_err[ind]
        max_ind = np.where(y==max(y))[0][0]
        max_ind_a = max_ind+1
        #par,err_cov = curve_fit(EM_polow_cool, x[0:int(max_ind)], y[0:int(max_ind)],p0=[1,1.0],bounds=([-np.inf,0], [np.inf,np.inf]),sigma=yerr[0:int(max_ind)],absolute_sigma=True)

        #Using linear corellation
        model = LinearRegression()
        model.fit(x[0:int(max_ind_a)].reshape((-1, 1)), y[0:int(max_ind_a)])
        r_sq = model.score(x[0:int(max_ind_a)].reshape((-1, 1)), y[0:int(max_ind_a)])
        print(f"coefficient of determination: {r_sq}")
        par = [model.coef_[0],model.intercept_]
        print('alpha=',par)
        y_model = EM2DEM(x[0:int(max_ind_a)],10**EM_polow_cool(x[0:int(max_ind_a)],*par))
        ax2.plot(x[0:int(max_ind_a)],y_model,color='b')
       
        #Fit the hot EM curve   
            
        model.fit(x[int(max_ind)::].reshape((-1, 1)), y[int(max_ind)::])
        r_sq = model.score(x[int(max_ind)::].reshape((-1, 1)), y[int(max_ind)::])
        par = [model.coef_[0],model.intercept_]
        print('beta',par)
        y_model = EM2DEM(x[int(max_ind)::],10**EM_polow_hot(x[int(max_ind)::],*par))
        ax2.plot(x[int(max_ind)::],y_model,color='b')
 
        fig2.canvas.draw() #redraw the figure       
         

def DemMap2EffT(DEM_logT,DEM_Map,Non_zero_pixcels):
    '''
    inputs:
        DEM_Map -> 3D array, dimension = [x,y,logT]
        DEM_logT -> DEM logT grids
        Non_zero_pixcels -> 2-column array, representing the indices of non-zero pixcels to which the counts is to be 
                          calculated, to reduce the computation time.
    outputs: Eff_Tmap --> EM-waighted temperature
    '''
    DEM_logT = DEM_logT[0::]
    DEM_Map = DEM_Map[:,:,0::]
    dT = (shift(DEM_logT, -1, cval=0.0) - shift(DEM_logT, 1, cval=0.0)) * 0.5
    ntemps = len(DEM_logT)
    dT[0] = DEM_logT[1] - DEM_logT[0]
    dT[ntemps-1] = (DEM_logT[ntemps-1]-DEM_logT[ntemps-2])

    #indices of non-zero pixcels
    ind_xx = Non_zero_pixcels[0][:]
    ind_yy = Non_zero_pixcels[1][:]

    Model_EM = DEM_Map[:,:,:]*0
    for i in range(len(ind_xx)):
        for j in range(len(DEM_logT)):
            dem = DEM_Map[ind_xx[i],ind_yy[i],j]
            Model_EM[ind_xx[i],ind_yy[i],j] += (dem * (10**DEM_logT[j]) *np.log(10.) * dT[j])

    Eff_Tmap = DEM_Map[:,:,0]*0
    for i in range(len(ind_xx)):
        EM_i = Model_EM[ind_xx[i],ind_yy[i],:]
        Eff_Tmap[ind_xx[i],ind_yy[i]] = sum(EM_i*(10**DEM_logT)) / sum(EM_i)

    return Eff_Tmap

for img in range(len(dem_files)):
    #load DEM_map data
    data = load_obj(dem_files[img][0:-4])
    ind = 0
    DEM_map = data['ObsDEMmap']    
    DEM_map_err = data['ObsDEMmap_Err']
    DEM_logT_mid = data['ObsDEMmap_logt_midbins']
    DEM_map_logT_err = data['ObsDEMmap_logtErr']
    DEM_map_chisq = data['chisq']
    Obs_Img = data['ObsImg'] #[x,y,channels]
    tresp = data['Tresp'] #[temperaure, channels]
    tresp_logt = data['Tresp_logt']

    #Get the location of non-zero pixcels
    dem_tot = np.sum(DEM_map,axis = 2)
    NonZeroPx = np.where(dem_tot > 0)

    DEM_logT_edge = np.array(DEM_logT_mid)-((DEM_logT_mid[2]-DEM_logT_mid[1])/2)
    Eff_Tmap = DemMap2EffT(DEM_logT_edge,DEM_map,NonZeroPx)

    
    fig1 = plt.figure(figsize=(4,4))
    fig2 = plt.figure(figsize=(6,6))
    wd_a = 0
    fnt_size = 10
    lev_size = 10
    ax1 = fig1.add_axes([0.1,0.1,0.8,0.8])
    ax2 = fig2.add_axes([0.15,0.15,0.8,0.8])
    #ax22 = fig2.add_axes([0.1,0.1,0.8,0.2])
    
    #ax1.imshow(data0.data,cmap=data0.cmap,origin='lower',interpolation='nearest',norm=colors.PowerNorm(gamma=0.5))
    pred_T = ax1.imshow(Eff_Tmap/1.0e6,origin='lower',interpolation='gaussian',cmap="gist_ncar",alpha=0.8)#hot_r")
    plt.colorbar(pred_T, ax=ax1)
    ax1.set_title('EM weighted T \n(right-click on pixcel)')
    cnv = fig1.canvas.mpl_connect('button_press_event', mouse_event)
    ax1.set_axis_off()
    plt.show()
  
    ''' 
    #indices of non-zero pixcels
    ind_xx = NonZeroPx[0][:]
    ind_yy = NonZeroPx[1][:]
    EM_slopeAlpha_map = Eff_Tmap[:,:]*0
    EM_slopeBeta_map = Eff_Tmap[:,:]*0
    #if store_outputs is True:
    #    #fa = open(dem_files[img][0:-4]+'_alpha.dat','w')
    #    #fb = open(dem_files[img][0:-4]+'_beta.dat','w')
    alpha_all = []
    beta_all = []
    for m in range(len(ind_xx)):
        pix_loc = [ind_xx[m],ind_yy[m]]
 
        #Fit the cool and hot EM curve
        DEM_logT_mid = np.array(DEM_logT_mid)
        x = DEM_logT_mid#-(0.05) #temperature edge
        y = np.log10(DEM2EM(x,DEM_map[pix_loc[0],pix_loc[1],:]))
        yerr = np.array(np.log10(DEM2EM(x,DEM_map[pix_loc[0],pix_loc[1],:]+DEM_map_err[pix_loc[0],pix_loc[1],:])) - y)
        ind = np.where(((x >= FitDEM_MinMaxlogT[0])&(x <= FitDEM_MinMaxlogT[1])))
        x=x[ind]
        y = y[ind]
        yerr = yerr[ind]
        xerr = DEM_map_logT_err[ind]
        max_ind = np.where(y==max(y))[0]
        if max_ind > 0: 
            max_ind = max_ind
            max_ind_a = max_ind+1
            #par,err_cov = curve_fit(EM_polow_cool, x[0:int(max_ind)], y[0:int(max_ind)],p0=[1,1.0],bounds=([-np.inf,0], [np.inf,np.inf]),sigma=yerr[0:int(max_ind)],absolute_sigma=True)

            #Using linear corellation
            try:
                #imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
                #y[0:int(max_ind)] = imputer.fit_transform(y[0:int(max_ind)].reshape(-1, 1))[:,0] #fill in the null values in y with nan
                model = LinearRegression()
                model.fit(x[0:int(max_ind_a)].reshape((-1, 1)), y[0:int(max_ind_a)])
                r_sq = model.score(x[0:int(max_ind_a)].reshape((-1, 1)), y[0:int(max_ind_a)])
                par = [model.coef_[0],model.intercept_]
                     
                if r_sq> 0.9 and r_sq < 1.0 : 
                    print('alpha',par)
                    #if store_outputs is True: fa.write('%0.2f\t'%par[0]) ; fa.write('%0.2f\n'%par[1])
                    alpha_all += [par[0]]
                    EM_slopeAlpha_map[pix_loc[0],pix_loc[1]] = par[0]
            except: pass

            #Fit the hot EM curve   
            try:
                model = LinearRegression()
                model.fit(x[int(max_ind)::].reshape((-1, 1)), y[int(max_ind)::])
                r_sq = model.score(x[int(max_ind)::].reshape((-1, 1)), y[int(max_ind)::])
                par = [model.coef_[0],model.intercept_]
                if r_sq> 0.9 and r_sq < 1.0 :
                    print('beta',par)
                    #if store_outputs is True: fb.write('%0.2f\t'%par[0]) ; fb.write('%0.2f\n'%par[1])
                    beta_all += [par[0]]
                    EM_slopeBeta_map[pix_loc[0],pix_loc[1]] = par[0]
            except: pass
    #fa.close()
    #fb.close()    
    plt.hist(alpha_all)
    plt.hist(beta_all)
    plt.show()
    if store_outputs is True: 
        results = {'DEM Slopes Map': 'code: HK_dem'}
        results['EM_weighted_Tmap'] = Eff_Tmap 
        results['EM_slopeAlpha_map'] = EM_slopeAlpha_map
        results['EM_slopeBeta_map'] = EM_slopeBeta_map
        name = dem_files[img][0:-4].split('/')[-1]
        save_obj(results, OutDir+'Slopes_'+name)
    '''
print("END Calculation--- %s seconds ---" % (time.time() - start_time),'\n')
