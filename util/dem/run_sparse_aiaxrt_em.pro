
PRO run_sparse_aiaxrt_em

;===================

Data_dir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/'

aia_files = ['AIA_rebin_level1.5_wvn094.fits','AIA_rebin_level1.5_wvn131.fits','AIA_rebin_level1.5_wvn171.fits','AIA_rebin_level1.5_wvn193.fits','AIA_rebin_level1.5_wvn211.fits','AIA_rebin_level1.5_wvn335.fits']

;aia_files = ['../AIA_level1.5_wvn094.fits','../AIA_level1.5_wvn131.fits','../AIA_level1.5_wvn171.fits','../AIA_level1.5_wvn193.fits','../AIA_level1.5_wvn211.fits','../AIA_level1.5_wvn335.fits']

xrt_files = ['XRT_rebin_XRTcomp_Be_thin-Open.fits']

OutDir = Data_dir+'DEMs/'

plate_scale = 2.8
crosscal_factor = 2
;===================

fits2map,Data_dir+'AIA/rebined/'+aia_files[0],a94map
fits2map,Data_dir+'AIA/rebined/'+aia_files[1],a131map
fits2map,Data_dir+'AIA/rebined/'+aia_files[2],a171map
fits2map,Data_dir+'AIA/rebined/'+aia_files[3],a193map
fits2map,Data_dir+'AIA/rebined/'+aia_files[4],a211map
fits2map,Data_dir+'AIA/rebined/'+aia_files[5],a335map

if n_elements(xrt_files) gt 0 then begin
    xrtf = Data_dir+'XRT/rebined/'+xrt_files
    read_xrt, xrtf, index, data
    index2map,index, data,xrtBeThinMap
    bethin_mag_rbin = xrtBeThinMap.data
    xrt_exptime = xrtBeThinMap.dur
endif

aia_exptime=[a94map.dur,a131map.dur,a171map.dur,a193map.dur,a211map.dur,a335map.dur]

;stop

aia_DX = a94map.DX
aia_DY = a94map.DY
a94_mag_rbin  = a94map.data
a131_mag_rbin = a131map.data 
a171_mag_rbin = a171map.data
a193_mag_rbin = a193map.data
a211_mag_rbin = a211map.data
a335_mag_rbin = a335map.data

nx_magpix=(size(a94_mag_rbin))[1]
ny_magpix=(size(a94_mag_rbin))[2]

n_chn = n_elements(aia_files) + n_elements(xrt_files)

aia_mag_cube = dblarr(nx_magpix,ny_magpix,6)
aia_mag_cube[*,*,0]=a94_mag_rbin
aia_mag_cube[*,*,1]=a131_mag_rbin
aia_mag_cube[*,*,2]=a171_mag_rbin
aia_mag_cube[*,*,3]=a193_mag_rbin
aia_mag_cube[*,*,4]=a211_mag_rbin
aia_mag_cube[*,*,5]=a335_mag_rbin

img_cube = dblarr(nx_magpix,ny_magpix,n_chn)

img_cube[*,*,0:5]=aia_mag_cube
if n_elements(xrt_files) gt 0 then img_cube[*,*,6]=bethin_mag_rbin
;do sparse solve only for aia
lgTmin = 5.6
dlgT = 0.1
nlgT = 15

date = strmid(a94map.TIME,7,4)+'-'+strmid(a94map.TIME,3,3)+'-'+strmid(a94map.TIME,0,2)

;===do sparse solve for AIA+XRT
chn_all=0
tresp_all=0
use_lgtaxis = findgen(nlgT)*dlgT+lgTmin
aia_sparse_em_init_mod, timedepend=date,/evenorm,use_lgtaxis=use_lgtaxis,xrt=['Be-thin'],plate_scale = plate_scale,crosscal_factor=crosscal_factor,chn_all=chn_all,tresp_all=tresp_all
exptime=aia_exptime[*]
exptimestr = '['+strjoin(string(aia_exptime[*],format='(F8.5)'),',')+']'
if n_elements(xrt_files) gt 0 then begin
    xrtexptimestr=[string(xrt_exptime,format='(F8.5)')] 
    tolfunc =  '[aia_bp_estimate_error(y[0:5]*'+exptimestr+$
    ',[94,131,171,193,211,335],num_images='+strtrim(string(1^2,format='(I4)'),2)+')/'+exptimestr+$
    ', sqrt((0.25*y[6])^2' + '+ (sqrt(y[6]*' + xrtexptimestr[0]+')/'+xrtexptimestr[0]+')^2)'+']'
endif else begin
    tolfunc =  '[aia_bp_estimate_error(y[0:5]*'+exptimestr+$
    ',[94,131,171,193,211,335],num_images='+strtrim(string(1^2,format='(I4)'),2)+')/'+exptimestr+']'
endelse

;joint DEM solver aia-xrt
aia_sparse_em_solve,img_cube, tolfunc=tolfunc, tolfac=1.4, oem=emcube, $
    status=status, coeff=coeff, adaptive_tolfac=adaptive_tolfac,tolmap=tolmap

;filter=['A94 A','131 A','171 A', '193 A', '211 A', '335 A', 'Be-thin']

;xrt_dem_iterative2,filter,dn_in,trmatrix_struct_new,$
    ;logt_out,dem_out_fox,obs_err=edn_in,mc_iter=100,$
    ;min_t=mint_range,max_t=maxt_range,dt=0.05,base_obs=base_obs_fox, $
    ;mod_obs=mod_obs_fox,chisq=chisq_fox,qabort=qab_fox,/quiet




;aia_sparse_dem_inspect, coeff, emcube, status, /ylog,img=aia_mag_cube

;solve only for aia
;aia_sparse_em_init, timedepend = date, /evenorm,use_lgtaxis=findgen(nlgT)*dlgT+lgTmin
;tolfunc = 'aia_bp_estimate_error(y*'+exptimestr+$
	;',[94,131,171,193,211,335],num_images='+strtrim(string(1^2,format='(I4)'),2)+')/'+exptimestr
;aia_sparse_em_solve,aia_mag_cube, tolfunc=tolfunc, tolfac=1.4, $
	;oem=emcube_aia, status=status_aia,coeff=coeff_aia,adaptive_tolfac=adaptive_tolfac_aia,$
	;tolmap=tolmap_aia

  ;whereis,status_aia,where(status_aia ne 0.),xx_aia,yy_aia
    ;b4_ip_emcube_aia=emcube_aia
    ;for i =0,n_elements(yy_aia)-1 do $
	;if ((xx_aia[i]-1 gt 1) and (yy_aia[i]-1 gt 1) and $
	;(xx_aia[i]+1 lt n_elements(emcube_aia[*,0,0])-1) and $
	;(yy_aia[i]+1 lt n_elements(emcube_aia[*,0,0])-1)) then $
	;for j=0,n_elements(emcube_aia[0,0,*])-1 do $
	;emcube_aia[xx_aia[i],yy_aia[i],j]=mean(b4_ip_emcube_aia[xx_aia[i]-2:xx_aia[i]+2,yy_aia[i]-2:yy_aia[i]+2,j])




st_nt=strtrim(string(nlgt),2)+'nlgT'

whereis,status,where(status ne 0.),xx,yy
b4_ip_emcube=emcube
;for i =0,n_elements(yy)-1 do $
    ;emcube[xx[i],yy[i],*]=emcube_aia[xx[i],yy[i],*]

emcube=b4_ip_emcube
p_smear=10
for i =0,n_elements(yy)-1 do $
    if ((xx[i]-p_smear gt 1) and (yy[i]-p_smear gt 1) and $
    (xx[i]+p_smear lt n_elements(emcube[*,0,0])-1) and $
    (yy[i]+p_smear lt n_elements(emcube[*,0,0])-1)) then $
    for j=0,n_elements(emcube[0,0,*])-1 do $
    emcube[xx[i],yy[i],j]=mean(b4_ip_emcube[xx[i]-p_smear:xx[i]+p_smear,yy[i]-p_smear:yy[i]+p_smear,j])
imgcube=img_cube
;save,imgcube,filename=OutDir+'/xrtshift_aiaxrt_imgcubes_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
;save,b4_ip_emcube,filename=OutDir+'/xrtshift_aiaxrt_b4_ip_emcubes_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
;save,emcube,filename=OutDir+'/xrtshift_aiaxrt_emcubes_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
;save,status,filename=OutDir+'/xrtshift_aiaxrt_status_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
;save,coeff,filename=OutDir+'/xrtshift_aiaxrt_coeff_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'

save,emcube,status,Data_dir,aia_files,xrt_files,OutDir,plate_scale,crosscal_factor,use_lgtaxis,chn_all,tresp_all,filename=OutDir+'/ObsDEM_AIA_XRT_dem.sav'

loadct,3
aia_sparse_dem_inspect, coeff, emcube, status, /ylog,img=imgcube

;a,b,c = get_required_outputs

stop
end


