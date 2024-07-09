; IDL procedure by Athiray
; Copyright (c) 2019, MaGIXS Sounding rocket mission, NASA Marshall Space Flight Center,  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 28 Sep 2020 21:48
; Last Mod 	: 20 Oct 2021 21:36

;-------------------  Details of the program --------------------------;



PRO make_truth_aiaxrt_em,ipdir,xrtl0l1=xrtl0l1

if ((ipdir eq '2021-07-27') or (ipdir eq '2021-07-28') or (ipdir eq '2021-07-29')) then $
;if ((ipdir eq '2021-07-27') or (ipdir eq '2021-07-28') or (ipdir eq '2021-07-29') or (ipdir eq '2021-07-30')) then $
    files=file_search(ipdir+'/aia*.fts') else $
    files=file_search(ipdir+'/aia*.fits')

files=[files[5],files[0:4]]

fits2map,files[0],a94map
fits2map,files[1],a131map
fits2map,files[2],a171map
fits2map,files[3],a193map
fits2map,files[4],a211map
fits2map,files[5],a335map


;mreadfits,files[0],a94in,/nodata
;mreadfits,files[1],a131in,/nodata
;mreadfits,files[2],a171in,/nodata
;mreadfits,files[3],a193in,/nodata
;mreadfits,files[4],a211in,/nodata
;mreadfits,files[5],a335in,/nodata


xrtf=file_search('/Users/spanchap/python_exercise/magixs_python/XRT_sunpy/'+ipdir+'/*.fits')
xrtf=xrtf[1]
read_xrt, xrtf, index, data
if keyword_set(xrtl0l1) then begin
    xrt_prep, index, data, index_out, data_out,despike_despot, /float
    if ipdir eq '2021-07-30' then begin
	data_out=shift(data_out,12,12)
	data_out[0:11,0:11]=0
    endif
    xrtmap=congrid(data_out/4.,4096,4096)
    xrt_exptime=index_out.exptime
endif else begin
    fits2map,xrtf,xrtm
    xrtmap=congrid(xrtm.data/4.,4096,4096)
    xrt_exptime=xrtm.dur
endelse

;xmap=make_map(xrtmap,dx=2.06/4.,dy=2.06/4.,xc=xrtm.xc,yc=xrtm.yc)
aia_exptime=[a94map.dur,a131map.dur,a171map.dur,a193map.dur,a211map.dur,a335map.dur]


stop

n_aia_pix=4096
a94_mag_rbin=congrid(a94map.data,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
a131_mag_rbin=congrid(a131map.data,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
a171_mag_rbin=congrid(a171map.data,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
a193_mag_rbin=congrid(a193map.data,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
a211_mag_rbin=congrid(a211map.data,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
a335_mag_rbin=congrid(a335map.data,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
;bethin_mag_rbin=congrid(xrtmap,n_aia_pix*(0.6/2.8),n_aia_pix*(0.6/2.8));,/interp)
be_mag_rbin=congrid(xrtmap,n_aia_pix*(2.06/(4*2.8)),n_aia_pix*(2.06/(4*2.8)));,/interp)
nbex=n_elements(be_mag_rbin[*,0])
nbey=n_elements(be_mag_rbin[0,*])

mag_si=size(a94_mag_rbin)
n_magpix=mag_si[1]

bethin_mag_rbin=dblarr(n_magpix,n_magpix)
bethin_mag_rbin[62:62+nbex-1,62:62+nbey-1]=be_mag_rbin



aia_mag_cube = dblarr(n_magpix,n_magpix,6)
aia_mag_cube[*,*,0]=a94_mag_rbin
aia_mag_cube[*,*,1]=a131_mag_rbin
aia_mag_cube[*,*,2]=a171_mag_rbin
aia_mag_cube[*,*,3]=a193_mag_rbin
aia_mag_cube[*,*,4]=a211_mag_rbin
aia_mag_cube[*,*,5]=a335_mag_rbin

img_cube = dblarr(n_magpix,n_magpix,7)

img_cube[*,*,0:5]=aia_mag_cube
img_cube[*,*,6]=bethin_mag_rbin
;do sparse solve only for aia
lgTmin = 5.6
dlgT = 0.1
nlgT = 15
;nlgT = 18
;dlgT = 0.05
;nlgT = 30
month=strmid(ipdir,5,2)
if month eq '01' then mn='Jan'
if month eq '02' then mn='Feb'
if month eq '03' then mn='Mar'
if month eq '04' then mn='Apr'
if month eq '05' then mn='May'
if month eq '06' then mn='Jun'
if month eq '07' then mn='Jul'
if month eq '08' then mn='Aug'
if month eq '09' then mn='Sep'
if month eq '10' then mn='Oct'

date='2021-'+mn+'-'+strmid(ipdir,8,2)
;===do sparse solve for AIA+XRT
aia_sparse_em_init, timedepend=date,/evenorm,$
    use_lgtaxis=findgen(nlgT)*dlgT+lgTmin,$
    xrt=['Be-thin'],/crosscalib
exptime=aia_exptime[*]
exptimestr = '['+strjoin(string(aia_exptime[*],format='(F8.5)'),',')+']'
xrtexptimestr=[string(xrt_exptime,format='(F8.5)')]
tolfunc =  '[aia_bp_estimate_error(y[0:5]*'+exptimestr+$
    ',[94,131,171,193,211,335],num_images='+strtrim(string(1^2,format='(I4)'),2)+')/'+exptimestr+$
    ', sqrt((0.25*y[6])^2' + '+ (sqrt(y[6]*' + xrtexptimestr[0]+')/'+xrtexptimestr[0]+')^2)'+']'
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
save,imgcube,filename=ipdir+'/xrtshift_aiaxrt_imgcubes_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
save,b4_ip_emcube,filename=ipdir+'/xrtshift_aiaxrt_b4_ip_emcubes_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
save,emcube,filename=ipdir+'/xrtshift_aiaxrt_emcubes_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
save,status,filename=ipdir+'/xrtshift_aiaxrt_status_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'
save,coeff,filename=ipdir+'/xrtshift_aiaxrt_coeff_stack_rbin2p8_'+ipdir+'_'+st_nt+'.sav'

loadct,3
aia_sparse_dem_inspect, coeff, emcube, status, /ylog,img=imgcube
stop
end


