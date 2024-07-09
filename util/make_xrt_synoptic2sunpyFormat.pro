pro make_xrt_synoptic2sunpyFormat

  ;------ Inputs ---------------
  data_dir = '/Volumes/BM_WD/BKG/NPP/MaGIXS-2/simulated_spec/20240619T1200/XRT'
  
  Out_dir = data_dir+'/prep/'
  ;-----------------------------

  ff18=file_search(data_dir+'/*.fits')
  ;ff18=[ff18[13]]
  for i=0,n_elements(ff18)-1 do begin
      name = (strsplit(ff18[i],'/',/extract))[-1]
      read_xrt,ff18[i],ind,data,/force
      ; grade_type=1 as these are summed data of 512x512 (so each pixel is sum of 4x4)
      xrt_prep,ind,data,indp,datap,/float,grade_map=gm,grade_type=1,/coalign
      ;xrt_jitter, indp, offsets
      ;data = image_translate(data, offsets)
      ; Change so compatiable with sunpy better
      indp.timesys='UTC'
      filtout=indp.ec_fw1_+'-'+indp.ec_fw2_
      resout=strcompress(string(indp.naxis1),/rem)
      fnameout='XRTcomp'+'_'+filtout+'.fits'
      write_xrt,indp,datap,outdir=Out_dir,outfile=fnameout,/ver
  endfor
  stop
end
