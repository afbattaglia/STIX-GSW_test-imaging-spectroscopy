;+
;
; PROJECT:
;   STIX
;
; NAME:
;   stx_plot_imaging_spectra
;
; PURPOSE:
;   Plot the spectra obtained with stx_imaging_spectroscopy, returns a structure
;   containing the fluxes of the different sources and the same structure is 
;   saved in a sav file. The spectra and the sav files are saved in the
;   path_data_folder. This procedure is meant to be used with the save files 
;   created by stx_imaging_spectroscopy.
;
; CALLING SEQUENCE:
;   stx_plot_imaging_spectra, path_data_folder, path_sci_file, path_bkg_file
;
; INPUTS:
;   path_data_folder   : path to the folder containing the save files generated
;                        by stx_imaging_spectroscopy. This should correspond to
;                        path_new_folder that is returned by stx_imaging_spectroscopy
;
;   path_sci_file      : path to the science fits file.
;                        It has to be the same as the cpd used for imaging
;
;   path_bkg_file      : path to the background fits file.
;                        It has to be the same as the background used for imaging
;
; OPTIONAL INPUTS AND KEYWORDS:
;   color_srcs         : array of colors for the different sources
;                        default: ['red', 'blue', 'orange', 'magenta', 'dark green', 'cyan', 'saddle brown', 'goldenrod']
;
;   observed_color     : color for the observed (spatially integrated) spectrum
;                        default: 'black'
;
;   stop_here          : if set, the procedure stops before exiting (for debugging)
;
; HISTORY:
;   July 2023, Battaglia A. F. (FHNW & ETHZ), initial release
;
; CONTACT:
;   andrea.battaglia@fhnw.ch
;
;-

pro stx_plot_imaging_spectra, path_data_folder, path_sci_file, path_bkg_file, $    ; inputs
  ;; --- Optional inputs and keywords
  color_srcs = color_srcs, $      
  observed_color = observed_color, $
  stop_here = stop_here
  


  ;;;;; Cosmetics plot
  default, color_srcs, ['red', 'blue', 'orange', 'magenta', 'dark green', 'cyan', 'saddle brown', 'goldenrod']  
  default, observed_color, 'black'
  
  
  ;;;;; Get the folder delimiter (different for different OS)
  folder_delimiter = path_sep()


  ;;;;; Check if path_data_folder is effectively a folder
  if not path_data_folder.endswith(folder_delimiter) then path_data_folder = path_data_folder.insert(folder_delimiter,path_data_folder.strlen())

  
  ;;;;; PS filenames
  ps_flnm_spectra = path_data_folder + 'spectra-sources-STIX.ps'
  
  
  ;;;;; Find the total number of sources in the series
  all_path_sav_stix = findfile(path_data_folder+'stix-imaging*.sav')
  array_nsources = []
  n_ebins = n_elements(all_path_sav_stix)
  for this_f=0,n_ebins-1 do begin
    restore,all_path_sav_stix[this_f],/ver
    array_nsources = [array_nsources, n_fwdfit_sources]
  endfor
  ;; Total number of sources
  n_tot_sources = fix(total(array_nsources[rem_dup(array_nsources)]))
  
  
  ;;;;; Get the indices when nebins change the number of sources
  ind_change_nsources = rem_dup(array_nsources)
  nchanges = n_elements(ind_change_nsources)
  
  ;stop
  ;;;;; Restore the imaging sav files and create the variables
  n_sources = n_tot_sources ;n_max_sources
  all_vis            = []
  ;all_clean_maps     = []
  ;total_clean_flux   = []
  ;all_mem_maps       = []
  ;total_mem_flux     = []
  ;all_fwdfit_maps    = []
  all_energy_ranges     = fltarr(2,n_ebins)
  ind_more_sources = 33
  ;; Distinguish the case of a single source or multiple sources
  if n_sources eq 1 then begin
    
    all_fwdfit_fluxes     = fltarr(n_ebins)
    all_fwdfit_err_fluxes = fltarr(n_ebins)
    all_fwdfit_pos        = fltarr(2,n_ebins)
    all_fwdfit_err_pos    = fltarr(2,n_ebins)
    
    ;; Loop on all energies
    for this_en=0,n_ebins-1 do begin
      
      restore,all_path_sav_stix[this_en],/ver
      all_vis                      = [all_vis, vis]
      all_energy_ranges[*,this_en] = this_energy_range
      
      all_fwdfit_fluxes[this_en]     = out_param_fwdfit.srcflux
      all_fwdfit_err_fluxes[this_en] = out_sigma_fwdfit.srcflux
      all_fwdfit_pos[*,this_en]      = [out_param_fwdfit.srcx, out_param_fwdfit.srcy]
      all_fwdfit_err_pos[*,this_en]  = [out_sigma_fwdfit.srcx, out_sigma_fwdfit.srcy]
      
    endfor
    
    
    ;;;;; Scale the fluxes at 1 AU
    all_fwdfit_fluxes_so = all_fwdfit_fluxes
    all_fwdfit_err_fluxes_so = all_fwdfit_err_fluxes
    total_flux_fwdfit = all_fwdfit_fluxes * distance^2
    total_flux_fwdfit_err = all_fwdfit_err_fluxes * distance^2
    

    ;;;;; Define the fwdfit axes for the plot
    e_axis_fwdfit_plot = [all_energy_ranges[0,0],mean(all_energy_ranges,dim=1),all_energy_ranges[1,-1]]
    e_axis_fwdfit = mean(all_energy_ranges,dim=1)
    total_flux_fwdfit_plot = [total_flux_fwdfit[0],total_flux_fwdfit,total_flux_fwdfit[-1]]
    
    
  endif else if (n_sources gt 1 and n_elements(ind_change_nsources) eq 1) then begin
    
    all_fwdfit_fluxes     = fltarr(2,n_ebins)
    all_fwdfit_err_fluxes = fltarr(2,n_ebins)
    all_fwdfit_pos        = fltarr(2,n_sources,n_ebins)
    all_fwdfit_err_pos    = fltarr(2,n_sources,n_ebins)

    ;; Loop on all energies
    for this_en=0,n_ebins-1 do begin

      restore,all_path_sav_stix[this_en],/ver
      all_vis                      = [all_vis, vis]
      all_energy_ranges[*,this_en] = this_energy_range
      
      ;; Loop on all sources
      for this_s=0,n_fwdfit_sources-1 do begin
        
        all_fwdfit_fluxes[this_s,this_en]     = out_param_fwdfit[this_s].srcflux
        all_fwdfit_err_fluxes[this_s,this_en] = out_sigma_fwdfit[this_s].srcflux
        all_fwdfit_pos[*,this_s,this_en]      = [out_param_fwdfit[this_s].srcx, out_param_fwdfit[this_s].srcy]
        all_fwdfit_err_pos[*,this_s,this_en]  = [out_sigma_fwdfit[this_s].srcx, out_sigma_fwdfit[this_s].srcy]
        
      endfor

    endfor
    
    
    ;;;;; Scale the fluxes at 1 AU
    all_fwdfit_fluxes_so = all_fwdfit_fluxes
    all_fwdfit_err_fluxes_so = all_fwdfit_err_fluxes
    all_fwdfit_fluxes *= distance^2
    all_fwdfit_err_fluxes *= distance^2


    ;;;;; Calculate the total flux and the related error
    total_flux_fwdfit = total(all_fwdfit_fluxes,1)
    total_flux_fwdfit_err = sqrt(total(all_fwdfit_err_fluxes^2,1))


    ;;;;; Define the fwdfit axes for the plot
    e_axis_fwdfit_plot = [all_energy_ranges[0,0],mean(all_energy_ranges,dim=1),all_energy_ranges[1,-1]]
    e_axis_fwdfit = mean(all_energy_ranges,dim=1)
    total_flux_fwdfit_plot = [total_flux_fwdfit[0],total_flux_fwdfit,total_flux_fwdfit[-1]]
    
    
  endif else begin
    
    
    all_fwdfit_fluxes     = fltarr(n_sources,n_ebins)
    all_fwdfit_err_fluxes = fltarr(n_sources,n_ebins)
    all_fwdfit_pos        = fltarr(2,n_sources,n_ebins)
    all_fwdfit_err_pos    = fltarr(2,n_sources,n_ebins)
    
    for this_en=0,ind_change_nsources[-1]-1 do begin
      restore,all_path_sav_stix[this_en],/ver
      all_vis            = [all_vis, vis]
      ;all_clean_maps     = [all_clean_maps, clean_map[0]]
      ;;;total_clean_flux   = [total_clean_flux, flux_clean]
      ;all_mem_maps       = [all_mem_maps, memge_map]
      ;;;total_mem_flux     = [total_mem_flux, flux_mem]
      ;all_fwdfit_maps    = [all_fwdfit_maps, fwdfit_map]
      ;all_fwdfit_fluxes[*,this_en] = out_param_fwdfit.srcflux
      ;all_fwdfit_err_fluxes[*,this_en] = out_sigma_fwdfit.srcflux
      if n_fwdfit_sources gt 1 then ind_more_sources = min([ind_more_sources, this_en])
      for this_s=0,n_fwdfit_sources-1 do begin
        ;stop
        if this_s+1 le n_fwdfit_sources then begin
          if n_fwdfit_sources gt 1 then begin
            all_fwdfit_fluxes[this_s,this_en] = out_param_fwdfit[this_s].srcflux
            all_fwdfit_err_fluxes[this_s,this_en] = out_sigma_fwdfit[this_s].srcflux
            all_fwdfit_pos[*,this_s,this_en] = [out_param_fwdfit[this_s].srcx, out_param_fwdfit[this_s].srcy]
            all_fwdfit_err_pos[*,this_s,this_en] = [out_sigma_fwdfit[this_s].srcx, out_sigma_fwdfit[this_s].srcy]
          endif else begin
            all_fwdfit_fluxes[this_s,this_en] = out_param_fwdfit.srcflux
            all_fwdfit_err_fluxes[this_s,this_en] = out_sigma_fwdfit.srcflux
            all_fwdfit_pos[*,this_s,this_en] = [out_param_fwdfit.srcx, out_param_fwdfit.srcy]
            all_fwdfit_err_pos[*,this_s,this_en] = [out_sigma_fwdfit.srcx, out_sigma_fwdfit.srcy]
          endelse
        endif else begin
          all_fwdfit_fluxes[this_s,this_en] = out_param_fwdfit.srcflux
          all_fwdfit_err_fluxes[this_s,this_en] = out_sigma_fwdfit.srcflux
          all_fwdfit_pos[*,this_s,this_en] = [0.,0.]
          all_fwdfit_err_pos[*,this_s,this_en] = [0., 0.]
        endelse
      endfor
      all_energy_ranges[*,this_en] = this_energy_range
    endfor

    add_this = array_nsources[ind_change_nsources[-1]-1]    ; at the point the number of sources increase, then this increases too
    for this_en=ind_change_nsources[-1],n_ebins-1 do begin
      restore,all_path_sav_stix[this_en];,/ver
      all_vis            = [all_vis, vis]
      ;all_clean_maps     = [all_clean_maps, clean_map[0]]
      ;;;total_clean_flux   = [total_clean_flux, flux_clean]
      ;all_mem_maps       = [all_mem_maps, memge_map]
      ;;;total_mem_flux     = [total_mem_flux, flux_mem]
      ;all_fwdfit_maps    = [all_fwdfit_maps, fwdfit_map]
      for this_s=add_this,n_fwdfit_sources+add_this-1 do begin
        all_fwdfit_fluxes[this_s,this_en] = out_param_fwdfit[this_s-add_this].srcflux
        all_fwdfit_err_fluxes[this_s,this_en] = out_sigma_fwdfit[this_s-add_this].srcflux
        all_fwdfit_pos[*,this_s,this_en] = [out_param_fwdfit[this_s-add_this].srcx, out_param_fwdfit[this_s-add_this].srcy]
        all_fwdfit_err_pos[*,this_s,this_en] = [out_sigma_fwdfit[this_s-add_this].srcx, out_sigma_fwdfit[this_s-add_this].srcy]
      endfor
      all_energy_ranges[*,this_en] = this_energy_range
    endfor
    
    ;;;;; Scale the fluxes at 1 AU
    all_fwdfit_fluxes_so = all_fwdfit_fluxes
    all_fwdfit_err_fluxes_so = all_fwdfit_err_fluxes
    all_fwdfit_fluxes *= distance^2
    all_fwdfit_err_fluxes *= distance^2


    ;;;;; Calculate the total flux and the related error
    total_flux_fwdfit = total(all_fwdfit_fluxes,1)
    total_flux_fwdfit_err = sqrt(total(all_fwdfit_err_fluxes^2,1))


    ;;;;; Define the fwdfit axes for the plot
    e_axis_fwdfit_plot = [all_energy_ranges[0,0],mean(all_energy_ranges,dim=1),all_energy_ranges[1,-1]]
    e_axis_fwdfit = mean(all_energy_ranges,dim=1)
    total_flux_fwdfit_plot = [total_flux_fwdfit[0],total_flux_fwdfit,total_flux_fwdfit[-1]]
    
  endelse
  
    
  ;IDL> help, all_fwdfit_fluxes
  ;FLOAT     = Array[2, 3]                     
  ;          = [number of sources, number of energy bins]
  ;IDL> help, all_fwdfit_pos
  ;ALL_FWDFIT_POS  FLOAT     = Array[2, 2, 3]  
  ;                          = [x and y, number of sources, number of energy bins]
  
  
  ;;;;; For the time_shift, check the convention used for imaging.
  this_time_shift = time_shift
  if earth_ut eq 0 then this_time_shift = 0
  this_time_range = stx_time2any(fwdfit_map.time_range) + this_time_shift
  
  
  ;;;;; Check if the folder with the specfile exists
  folder_specfile = path_data_folder+'ospex_fits'
  if not file_test(folder_specfile,/dir) then file_mkdir, folder_specfile
  
  ;;;;; Change working directory
  this_wd = curdir()
  cd,folder_specfile
  folder_specfile += folder_delimiter
  
  
  set_logenv, 'OSPEX_NOINTERACTIVE', '1'
  
  
  ;;;;; Here we have to extract the observed spectrum
  time_shift = this_time_shift
  stx_convert_pixel_data, $
    fits_path_data = path_sci_file,$
    fits_path_bk = path_bkg_file, $
    distance = distance, $
    time_shift = time_shift, $
    flare_location = xy_flare_stix, $
    ospex_obj = ospex_obj, $
    background_data = background_data, $
    /sav_srm ; --> this keyword is not officially implemented!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  wdelete
 

  ;;;; Find the OSPEX spectrum file
  flspec = findfile('stx_spectrum*.fits')

  
  ;;;;;; Out of the spectrum file extract the observed spectrum
  !null = mrdfits(flspec, 0, primary_header)
  out_ospex = mrdfits(flspec, 1, ospex_header)
  out_energy = mrdfits(flspec, 2, energy_header)
  
  ;; Detector area in cm^2 at 1 AU
  det_area = sxpar(ospex_header,'GEOAREA')
  
  ;; Energy axis
  e_min_ospex = out_energy.e_min
  e_max_ospex = out_energy.e_max
  e_bin_width = e_max_ospex - e_min_ospex
  
  ;; Time axis
  time_shift = this_time_shift
  stx_read_pixel_data_fits_file, path_sci_file, time_shift, t_axis = t_axis, /silent
  t_axis_ospex = stx_time2any(t_axis.mean)
  t_ospex_start = stx_time2any(t_axis.time_start)
  t_ospex_end = stx_time2any(t_axis.time_end)
  tot_duration = total(t_axis.duration)
  
  ;; Find the indices for the time range
  ind_low_time = (where(abs(this_time_range[0]-t_ospex_start) eq min(abs(this_time_range[0]-t_ospex_start))))
  ind_high_time = (where(abs(this_time_range[1]-t_ospex_end) eq min(abs(this_time_range[1]-t_ospex_end))))
  
  ;; Observed STIX spectrum (s-1 cm-2 keV-1)
  ; OSPEX gives back only the count rate (counts/sec)
  ospex_rate = out_ospex.rate
  ospex_err  = out_ospex.stat_err
  if ind_low_time eq ind_high_time then begin
    tmp_observed_spectrum = ospex_rate[*,ind_low_time:ind_high_time] / e_bin_width / det_area
    ;tmp_observed_spectrum_error = sqrt(total(ospex_err[*,ind_low_time:ind_high_time]^2,2)) / e_bin_width / det_area
    tmp_observed_spectrum_error = ospex_err[*,ind_low_time:ind_high_time] / e_bin_width / det_area

  endif else begin
    tmp_observed_spectrum = mean(ospex_rate[*,ind_low_time:ind_high_time],dim=2) / e_bin_width / det_area
    ;tmp_observed_spectrum_error = sqrt(total(ospex_err[*,ind_low_time:ind_high_time]^2,2)) / e_bin_width / det_area
    tmp_observed_spectrum_error = mean(ospex_err[*,ind_low_time:ind_high_time],dim=2) / e_bin_width / det_area
  endelse
  
  ;; Observed background spectrum
  ; OSPEX gives back only the counts
  bk_counts = background_data.COUNTS
  bk_error = background_data.ERROR
  tmp_observed_bk = total(bk_counts,2) / e_bin_width / det_area / tot_duration
  tmp_observed_bk_error = mean(bk_error,dim=2) / e_bin_width / det_area / tot_duration + 0.05 * tmp_observed_bk
  
  ;; The observed spectrum has to be compared with the same energy binning as the fwdfit
  e_min_ospex_axis = []
  e_max_ospex_axis = []
  observed_spectrum = []
  observed_spectrum_error = []
  observed_bk_spectrum = []
  observed_bk_spectrum_error = []
  for this_en=0,n_ebins-1 do begin
    ter = all_energy_ranges[*,this_en]
    ile = (where(abs(ter[0]-e_min_ospex) eq min(abs(ter[0]-e_min_ospex))))[0]
    ihe = (where(abs(ter[1]-e_max_ospex) eq min(abs(ter[1]-e_max_ospex))))[-1]
    
    e_min_ospex_axis = [e_min_ospex_axis, ter[0]]
    e_max_ospex_axis = [e_max_ospex_axis, ter[1]]
    observed_spectrum = [observed_spectrum, mean(tmp_observed_spectrum[ile:ihe])]
    ;observed_spectrum_error = [observed_spectrum_error, sqrt(total(tmp_observed_spectrum_error[ile:ihe]^2))]
    observed_spectrum_error = [observed_spectrum_error, mean(tmp_observed_spectrum_error[ile:ihe])]
    observed_bk_spectrum = [observed_bk_spectrum, mean(tmp_observed_bk[ile:ihe])]
    observed_bk_spectrum_error = [observed_bk_spectrum_error, mean(tmp_observed_bk_error[ile:ihe])]
  end
  observed_spectrum_plot = [observed_spectrum[0], observed_spectrum, observed_spectrum[-1]]
  observed_bk_spectrum_plot = [observed_bk_spectrum[0], observed_bk_spectrum, observed_bk_spectrum[-1]]
  e_axis_ospex = (e_min_ospex_axis + e_max_ospex_axis) / 2
  e_axis_ospex_plot = [e_min_ospex_axis[0], e_axis_ospex, e_max_ospex_axis[-1]]
  
  
  ;;;;; Calculate the residuals as (observed - fwdfit) / error_fwdfit
  residuals = (observed_spectrum - total_flux_fwdfit) / total_flux_fwdfit_err
  residuals_plot = [residuals[0], residuals, residuals[-1]]
  
  
  ;;;;; Percent difference
  perc_diff = observed_spectrum / total_flux_fwdfit
  perc_diff_plot = [perc_diff[0], perc_diff, perc_diff[-1]]
  
  
  ;;;;; Plot the spectra as a ps file
  ;; Cosmetics plot
  xr = [all_energy_ranges[0,0],all_energy_ranges[1,-1]]
  yr = [min(total_flux_fwdfit)/10.,max(total_flux_fwdfit)*10.]
  
  ;; Parameters page and plot
  page_width=8.5
  page_height=11.5
  xsize=4.5
  ysize=6.5
  
  chart = 3
  ch_sz = 0.9
  
  !p.thick = 3
  !x.thick = 3
  !y.thick = 3
  !p.charthick = 2
  !p.charsize = ch_sz
  !x.charsize = ch_sz
  !y.charsize = ch_sz
  !p.ticklen = 0.04
  
  ;; Position plots
  ; Positions spectra
  p1 = [0.15 , 0.27, 0.91 , 0.95]
  ; Position residuals (observed - fwdfit) / error_fwdfit
  p11 = [0.15 , 0.06, 0.91 , 0.27]
  
  ;; Setup the postscript
  mydevice = !D.NAME
  set_plot,'PS'
  !p.font=0
  device,set_font='Helvetica'    ;HELVETICA
  device ,/isolatin1
  device,filename=ps_flnm_spectra, landscape=0,preview=1,/encapsulated
  device,xsize=xsize,ysize=ysize,/inches,font_index=6,/color,bits_per_pixel=24
  !p.multi=[0,2,1]
  
  loadct,0
  
  if this_time_shift eq 0 then begin 
    time_interval = time_range_so
    txt_ref_time = ' SolarOrbiter-UT'
  endif else begin
    time_interval = time_range_earth
    txt_ref_time = ' Earth-UT'
  endelse
  this_title = anytim(time_interval[0],/stime,/date_o) + ' ' + $
    STRMID(atime( time_interval[0]), 10, 8)+'-'+STRMID(atime(time_interval[1]), 10, 8)+txt_ref_time
  
  plot,e_axis_ospex,observed_spectrum,/xlog,/ylog,/nodata,xtitle='',$
    ytitle='!6Count Flux [cts s!U-1!N cm!U-2!N keV!U-1!N]',xst=1,xr=xr, $
    position=p1,yst=1,yr=yr,$; ytickformat='exponentlab',$
    title=this_title, xtickname=make_array(2,/string,value=' '), $
    CLIP=p1, /NORM, NOCLIP=1,chars=ch_sz-0.1
  
  ; Observed flare spectrum
  oplot,e_axis_ospex_plot,observed_spectrum_plot,psym=10,thick=5,color =cgcolor(observed_color),linest=0
  errplot,e_axis_ospex,(observed_spectrum-observed_spectrum_error)>1.d-15, observed_spectrum+observed_spectrum_error,$
    thick=2,color=0,CLIP=10.^[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]],NOCLIP=0,width=0
  
  ; Observed bk spectrum
  oplot,e_axis_ospex_plot,observed_bk_spectrum_plot,psym=10,thick=2.5,color =cgcolor(observed_color),linest=1
  errplot,e_axis_ospex,(observed_bk_spectrum-observed_bk_spectrum_error)>1.d-15, observed_bk_spectrum+observed_bk_spectrum_error,$
    thick=1,color=0,CLIP=10.^[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]],NOCLIP=0,width=0
  
  ; fwdfit spectrum
  oplot,e_axis_fwdfit_plot,total_flux_fwdfit_plot,thick=chart,color=cgcolor('slate gray'),psym=10,NOCLIP=0
  errplot,e_axis_fwdfit,(total_flux_fwdfit-total_flux_fwdfit_err)>1.d-15,total_flux_fwdfit+total_flux_fwdfit_err,$
    thick=2,color=cgcolor('slate gray'),CLIP=10.^[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]],NOCLIP=0,width=0
    
  text_legend = ['Observed flare spectrum', 'Observed bk spectrum', 'Total FWDFIT spectrum']
  color_legend = [cgcolor('black'), cgcolor('black'), cgcolor('slate gray')]
  chars_legend = [ch_sz, ch_sz, ch_sz]
  linest_legend = [0, 1, 0]
  
;  ; clean spectrum
;  ytit_res = ''
;  perc_diff_clean = total_clean_flux * 0
;  if total(total_clean_flux) gt 1. then begin
;    oplot,e_axis_fwdfit_plot,[total_clean_flux[0],total_clean_flux,total_clean_flux[-1]],thick=chart,color=cgcolor('blue'),psym=10,NOCLIP=0
;    
;    perc_diff_clean = observed_spectrum / total_clean_flux
;    perc_diff_clean_plot = [perc_diff_clean[0], perc_diff_clean, perc_diff_clean[-1]]
;    
;    text_legend = [text_legend, 'Total CLEAN spectrum']
;    color_legend = [color_legend, cgcolor('blue')]
;    chars_legend = [chars_legend, ch_sz]
;    linest_legend = [linest_legend, 0]
;    ytit_res += ',clean'
;  endif
;  
;  ; mem spectrum
;  perc_diff_mem = total_mem_flux * 0
;  if total(total_mem_flux) gt 1. then begin
;    oplot,e_axis_fwdfit_plot,[total_mem_flux[0],total_mem_flux,total_mem_flux[-1]],thick=chart,color=cgcolor('sea green'),psym=10,NOCLIP=0
;    
;    perc_diff_mem = observed_spectrum / total_mem_flux
;    perc_diff_mem_plot = [perc_diff_mem[0], perc_diff_mem, perc_diff_mem[-1]]
;
;    text_legend = [text_legend, 'Total MEM_GE spectrum'] 
;    color_legend = [color_legend, cgcolor('sea green')]
;    chars_legend = [chars_legend, ch_sz]
;    linest_legend = [linest_legend, 0]
;    ytit_res += ',mem_ge'
;  endif
  
  
  if n_sources gt 1 then begin
    for bb=0,nchanges-1 do begin
      for ii=0,n_sources-1 do begin
        imin = ind_change_nsources[bb]
        if bb eq nchanges-1 then imax = n_ebins-1 else imax = ind_change_nsources[bb+1]-1
        if ALL_FWDFIT_POS[0,ii,imin] ne 0. and ALL_FWDFIT_POS[1,ii,imin] ne 0 then begin
          this_flux = [all_fwdfit_fluxes[ii,imin],transpose(all_fwdfit_fluxes[ii,imin:imax]),all_fwdfit_fluxes[ii,imax]]
          this_eflux_p = transpose(all_fwdfit_fluxes[ii,imin:imax]) + transpose(all_fwdfit_err_fluxes[ii,imin:imax])
          this_eflux_m = transpose(all_fwdfit_fluxes[ii,imin:imax]) - transpose(all_fwdfit_err_fluxes[ii,imin:imax])
          this_eaxis = [all_energy_ranges[0,imin],mean(all_energy_ranges[*,imin:imax],dim=1),all_energy_ranges[1,imax]]
          this_eaxis_err = mean(all_energy_ranges[*,imin:imax],dim=1)
          
          oplot, this_eaxis, this_flux, color=cgcolor(color_srcs[ii+bb]), thick=chart, psym=10,NOCLIP=0
          errplot,this_eaxis_err,this_eflux_m>1.d-15,this_eflux_p,thick=2,color=cgcolor(color_srcs[ii+bb]),$
            CLIP=10.^[!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]],NOCLIP=0,width=0
    
          this_pos_txt = 'x='+num2str(mean(ALL_FWDFIT_POS[0,ii,imin:imax]),format='(I10.0)')+$
            '; y='+num2str(mean(ALL_FWDFIT_POS[1,ii,imin:imax]),format='(I10.0)')
          text_legend = [text_legend, this_pos_txt]
          color_legend = [color_legend, cgcolor(color_srcs[ii+bb])]
          chars_legend = [chars_legend, ch_sz]
          linest_legend = [linest_legend, 0]
        endif
      endfor
    endfor
  endif
  
  ssw_legend,text_legend,/top,/right,chars=ch_sz,colors=color_legend,box=0,linest=linest_legend,thick=chart
  
  plot,e_axis_fwdfit_plot,residuals_plot,/xlog,position=p11,xst=1,psym=10,xr=xr,xtitle='Energy [keV]', $
    ytitle='(obs - fwdfit)/err_fwdfit' , yst=9, yr=[-10,10], $
    ;yticks=7,ytickv=[-5,-3,-1, 0, 1, 3, 5],ytickname=['-5','-3','-1','0','1','3','5'], $
    xticklen=0.1,yticklen=0.03,chars=ch_sz-0.1, thick=4
  
  oplot, e_axis_fwdfit_plot, e_axis_fwdfit_plot*0., color=cgcolor('black'), thick=2, linest=1,NOCLIP=0
  
  cgAxis, YAxis=1, ystyle=1, YRange=[0.5,1.5], $
    ytitle='obs / fwdfit', /Save, color=cgcolor('slate gray'),chars=ch_sz-0.1
  
  oplot,e_axis_fwdfit_plot, perc_diff_plot, color=cgcolor('slate gray'), thick=4, psym=10,NOCLIP=0
;  if total(perc_diff_clean) gt 0.1 then oplot,e_axis_fwdfit_plot, perc_diff_clean_plot, color=cgcolor('blue'), thick=4, psym=10,NOCLIP=0
;  if total(perc_diff_mem) gt 0.1 then oplot,e_axis_fwdfit_plot, perc_diff_mem_plot, color=cgcolor('sea green'), thick=4, psym=10,NOCLIP=0
  oplot,e_axis_fwdfit_plot, perc_diff_plot*0.+1.2, color=cgcolor('slate gray'), thick=2, linest=1,NOCLIP=0
  oplot,e_axis_fwdfit_plot, perc_diff_plot*0.+0.8, color=cgcolor('slate gray'), thick=2, linest=1,NOCLIP=0
  
  
  device,/close_file
  set_plot,mydevice
  !p.thick = 1
  !x.thick = 1
  !y.thick = 1
  !p.charthick = 1
  !p.charsize = 1
  !x.charsize = 1
  !y.charsize = 1
  !p.ticklen = 0.02
  
  
  
  ;;obj_destroy, ospex_obj
  
  set_logenv, 'OSPEX_NOINTERACTIVE', '0'
  
;  ;;;;; Restore the SRM IDL file
;  path_srm_save = findfile('stx_srm*IDL-save.sav')
;  restore,path_srm_save,/ver
;  
;  cd,this_wd
;  
;  ;;;;; Store the fwdfit fluxes and the observed spectrum in a save file
;
;  ;;;;; Create the output structure and save it in a save file
;  total_fwdfit_flux = total_flux_fwdfit
;  total_fwdfit_flux_err = total_flux_fwdfit_err
;  readme = 'fwdfit flux already scaled to 1 AU'
;  flnm_fluxes = path_data_folder + 'fwdfit-and-observed-spectra.sav'
;  
;  ff_src = all_fwdfit_fluxes
;  ff_src_err = all_fwdfit_err_fluxes
;  
;  flux_str = {$
;    info            : 'Struct with fwdfit fluxes', $
;    readme          : readme, $
;    units           : 'cnts/s/keV/cm^2', $
;    obs             : observed_spectrum, $
;    obs_err         : observed_spectrum_error, $
;    ff_tot          : total_fwdfit_flux, $
;    ff_tot_err      : total_fwdfit_flux_err, $
;    ff_sources      : ff_src, $
;    ff_sources_err  : ff_src_err, $
;    e_min           : all_energy_ranges[0,*], $
;    e_max           : all_energy_ranges[1,*], $
;    time_range      : this_time_range, $
;    time_shift      : this_time_shift, $
;    det_area        : det_area, $
;    srm             : srm}
;  
;  save,filename=flnm_fluxes,flux_str,e_axis_ospex,residuals,perc_diff,observed_bk_spectrum,$
;    observed_bk_spectrum_error;,total_clean_flux,total_mem_flux
;
;
; OPTIONAL OUTPUTS:
;   flux_str           : structure containing the fluxes of the different sources.
;                        We note that obs stands for the spatially integrated spectrum
;                        and ff stands for the fwdfit (imaging) spectrum.
;                        It contains the following keywords:
;                          flux_str = {$
;                                  info            : 'Struct with fwdfit fluxes', $
;                                  readme          : readme, $
;                                  units           : 'cnts/s/keV/cm^2', $
;                                  obs             : observed_spectrum, $
;                                  obs_err         : observed_spectrum_error, $
;                                  ff_tot          : total_fwdfit_flux, $
;                                  ff_tot_err      : total_fwdfit_flux_err, $
;                                  ff_sources      : ff_src, $
;                                  ff_sources_err  : ff_src_err, $
;                                  e_min           : all_energy_ranges[0,*], $
;                                  e_max           : all_energy_ranges[1,*], $
;                                  time_range      : this_time_range, $
;                                  time_shift      : this_time_shift, $
;                                  det_area        : det_area, $
;                                  srm             : srm}
;
;   ospex_obj          : ospex object containing the spatially integrated spectrum
  

  if keyword_set(stop_here) then stop

end ; End of the script

