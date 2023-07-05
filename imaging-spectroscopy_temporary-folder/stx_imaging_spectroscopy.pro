;+
;
; PROJECT:
;   STIX
;
; NAME:
;   stx_imaging_spectroscopy
;
; PURPOSE:
;   Perform imaging spectroscopy with STIX. This procedures stores the results in a folder
;   and creates a save file with the results of the forward fitting algorithm.
;   Idea of the procedure: if you do not know a priori the configuration of the
;   reconstructed image (in the considered energy range), then at each energy bin you have 
;   to specify the number of sources. If you know it, by specifying configuration_fwdfit, you
;   do not have to select the sources each time.
;   By default: 
;       - if we have less than 1000 counts in a given energy range, it stops iterating.
;       - if configuration_fwdfit is not specified, the script assumes circular Gaussians. 
;         Set the keyword "/ellipse_shape" to use elliptical Gaussians
;       - if the location of the sources is not specified, then they are fitted with the
;         fwdfit algorithm. If the keyword "/select_location" is set, then the user can
;         select the location of the sources on the screen (for each energy bin) using the
;         cursor. If the keyword "/select_box_location" is set, then the user can select
;         boxes on the screen to constrain the location of the sources for the fwdfit algorithm.
;       - if the keyword "/select_box_location" is set, then the user can eventually fix the
;         location of some of the sources by double clicking on the same location. This allows
;         to fix the location of the sources that are known (e.g. footpoints) and fit with boxes 
;         the location of the sources that are not known (e.g. nonthermal loop-top sources).
;       - if the keyword "/select_box_location" is set and the user manually fixes the location
;         (by double clicking on the same location), then the shape will automatically change in
;         circular Gaussian
;
; CALLING SEQUENCE:
;   stx_imaging_spectroscopy, path_sci_file, path_bkg_file, aux_fits_file, time_range, energy_max_inversion
;
; INPUTS:
;   path_sci_file        : string containing the path of the science L1 fits file of the event
;
;   path_bkg_file        : string containing the path of the backgroung L1 fits file
;
;   path_aux_file        : string containing the path of the auxiliary fits file of the day
;   
;   time_range           : array containing the selected start and the end time 
;
;   energy_max_inversion : maximum energy to use for the inversion to calculate the regularized
;                          visibilities. If the keyword "/standard_vis" is set, then this is ignored.
;
; OPTIONAL OUTPUTS:
;   path_new_folder      : string containing the path of the newly created folder where the 
;                          results are stored
;
; OPTIONAL INPUTS AND KEYWORDS:
;   standard_vis         : if set, then the "standard visibilities" are used. In this case, the
;                          energy_max_inversion is ignored.
;
;   energy_min_inversion : minimum energy to use for the inversion to calculate the regularized
;                          visibilities. If the keyword "/standard_vis" is set, then this is ignored.
;                          Default is 9 keV (we cannot go lower with the regularized visibilities).
;
;   configuration_fwdfit : array containing the configuration of the forward fitting algorithm
;
;   source_loc           : array containing the location of the source (if known a priori).
;                          This has to have the same number of elements as configuration_fwdfit.
;                          This keyword is mutually exclusive with /select_location.
;
;   source_fwhm          : array containing the fwhm of the source (if known a priori).
;                          This has to have the same number of elements as configuration_fwdfit.
;
;   min_fwhm             : it allows to constrain the minimum size of the FWHM of the fwdfit source.
;                          A good choice is to set it equals the resolution of the finest sc used.
;                          If this is a single value (and not an array) the same min FWHM will be used
;                          for all the sources. If it is an array, then it has to have the same number
;                          of elements as configuration_fwdfit.
;
;   max_fwhm             : the same as min_fwhm, but for the maximum FWHM.
;
;   uncertainty          : set to 0 to avoid the computation of the parameters uncertainty of the fwdfit 
;                          algorithm (confidence strip approach). Default is 1 (i.e., compute the 
;                          uncertainties).
;
;   path_sav_folder      : folder where to store the sav files. A new folder will be created inside
;                          this folder, with the date and integration time of the observation as name.
;                          If not specified, then the current directory is used.
;                          If the folder with the date and integration time already exists, then the
;                          new results will be stored inside the same folder (the files with the same
;                          energy bins will be overwritten!)
;
;   suffix_folder        : suffix to add to the folder name when it is created
;
;   energy_range         : array containing the energy range to image. This is a two-elements array,
;                          and the procedure will image all native energies in this range.
;
;   energy_low           : array containing the lower energy edges to image. This is useful if you want
;                          to image energy bins that are different from the native ones. If this is set,
;                          then energy_high has to be set as well.
;
;   energy_high          : the same as energy_low, but for the upper energy edges.
;
;   earth_ut             : if set, then the time_range (given as input) is in Earth UT, otherwise it is 
;                          in Solar Orbiter UT (default).
;
;   stop_here            : if set, then the procedure stops after every iteration in energy (for debugging)
;
;   min_counts           : minimum number of counts to continue iterating. If the number of counts in a given
;                          energy bin is lower than this value, then the procedure stops iterating.
;
;   contour_clean_box    : threshold to define the CLEAN box. The box is defined as the region where the
;                          BPROJ map is above this threshold times the maximum of the BPROJ map.
;
;   select_location      : if set, then the user can select the location of the sources on the screen (for
;                          each energy bin) using the cursor. If this is set, then source_loc is discarded.
;
;   select_box_location  : instead of fixing the location of the sources, it is possible to select boxes
;                          around the sources in order to constrain the location of the fwdfit sources
;
;   box_location         : instead of selecting the boxes on the screen, one can pass it and then it will
;                          be used for all the energy bins. One has to give the coordinates of the bottom left
;                          and top-right corners of the box as
;                                     [[bl_x, bl_y],          bottom-left (arcsec)
;                                      [tr_x, tr_y]]          top-right (arcsec)
;                          If one needs to give more boxes, then it has to be of the following format: [2,2,nb_sources]
;                          For example:
;                                     [[[bl_x1, bl_y1],[tr_x1, tr_y1]],          ; bottom-left and top-right(arcsec) box source 1
;                                      [[bl_x2, bl_y2],[tr_x2, tr_y2]]]          ; bottom-left and top-right(arcsec) box source 2
;
;   ellipse_shape        : if set, then the sources are elliptical Gaussians, otherwise they are circular.
;                          If the keyword select_box_location is set and the user manually fixes the location, 
;                          then the shape will automatically change in circular Gaussian
;   
;   subc_index           : array containing the subcollimator indices to use. If not specified, then all
;                          the TOP24 subcollimators are used.
;
;   pixels               : string containing the pixels to use. If not specified, then the TOP+BOT pixels are used.
;
;   imsize               : array containing the image size. If not specified, then [257,257] is used.
;
;   pixel                : array containing the pixel size. If not specified, then [1.,1.] is used.
;
;   niter                : number of iterations for the CLEAN algorithm. If not specified, then 200 is used.
;
;   gain                 : gain for the CLEAN algorithm. If not specified, then 0.05 is used.
;
;   nmap                 : number of maps to create for the CLEAN algorithm. If not specified, then 200 is used.
;
;   beam_width           : beam width for the CLEAN algorithm. If not specified, then 14.8 is used (resolution of
;                          the subcollimator 3).
;
; HISTORY: 
;   July 2023, Battaglia A. F. (FHNW & ETHZ), initial release
;
; CONTACT:
;   andrea.battaglia@fhnw.ch
;
;-

pro stx_imaging_spectroscopy, path_sci_file, path_bkg_file, aux_fits_file, time_range, energy_max_inversion,$     ; inputs
  ;; --- Optional inputs and keywords
  standard_vis = standard_vis, $
  energy_min_inversion = energy_min_inversion, $
  configuration_fwdfit = configuration_fwdfit, $
  source_loc = source_loc, $
  source_fwhm = source_fwhm, $
  min_fwhm = min_fwhm, $
  max_fwhm = max_fwhm, $
  uncertainty = uncertainty, $
  path_sav_folder = path_sav_folder, $
  suffix_folder = suffix_folder, $
  energy_range = energy_range, $
  energy_low = energy_low, $
  energy_high = energy_high, $
  earth_ut = earth_ut, $
  stop_here = stop_here, $
  min_counts = min_counts, $
  contour_clean_box = contour_clean_box, $
  select_location = select_location, $
  select_box_location = select_box_location, $
  box_location = box_location, $
  ellipse_shape = ellipse_shape, $
  ;; --- Standard imaging inputs
  subc_index = subc_index, $
  pixels = pixels, $
  imsize = imsize, $
  pixel = pixel, $
  niter = niter, $
  gain = gain, $
  nmap = nmap, $
  beam_width = beam_width, $
  ;; --- Optional output
  path_new_folder = path_new_folder
  
  
  ;;;;; Cosmetics plot
  ch_sz = 3
  grid = 5
  th_sz_txt = 1.5
  
  
  ;;;;; Default parameters
  default, path_sav_folder, curdir()
  default, energy_range, [4,28]
  default, energy_min_inversion, 9
  default, suffix_folder, ''
  default, subc_index, stx_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c','6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c'])
  default, pixels, 'TOP+BOT'
  default, imsize, [257,257]
  default, pixel, [1.,1.]
  default, niter, 200
  default, gain, 0.05
  default, nmap, 200
  default, beam_width, 14.8
  default, contour_clean_box, 0.5
  default, min_counts, 1d3
  default, uncertainty, 1
  default, energy_max_inversion, 0.
  if keyword_set(earth_ut) then earth_ut = 1 else earth_ut = 0

  
  
  ;;;;; Some error messages
  if (keyword_set(energy_low) and not keyword_set(energy_high)) or (keyword_set(energy_high) and not keyword_set(energy_low)) then begin
    message,'---> To manually define the energy edges, both energy_high and energy_low have to be defined.'
  endif
  if (keyword_set(source_loc) or keyword_set(select_location)) and (keyword_set(select_box_location) or keyword_set(box_location)) then begin
    message,'---> It is not possible to define a box (/select_box_location or /box_location) for the locations and at the same time fix them (/select_location or /source_loc). Please, review your keywords.'
  endif
  if not keyword_set(standard_vis) and energy_max_inversion eq 0 then begin
    print,''
    print,'By default, regularized visibilities are used. Therefore, energy_max_inversion has to be defined!'
    message,'---> If you want to use the standard visibilities, then set the keyword /standard_vis'
  endif

  
  ;;;;; Some variables that we need to initialize here
  old_nsources = 0.
  all_fixedpos = []
  all_fixedfwhm = []
  all_minfwhm = []
  all_maxfwhm = []
  all_boxes = []
  all_box_or_fix = []
  
  
  ;;;;; Get the folder delimiter (different for different OS)
  folder_delimiter = path_sep()
  
  
  ;;;;; Check if path_sav_folder is effectively a folder
  if not path_sav_folder.endswith(folder_delimiter) then path_sav_folder = path_sav_folder.insert(folder_delimiter,path_sav_folder.strlen())


  ;;;;; Get the light travel time
  stx_get_header_corrections,path_sci_file,time_shift=time_shift,distance=distance


  ;;;;; Get the energy axis in the file
  stx_read_pixel_data_fits_file, path_sci_file, 0., energy_str = energy_str, $
    energy_header = energy_header, e_axis = e_axis, t_axis = t_axis_so, use_discriminators = 0
  t_st_so = stx_time2any(t_axis_so.time_start)
  t_en_so = stx_time2any(t_axis_so.time_end)


  ;;;;; Conversion of the times (in case earth_ut is set)
  if keyword_set(earth_ut) then begin

    time_range_so_user = anytim(time_range) - time_shift

    ind_st_so = where(abs(time_range_so_user[0]-t_st_so) eq min(abs(abs(time_range_so_user[0]-t_st_so))))
    ind_en_so = where(abs(time_range_so_user[1]-t_en_so) eq min(abs(abs(time_range_so_user[1]-t_en_so))))

    time_range_so = [t_st_so[ind_st_so], t_en_so[ind_en_so]]
    time_range_earth = time_range_so + time_shift
    tmp_time = anytim(time_range_earth[0],/ex,/trunc,/time_o)
    this_time_start = num2str(tmp_time[0])+num2str(tmp_time[1])+num2str(tmp_time[2])
    tmp_time = anytim(time_range_earth[1],/ex,/trunc,/time_o)
    this_time_end = num2str(tmp_time[0])+num2str(tmp_time[1])+num2str(tmp_time[2])
    text_ref_time = 'Earth-UT'

  endif else begin

    time_range_so_user = anytim(time_range)

    ind_st_so = where(abs(time_range_so_user[0]-t_st_so) eq min(abs(abs(time_range_so_user[0]-t_st_so))))
    ind_en_so = where(abs(time_range_so_user[1]-t_en_so) eq min(abs(abs(time_range_so_user[1]-t_en_so))))

    time_range_so = [t_st_so[ind_st_so], t_en_so[ind_en_so]]
    time_range_earth = time_range_so + time_shift
    
    
    tmp_time = anytim(time_range_so[0],/ex,/trunc,/time_o)
    this_time_start = num2str(tmp_time[0])+num2str(tmp_time[1])+num2str(tmp_time[2])
    tmp_time = anytim(time_range_so[1],/ex,/trunc,/time_o)
    this_time_end = num2str(tmp_time[0])+num2str(tmp_time[1])+num2str(tmp_time[2])
    text_ref_time = 'SolarOrbiter-UT'

  endelse
  this_date = anytim(time_range_earth[0],/ccsds,/date_o)
  
  
  ;;;;; Create aux data
  aux_data = stx_create_auxiliary_data(aux_fits_file, time_range_so)
  
  
  ;;;;; Estimate flare location, if not already specified
  if not keyword_set(source_loc) then begin
    stx_estimate_flare_location, path_sci_file, time_range_so, aux_data, flare_loc=flare_loc, path_bkg_file=path_bkg_file
    xy_flare_stix = flare_loc
    mapcenter = flare_loc
  endif else begin
    xy_flare_stix = source_loc[*,0]
    mapcenter = xy_flare_stix
  endelse
  
  
  ;;;;; Coordinate transformaion
  ; from Helioprojective Cartesian to STIX coordinate frame
  mapcenter_stix = stx_hpc2stx_coord(mapcenter, aux_data)
  xy_flare_stix  = stx_hpc2stx_coord(xy_flare_stix, aux_data)
  
  ;;;;; Distinguish the case of standard or regularized visibilities
  if keyword_set(standard_vis) then begin

    ;;;;; String for the type of visibility used
    vis_type = 'standard visibilities'

    ;;;;; Define the energy range
    erange_low = e_axis.low
    erange_high = e_axis.high
    ind_elow = where(abs(energy_range[0]-erange_low) eq min(abs(energy_range[0]-erange_low)))
    ind_ehigh = where(abs(energy_range[1]-erange_high) eq min(abs(energy_range[1]-erange_high)))
    erange_low = erange_low[ind_elow:ind_ehigh]
    erange_high = erange_high[ind_elow:ind_ehigh]

    ;;;;; Number of energies to image
    if keyword_set(energy_low) then n_energies = n_elements(energy_low) else n_energies = n_elements(erange_low)

  endif else begin

    ;;;;; String for the type of visibility used
    vis_type = 'regularized visibilities'

    ;;;;; Define the energy range
    erange_low_all = e_axis.low
    erange_high_all = e_axis.high
    ind_elow = where(abs(energy_range[0]-erange_low_all) eq min(abs(energy_range[0]-erange_low_all)))
    ind_ehigh = where(abs(energy_range[1]-erange_high_all) eq min(abs(energy_range[1]-erange_high_all)))
    ind_elow_elec  = where(abs(energy_min_inversion-erange_low_all) eq min(abs(energy_min_inversion-erange_low_all)))
    ind_ehigh_elec = where(abs(energy_max_inversion-erange_high_all) eq min(abs(energy_max_inversion-erange_high_all)))
    erange_low = erange_low_all[ind_elow:ind_ehigh]
    erange_high = erange_high_all[ind_elow:ind_ehigh]


    ;;;;; Number of energies to image
    if keyword_set(energy_low) then begin
      n_energies = n_elements(energy_low) 
      this_elow = energy_low[0]
    endif else begin
      n_energies = n_elements(erange_low)
      this_elow = erange_low[0]
    endelse
    
    
    ;;;;; Number of energies to calculate the inversion (for the electron maps)
    erange_low_elec  = erange_low_all[ind_elow_elec:ind_ehigh_elec]
    erange_high_elec = erange_high_all[ind_elow_elec:ind_ehigh_elec]
    n_energies_elec  = n_elements(erange_low_elec)
    
    
    ;;;;; Loop over the energies to create the visibilities for each subcollimator at each energy
    ;; This will then be used for electron maps
    vis_tot = []
    for i=0,n_energies_elec-1 do begin
      ;;;;; Energy range
      this_energy_range = [erange_low_elec[i],erange_high_elec[i]]

      ;;;;; Construct the visibilities
      vis = stx_construct_calibrated_visibility(path_sci_file, time_range_so, this_energy_range, mapcenter_stix, subc_index=subc_index, $
        path_bkg_file=path_bkg_file, xy_flare=xy_flare_stix, no_small=no_small)
      
      ;;;;; Concatenate
      vis_tot = [vis_tot, vis]
    end


    ;;;;; Perform the inversion: compute the electron visibilities from STIX visibilities
    ;; Solar orbiter distance in cm
    dist_solo_sun = distance * 1.496d13
    stx_visibilities_inversion, vis_tot, dist_solo_sun, attenuator, reg_el_vis, orig_ph_vis, reg_ph_vis

  endelse
  
  ;;;;; Loop over all energies for imaging-spectroscopy
  for this_en = 0,n_energies-1 do begin
    
    
    ;;;;; Set the energy range
    if keyword_set(energy_low) then begin
      this_energy_range = [energy_low[this_en],energy_high[this_en]]
    endif else begin
      this_energy_range = [erange_low[this_en], erange_high[this_en]]
    endelse
    
    ;;;;; Distinguish the case of standard or regularized visibilities
    ;; to create the visibility structure
    if keyword_set(standard_vis) then begin

      ;;;;; Create visibility structure
      vis = stx_construct_calibrated_visibility(path_sci_file, time_range_so, this_energy_range, mapcenter_stix, $
        path_bkg_file=path_bkg_file, xy_flare=xy_flare_stix, no_small=no_small, $
        subc_index=subc_index)

      ;;;;; Get the number of counts
      ncounts = mean(vis.tot_counts - vis.tot_counts_bkg)

    endif else begin

      ;;;;; Find the indices of the regularized visibilities for the considered energy range
      loc = where(reg_ph_vis.energy_range[0] eq this_energy_range[0], count)
      vis = replicate(stx_visibility(),  count)
      vis = reg_ph_vis[loc]

      ;;;;; Create visibility structure (just for the # of counts)
      vis_orig = stx_construct_calibrated_visibility(path_sci_file, time_range_so, this_energy_range, mapcenter_stix, $
        path_bkg_file=path_bkg_file, xy_flare=xy_flare_stix, no_small=no_small, $
        subc_index=subc_index)

      ;;;;; Get the number of counts
      ncounts = mean(vis_orig.tot_counts - vis_orig.tot_counts_bkg)

    endelse
    
    ;;;;; If the number of counts are lower than min_counts, then stop iterating and go to the plotting procedure
    if ncounts lt min_counts then begin
      txt1e = num2str(this_energy_range[0],format='(I10.0)')
      txt2e = num2str(this_energy_range[1],format='(I10.0)')
      txtcn = num2str(min_counts,format='(I10.0)')
      print,''
      print,''
      print,'********************************************************************************'
      print,'Number of counts in the energy bin '+txt1e+'-'+txt2e+' keV lower than '+txtcn+'.'
      print,'No further images are created. Abort the loop'
      print,'********************************************************************************'
      print,''
      print,''
      break
    endif
    
    
    ;;;;; BPROJ
    bp_map = stx_bproj(vis, imsize, pixel, aux_data)
    
    
    ;;;;; CLEAN
    ;; Get the CLEAN box
    clean_box = where(rotate(bp_map.data,3) ge contour_clean_box*max(bp_map.data))
    ;; Do CLEAN
    clean_map = stx_vis_clean(vis, aux_data, niter=niter, image_dim=imsize[0], PIXEL=pixel[0], $
      gain=gain, nmap=nmap, /plot, beam_width=beam_width, clean_box=clean_box)
      
    
    ;;;;; MEM
    memge_map = stx_mem_ge(vis,imsize,pixel,aux_data,total_flux=max(abs(vis.obsvis)))
    
    
    ;;;;; Plot the maps and the number of counts
    txt1e = num2str(this_energy_range[0],format='(I10.0)')
    txt2e = num2str(this_energy_range[1],format='(I10.0)')
    print,''
    print,''
    print,'********************************************************************************'
    print,'Number of counts in this energy bin ('+txt1e+'-'+txt2e+' keV):'
    print,ncounts
    print,'********************************************************************************'
    print,''
    print,''
    
    window,10,xsize=1700,ysize=600
    !p.multi=[0,3,1]
    loadct,5,/sil
    plot_map,bp_map,/limb,grid=grid,chars=ch_sz
    plot_map,clean_map[0],/limb,grid=grid,chars=ch_sz
    plot_map,clean_map[0],/over,/perc,levels=[30,50],color=cgcolor('cyan'),thick=2
    plot_map,clean_map[0],/over,levels=[max(clean_map[2].data)],color=cgcolor('dark green'),thick=2
    xyouts,0.41,0.75,'Peak residual',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('dark green')
    xyouts,0.41,0.81,'[30,50]%',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('cyan')
    loadct,5,/sil
    plot_map,memge_map,/limb,grid=grid,chars=ch_sz
    plot_map,memge_map,/over,/perc,levels=[30,50],color=cgcolor('cyan'),thick=2
    xyouts,0.74,0.81,'[30,50]%',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('cyan')

    
    ;;;;; If configuration_fwdfit is not set, then select the number of sources
    if not keyword_set(configuration_fwdfit) then begin
      ;; Plot CLEAN again to select the sources on the screen
      loadct,5,/silent
      !p.multi = 0
      window,12,xsize=1300,ysize=1300;,xpos=2300
      plot_map,clean_map[4],/limb,grid=grid,chars=ch_sz
      plot_map,clean_map[4],/over,/perc,levels=[30,50],color=cgcolor('cyan'),thick=2
      ;plot_map,clean_map[4],/over,levels=[max(clean_map[2].data)],color=cgcolor('dark green'),thick=2
      ;xyouts,0.21,0.75,'Peak residual',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('dark green')
      xyouts,0.21,0.8,'[30,50]%',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('cyan')
      
      ;; Overplot the center of the previously selected sources
      if old_nsources gt 0 and keyword_set(source_loc) then begin
        for i=0,n_elements(all_fixedpos[0,*])-1 do plots,all_fixedpos[0,i],all_fixedpos[1,i],/data,psym=1,syms=5,color=cgcolor('Slate Gray'),thick=3
        for i=0,n_elements(old_fixedpos[0,*])-1 do plots,old_fixedpos[0,i],old_fixedpos[1,i],/data,psym=1,syms=5,color=cgcolor('magenta'),thick=5
        xyouts,0.21,0.77,'All previously selected locations',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('Slate Gray')
        xyouts,0.21,0.74,'Latest selected locations',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('magenta')
      endif

      ;; Overplot the previously selected boxes
      if old_nsources gt 0 and keyword_set(box_location) then begin
        plot_boxes, all_boxes, all_box_or_fix
        plot_boxes, old_boxes, box_or_fix, color='magenta'
        xyouts,0.21,0.77,'All previously selected boxes',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('Slate Gray')
        xyouts,0.21,0.74,'Latest selected boxes',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('magenta')
      endif
      
      ;; Print a message and ask the user to input the number of sources
      print,''
      print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      read, these_sources, prompt='Enter number of sources and press ENTER (0: abort): '
      n_fwdfit_sources = fix(these_sources)
      if n_fwdfit_sources eq 0 then break
      
      ;; Then, the configuration is n_fwdfit_sources times circular gaussians
      this_configuration_fwdfit = []
      this_shape = 'circle'
      if keyword_set(ellipse_shape) then this_shape = 'ellipse'
      for i=0,n_fwdfit_sources-1 do this_configuration_fwdfit = [this_configuration_fwdfit, this_shape]
    
    endif else begin
      ;; If the locations are known a priori, then
      n_fwdfit_sources = n_elements(configuration_fwdfit)
      this_configuration_fwdfit = configuration_fwdfit
      old_fixedpos = source_loc
    endelse
    
    
    ;;;;; If the user wants to select the location, loop on the number of locations to fix
    if keyword_set(select_location) then begin
      source_loc = fltarr(2,n_fwdfit_sources)

      ;; Ask if the user wants to use the previously selected boxes
      use_the_same = 0
      if old_nsources eq n_fwdfit_sources then begin
        print,''
        print,''
        print,'Previously selected locations (columns: x and y): '
        print, old_fixedpos
        print,''
        print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        read, use_the_same, prompt='Do you want to use the same locations? (1: yes, 0: no): '
      endif

      ;; If yes, use the same then
      if use_the_same eq 1 then begin
        source_loc = old_fixedpos 
        this_configuration_fwdfit = old_config_fwdfit

      endif else begin
        for ss=0,n_fwdfit_sources-1 do begin
          print,''
          print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          print,'Select on the CLEAN map the location of the source '+strtrim(ss,2)+': '
          xy, pos=tmp_pos

          source_loc[*,ss] = [tmp_pos[0], tmp_pos[1]]
        endfor
        old_fixedpos = source_loc
      endelse
    endif
    
    
    ;;;;; If the user wants to constrain the location, loop on the number of locations to constrain
    box_or_fix = fltarr(n_fwdfit_sources)    ; if 0: BOX, if 1: FIXED
    if keyword_set(select_box_location) then begin
      box_location = fltarr(2,2,n_fwdfit_sources)   ; [bl or tr, x or y, number of sources]

      ;; Ask if the user wants to use the previously selected locations
      use_the_same_box = 0
      if old_nsources eq n_fwdfit_sources then begin
        print,''
        print,''
        print,'Previously constrained locations (columns: x and y): '
        print, old_boxes
        print,''
        print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        read, use_the_same_box, prompt='Do you want to use the same boxes? (1: yes, 0: no): '
      endif

      ;; If yes, use the same then
      if use_the_same_box eq 1 then begin
        box_location = old_boxes
        this_configuration_fwdfit = old_config_fwdfit
        box_or_fix = all_box_or_fix[-n_fwdfit_sources:-1]

      endif else begin
        for ss=0,n_fwdfit_sources-1 do begin
          print,''
          print,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
          print,'Select on the CLEAN map the bottom-left corner of the box '+strtrim(ss,2)+': '
          print,'(if you select twice the same location, then it will be fixed!)'
          xy, pos=tmp_bl
          print,'Select on the CLEAN map the top-right corner of the box '+strtrim(ss,2)+': '
          xy, pos=tmp_tr

          box_location[0,*,ss] = [tmp_bl[0], tmp_bl[1]]
          box_location[1,*,ss] = [tmp_tr[0], tmp_tr[1]]
          threshold_fix = 2.
          if abs(box_location[0,0,ss] - box_location[1,0,ss]) lt threshold_fix and abs(box_location[0,1,ss] - box_location[1,1,ss]) lt threshold_fix then begin
            box_or_fix[ss] = 1
            print,''
            print,' --> The location of the source '+strtrim(ss,2)+' will be fixed!'
            print,''
            if keyword_set(ellipse_shape) then begin
              this_configuration_fwdfit[ss] = 'circle'
              print,''
              print,' --> The source '+strtrim(ss,2)+' is now a circular Gaussian'
              print,''
            endif
          endif
        endfor
        old_boxes = box_location
      endelse
    endif
    
    
    ;;;;; Initialize the structure needed by fwdfit
    srcin = VIS_FWDFIT_PSO_MULTIPLE_SRC_CREATE(vis, this_configuration_fwdfit)
    
    
    ;;;;; If the source location has been fixed or selected, then fix it
    ;; otherwise define the boundaries for the fwfit location (if specified)
    if keyword_set(source_loc) then begin
      counter_ellipses = 0
      for this_circ = 0,n_fwdfit_sources-1 do begin
        if this_configuration_fwdfit[this_circ] eq 'ellipse' then begin
          srcin.ellipse[counter_ellipses].param_opt.param_x = source_loc[0,this_circ]
          srcin.ellipse[counter_ellipses].param_opt.param_y = source_loc[1,this_circ]
          counter_ellipses += 1
        endif else begin
          srcin.circle[this_circ-counter_ellipses].param_opt.param_x = source_loc[0,this_circ]
          srcin.circle[this_circ-counter_ellipses].param_opt.param_y = source_loc[1,this_circ]
        endelse
        
      endfor
      all_fixedpos = [[all_fixedpos], [source_loc]]
      
    endif else if keyword_set(box_location) then begin
      counter_ellipses = 0
      for this_circ = 0,n_fwdfit_sources-1 do begin
        ;;; if (x,y) bottom-left and top-right are equal (or below ~2 arcsec), 
        ;;; then fix the location, if not, do boxes
        if abs(box_location[0,0,this_circ] - box_location[1,0,this_circ]) lt threshold_fix and abs(box_location[0,1,this_circ] - box_location[1,1,this_circ]) lt threshold_fix then begin
          srcin.circle[this_circ-counter_ellipses].param_opt.param_x = box_location[0,0,this_circ]
          srcin.circle[this_circ-counter_ellipses].param_opt.param_y = box_location[0,1,this_circ]
        endif else begin
          if this_configuration_fwdfit[this_circ] eq 'ellipse' then begin
            srcin.ellipse[counter_ellipses].lower_bound.l_b_x = box_location[0,0,this_circ] - mapcenter[0]
            srcin.ellipse[counter_ellipses].lower_bound.l_b_y = box_location[0,1,this_circ] - mapcenter[1]
            srcin.ellipse[counter_ellipses].upper_bound.u_b_x = box_location[1,0,this_circ] - mapcenter[0]
            srcin.ellipse[counter_ellipses].upper_bound.u_b_y = box_location[1,1,this_circ] - mapcenter[1]
            counter_ellipses += 1
          endif else begin
            srcin.circle[this_circ-counter_ellipses].lower_bound.l_b_x = box_location[0,0,this_circ] - mapcenter[0]
            srcin.circle[this_circ-counter_ellipses].lower_bound.l_b_y = box_location[0,1,this_circ] - mapcenter[1]
            srcin.circle[this_circ-counter_ellipses].upper_bound.u_b_x = box_location[1,0,this_circ] - mapcenter[0]
            srcin.circle[this_circ-counter_ellipses].upper_bound.u_b_y = box_location[1,1,this_circ] - mapcenter[1]
          endelse
          
        endelse
        
      endfor
      all_boxes = [[[all_boxes]], [[box_location]]]
    endif
    
    
    ;;;;; If the user pre-defined the source fwhm, then fix it
    if keyword_set(source_fwhm) then begin
      for this_circ = 0,n_fwdfit_sources-1 do begin
        if n_elements(source_fwhm) eq 1 then begin 
          srcin.circle[this_circ].param_opt.param_fwhm = source_fwhm 
        endif else begin
          srcin.circle[this_circ].param_opt.param_fwhm = source_fwhm[this_circ]
        endelse
        
      endfor
      all_fixedfwhm = [all_fixedfwhm, source_fwhm]
    endif
    
    
    ;;;;; If the user pre-defined the min fwhm of the source, then set it
    if keyword_set(min_fwhm) then begin
      counter_ellipses = 0
      for this_circ = 0,n_fwdfit_sources-1 do begin
        if n_elements(min_fwhm) eq 1 then begin
          if this_configuration_fwdfit[this_circ] eq 'ellipse' then begin
            srcin.ellipse[counter_ellipses].lower_bound.l_b_fwhm = min_fwhm
            counter_ellipses += 1
          endif else begin
            srcin.circle[this_circ-counter_ellipses].lower_bound.l_b_fwhm = min_fwhm
          endelse
           
        endif else begin
          if this_configuration_fwdfit[this_circ] eq 'ellipse' then begin
            srcin.ellipse[counter_ellipses].lower_bound.l_b_fwhm = min_fwhm[this_circ]
            counter_ellipses += 1
          endif else begin
            srcin.circle[this_circ-counter_ellipses].lower_bound.l_b_fwhm = min_fwhm[this_circ-counter_ellipses]
          endelse
          
        endelse
        
      endfor
      all_minfwhm = [all_minfwhm, min_fwhm]
    endif
    
    
    ;;;;; If the user pre-defined the max fwhm of the source, then set it
    if keyword_set(max_fwhm) then begin
      for this_circ = 0,n_fwdfit_sources-1 do begin
        if n_elements(max_fwhm) eq 1 then srcin.circle[this_circ].upper_bound.u_b_fwhm = max_fwhm else srcin.circle[this_circ].upper_bound.u_b_fwhm = max_fwhm[this_circ]
      endfor
      all_maxfwhm = [all_maxfwhm, max_fwhm]
    endif
    
    
    ;;;;; Create a copy for the sve-file
    srcin_fwdfit = srcin
    
    
    ;;;;; FWDFIT
    fwdfit_map = stx_vis_fwdfit_pso(this_configuration_fwdfit,vis,aux_data,srcin=srcin,imsize=imsize,$
      pixel=pixel,srcstr=out_param_fwdfit,fitsigmas=out_sigma_fwdfit,uncertainty=uncertainty,redchisq=chi2_fwdfit)


    ;;;;; Do visibility subtraction: observed visibilities - predicted visibilities by fwdfit
    vis_diff = vis
    vis_diff.obsvis = vis.obsvis - 1 * fwdfit_map.pred_vis.mapvis


    ;;;;; Create the BPROJ map out of the subtracted visibilities
    bp_diff_map = stx_bproj(vis_diff, imsize, pixel, aux_data)
    bp_diff_map.id = 'STX BPROJ (diff): '


    ;;;;; Plot visibility amplitude and phase
    stx_plot_fit_map, fwdfit_map


    ;;;;; Prepare some text for plots and save files
    txt_en_low = num2str(this_energy_range[0],format='(I10.0)')
    txt_en_high = num2str(this_energy_range[1],format='(I10.0)')
    txt_en_low2 = num2str(this_energy_range[0],format='(I3.3)')
    txt_en_high2 = num2str(this_energy_range[1],format='(I3.3)')
    text_erange = txt_en_low+'-'+txt_en_high
    text_erange_sav = txt_en_low2+'-'+txt_en_high2


    ;;;;; Plot the 5 maps
    window,11,xsize=1700,ysize=1200
    !p.multi=[0,3,2]
    loadct,5,/sil
    plot_map,clean_map[0],/limb,grid=grid,chars=ch_sz,title=title_clean
    if keyword_set(source_loc) then for i=0,n_elements(old_fixedpos[0,*])-1 do plots,old_fixedpos[0,i],old_fixedpos[1,i],/data,psym=1,syms=3,color=cgcolor('slate gray'),thick=3
    if keyword_set(box_location) then plot_boxes, old_boxes, box_or_fix
    loadct,5,/sil
    plot_map,memge_map,/limb,grid=grid,chars=ch_sz,title=title_mem
    if keyword_set(source_loc) then for i=0,n_elements(old_fixedpos[0,*])-1 do plots,old_fixedpos[0,i],old_fixedpos[1,i],/data,psym=1,syms=3,color=cgcolor('slate gray'),thick=3
    if keyword_set(box_location) then plot_boxes, old_boxes, box_or_fix
    loadct,5,/sil
    plot_map,fwdfit_map,/limb,grid=grid,chars=ch_sz
    if keyword_set(source_loc) then for i=0,n_elements(old_fixedpos[0,*])-1 do plots,old_fixedpos[0,i],old_fixedpos[1,i],/data,psym=1,syms=3,color=cgcolor('slate gray'),thick=3
    if keyword_set(box_location) then plot_boxes, old_boxes, box_or_fix
    loadct,5,/sil
    plot_map,bp_map,/limb,grid=grid,chars=ch_sz
    if keyword_set(source_loc) then for i=0,n_elements(old_fixedpos[0,*])-1 do plots,old_fixedpos[0,i],old_fixedpos[1,i],/data,psym=1,syms=3,color=cgcolor('slate gray'),thick=3
    if keyword_set(box_location) then plot_boxes, old_boxes, box_or_fix
    loadct,5,/sil
    plot_map,bp_diff_map,/limb,grid=grid,chars=ch_sz,dmax=max(bp_map.data)
    if keyword_set(source_loc) then for i=0,n_elements(old_fixedpos[0,*])-1 do plots,old_fixedpos[0,i],old_fixedpos[1,i],/data,psym=1,syms=3,color=cgcolor('slate gray'),thick=3
    if keyword_set(box_location) then plot_boxes, old_boxes, box_or_fix
    plot_map,bp_map,/nodata

    ;; Print some useful text
    xstart = 0.7
    ystart = 0.45
    dy = 0.06
    dyy = dy / 2

    if path_bkg_file eq '' then txt_counts = 'Tot counts: ' else txt_counts = 'Tot counts (bk sub): '
    txt_chi2 = num2str(chi2_fwdfit,format='(F10.2)')
    txt_nsrc = num2str(n_fwdfit_sources,format='(I10.0)')
    if keyword_set(source_loc) then txt_pos = 'fixed' else txt_pos = 'free'
    if keyword_set(box_location) then txt_pos = 'constrained'
    if total(box_or_fix) gt 0 then txt_pos = 'mix'
    if keyword_set(min_fwhm) or keyword_set(max_fwhm) then txt_fwhm = 'constrained' else txt_fwhm = 'free'
    if keyword_set(source_fwhm) then txt_fwhm = 'fixed'

    xyouts,xstart,ystart,'Energy range: '+text_erange+' keV',/norm,chars=ch_sz,chart=th_sz_txt
    xyouts,xstart,ystart-dy,txt_counts+num2str(ncounts,format='(I10.0)'),/norm,chars=ch_sz,chart=th_sz_txt
    xyouts,xstart,ystart-2*dy,'Tot fwdfit sources: '+txt_nsrc,/norm,chars=ch_sz,chart=th_sz_txt
    xyouts,xstart,ystart-2*dy-dyy,'Reduced chisq fwdfit: '+txt_chi2,/norm,chars=ch_sz,chart=th_sz_txt
    xyouts,xstart,ystart-3*dy-dyy,'fwdfit position:  '+txt_pos,/norm,chars=ch_sz,chart=th_sz_txt
    xyouts,xstart,ystart-3*dy-2*dyy,'fwdfit FWHM:    '+txt_fwhm,/norm,chars=ch_sz,chart=th_sz_txt
    xyouts,xstart,ystart-6*dy,vis_type,/norm,chars=ch_sz,chart=th_sz_txt

    if keyword_set(stop_here) then stop
    
    
    ;;;;; Store the data in a sav-file
    ; Check if the folder with the time range exists, othervise create it
    this_time_start_str = anytim(this_time_start)
    tmp_path = this_date+'_'+this_time_start+'-'+this_time_end+'_'+text_ref_time+suffix_folder
    filename_final_folder = path_sav_folder + tmp_path
    ;stop
    if not file_test(filename_final_folder,/dir) then file_mkdir, filename_final_folder
    
    ; Names of the save and png files
    filename_sav = folder_delimiter + 'stix-imaging-spectroscopy_'+text_erange_sav+'-keV_'+this_date+'_'+this_time_start+'-'+this_time_end+'_'+text_ref_time+'.sav'
    filename_png = folder_delimiter + 'stix-imaging-spectroscopy_'+text_erange_sav+'-keV_'+this_date+'_'+this_time_start+'-'+this_time_end+'_'+text_ref_time+'.png'

    save,filename=filename_final_folder+filename_sav,vis,clean_map,memge_map,fwdfit_map,$
      time_range_earth,time_range_so,this_configuration_fwdfit,out_param_fwdfit,$
      out_sigma_fwdfit,aux_data,this_energy_range,subc_index,beam_width,time_shift,distance,$
      chi2_fwdfit,text_ref_time,earth_ut,n_fwdfit_sources,srcin_fwdfit

    write_png,filename_final_folder+filename_png,tvrd(/true)

    !p.multi=0
    
    path_new_folder = filename_final_folder
    
    
    ;;;;; Store some values for the next loop-step
    old_nsources = n_fwdfit_sources
    old_config_fwdfit = this_configuration_fwdfit
    all_box_or_fix = [all_box_or_fix, box_or_fix]
    
    
  endfor
  
  
  ;; Plot CLEAN again to show all selected locations
  if keyword_set(source_loc) then begin
    loadct,5,/silent
    !p.multi = 0
    window,12,xsize=1300,ysize=1300;,xpos=2300
    plot_map,clean_map[4],/limb,grid=grid,chars=ch_sz
    plot_map,clean_map[4],/over,/perc,levels=[30,50],color=cgcolor('cyan'),thick=2
    ;plot_map,clean_map[4],/over,levels=[max(clean_map[2].data)],color=cgcolor('dark green'),thick=2
    ;xyouts,0.21,0.75,'Peak residual',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('dark green')
    xyouts,0.21,0.8,'[30,50]%',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('cyan')
    for i=0,n_elements(all_fixedpos[0,*])-1 do plots,all_fixedpos[0,i],all_fixedpos[1,i],/data,psym=1,syms=5,color=cgcolor('Slate Gray'),thick=3
    xyouts,0.21,0.77,'All selected locations',/normal,charsize=ch_sz,chart=2.2,color=cgcolor('Slate Gray')
  endif
  
  
  
end   ; End of stx_imaging_spectroscopy