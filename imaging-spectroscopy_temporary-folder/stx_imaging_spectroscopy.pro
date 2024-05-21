;+
;
; PROJECT:
;   STIX
;
; NAME:
;   stx_imaging_spectroscopy
;
; PURPOSE:
;   Perform imaging spectroscopy with STIX. This procedures calls stx_imaging spectroscopy_photon
;   or stx_imaging_spectroscopy_electron to produce the maps and files. Depending on the user, 
;   the photon visibilities (observed or regularized) or the electron visibilities can be used.
;   By default:
;       - the regularized photon maps are created
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
;                          visibilities. If the keyword "/observed_vis" is set, then this is ignored.
;
; OPTIONAL OUTPUTS:
;   path_new_folder      : string containing the path of the newly created folder where the
;                          results are stored
;
; OPTIONAL INPUTS AND KEYWORDS:
;   electron_maps        : if set, the electron maps are created; else by default the regularized photon maps are created
;   
;   observed_vis         : if set, then the observed visibilities are used. In this case, the
;                          energy_max_inversion is ignored.
;
;   energy_min_inversion : minimum energy to use for the inversion to calculate the regularized
;                          visibilities. If the keyword "/observed_vis" is set, then this is ignored.
;                          Default is 9 keV (we cannot go lower with the regularized visibilities).
;
;   mapcenter            : array with x,y coordinates for the mapcenter (if known)
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
;                                     [[[bl_x1, tr_x1],[bl_y1, tr_y1]],          ; bottom-left and top-right(arcsec) box source 1
;                                      [[bl_x2, tr_x2],[bl_y2, tr_y2]]]          ; bottom-left and top-right(arcsec) box source 2
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
;   no_mem_ge            : if set, it does not produce any MEM_GE maps
;
;   no_clean             : if set, it does not produce any CLEAN maps
;
;   no_fwdfit            : if set, it does not produce any FWDFIT maps
;
; HISTORY:
;   February 2024, Stiefel M. Z. (FHNW & ETHZ), initial release
;
; CONTACT:
;   muriel.stiefel@fhnw.ch
;
;-

pro stx_imaging_spectroscopy, path_sci_file, path_bkg_file, aux_fits_file, time_range, energy_max_inversion,$     ; inputs
  ;; --- Optional inputs and keywords
  electron_maps = electron_maps, $
  observed_vis = observed_vis, $
  energy_min_inversion = energy_min_inversion, $
  configuration_fwdfit = configuration_fwdfit, $
  mapcenter = mapcenter, $
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
  no_mem_ge = no_mem_ge, $
  no_clean = no_clean, $
  no_fwdfit = no_fwdfit, $
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
  
  if keyword_set(electron_maps) and keyword_set(observed_vis) then begin
    print,'Keyword electron_maps and observed_vis can not be used at the same time'
    message,'---> Please, review your keywords!'
  endif
  
  if keyword_set(electron_maps) then begin
    ;Create the electron maps
    
    stx_imaging_spectroscopy_electron, path_sci_file, path_bkg_file, aux_fits_file, time_range, energy_max_inversion,$     ; inputs
      ;; --- Optional inputs and keywords
      energy_min_inversion = energy_min_inversion, $
      configuration_fwdfit = configuration_fwdfit, $
      mapcenter = mapcenter, $
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
      no_mem_ge = no_mem_ge, $
      no_clean = no_clean, $
      no_fwdfit = no_fwdfit, $
      ;; --- Standard imaging inputs
      subc_index = subc_index, $
      pixels = pixels, $
      imsize = imsize, $
      ;pixel = pixel, $
      niter = niter, $
      gain = gain, $
      nmap = nmap, $
      beam_width = beam_width, $
      ;; --- Optional output
      path_new_folder = path_new_folder
      
  endif else begin
    ;Create the photon maps
    
    stx_imaging_spectroscopy_photon, path_sci_file, path_bkg_file, aux_fits_file, time_range, energy_max_inversion,$     ; inputs
      ;; --- Optional inputs and keywords
      observed_vis = observed_vis, $
      energy_min_inversion = energy_min_inversion, $
      configuration_fwdfit = configuration_fwdfit, $
      mapcenter = mapcenter, $
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
      no_mem_ge = no_mem_ge, $
      no_clean = no_clean, $
      no_fwdfit = no_fwdfit, $
      ;; --- Standard imaging inputs
      subc_index = subc_index, $
      pixels = pixels, $
      imsize = imsize, $
      ;pixel = pixel, $
      niter = niter, $
      gain = gain, $
      nmap = nmap, $
      beam_width = beam_width, $
      ;; --- Optional output
      path_new_folder = path_new_folder
    
  endelse
  
  
end
