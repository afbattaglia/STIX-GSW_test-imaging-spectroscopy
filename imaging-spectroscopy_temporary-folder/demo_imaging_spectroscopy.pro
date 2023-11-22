; **************************************************************
; Demo for the imaging spectroscopy software.
;
; ---------
; Please, read the following document for more information:
;    User-Guide-to-STIX-Imaging-Spectroscopy_Status-September2023.pdf
; ---------
; 
; Andrea Francesco Battaglia
; Last modification date: 27-Sep-2023
; **************************************************************


; ********** Example with the 2023 May 16 17:20 flare **********

;;;;; Path where to store the save and png files that will be generated
path_sav_folder = 'path_to_folder'

;;;;; Path to the science, bkg and auxiliary files
;; It only works by giving the ABSOLUTE PATH!!!
;; You can find the data in the following folder:
;;          imaging-spectroscopy_temporary-folder/data4demo/
path_data4demo = 'path_to_the_folder/data4demo/'
this_path_sci_file = path_data4demo + 'cpd/solo_L1_stix-sci-xray-cpd_20230516T170425-20230516T173914_V01_2305162947-52420.fits'
this_path_bkg_file = path_data4demo + 'bkg/solo_L1_stix-sci-xray-cpd_20230517T061006-20230517T065706_V01_2305176702-52442.fits'
this_aux_fits_file = path_data4demo + 'aux/solo_L2_stix-aux-ephemeris_20230516_V01.fits'

;;;;; Time range
; It can be given in Solar Orbiter UT or Earth UT
; If you want to use Earth UT times, then you MUST set the keyword earth_ut
time_range = ['16-May-2023 17:20:30', '16-May-2023 17:21:50']

;;;;; Energy range
;; By setting this, the script will loop all native energy bins within the specified range
energy_range = [18, 50]
;; By setting this, you can bin in energy
;energy_low = [32,40,50]
;energy_high = [40,50,70]

;;;;; Maximum energy to use for the inversion to calculate the regularized visibilities
energy_max_inversion = 100

;;;;; Suffix to append at the end of the newly created folder
;; This allows to run the imaging spectroscopy software with different settings 
;; by creating different folders, instead of overwrite them
suffix_folder = '_demo'

;;;;; Set the minimum/maximum size of the source FWHM
min_fwhm = 14.6  ; corresponding to the resolution of sc3
max_fwhm = 178.6  ; corresponding to the resolution of sc10

;;;;; If to calculate the uncertainty on the fwdfit parameters
;; Default: uncertainty = 1 (i.e., calculate the uncertainties)
uncertainty = 1

;;;;; If the configuration is known a priori, set this, otherwise
;; undo the configuration_fwdfit keyword
;configuration_fwdfit = ['circle','circle']

;;;;; If you want to set the location of the boxes
;   box_location         : instead of selecting the boxes on the screen, one can pass it and then it will
;                          be used for all the energy bins. One has to give the coordinates of the bottom left
;                          and top-right corners of the box as
;                                     [[bl_x, bl_y],          bottom-left (arcsec)
;                                      [tr_x, tr_y]]          top-right (arcsec)
;                          If one needs to give more boxes, then it has to be of the following format: [2,2,nb_sources]
;                          For example:
;                                     [[[bl_x1, tr_x1],[bl_y1, tr_y1]],          ; bottom-left and top-right(arcsec) box source 1
;                                      [[bl_x2, tr_x2],[bl_y2, tr_y2]]]          ; bottom-left and top-right(arcsec) box source 2
;box_location = [[[-811.740,-767.401],[215.066,265.086]],$
;  [[-765.989,-724.191],[192.787,237.691]]]

; ******************************************************************************************





;; STEP 1: run the imaging spectroscopy tool
stx_imaging_spectroscopy, $
  ;; --- Necessary inputs
  this_path_sci_file, $
  this_path_bkg_file, $
  this_aux_fits_file, $
  time_range, $
  energy_max_inversion, $
  ;; --- Optional inputs and keywords
  ;/select_loc, $
  ;/observed_vis, $
  path_sav_folder = path_sav_folder, $
  ;/stop_here, $
  ;/earth_ut, $
  ;source_fwhm = source_fwhm
  min_fwhm = min_fwhm, $
  max_fwhm = max_fwhm, $
  ;/ellipse_shape, $
  /select_box, $
  ;energy_low = energy_low, $
  ;energy_high = energy_high, $
  ;configuration_fwdfit = configuration_fwdfit, $
  ;source_loc = source_loc, $
  ;box_location = box_location, $
  uncertainty = uncertainty, $
  energy_range = energy_range, $
  suffix_folder = suffix_folder, $
  ;/no_mem, $
  ;/no_clean, $
  ;/no_fwdfit, $
  ;; --- Output
  path_new_folder = path_new_folder


;; OPTIONAL STEP: plot the fwdfit and spatially integrated spectra
stx_plot_imaging_spectra, $
  ;; -- Necessary inputs
  path_new_folder, $     ; this is an output of stx_imaging_spectroscopy
  this_path_sci_file, $
  this_path_bkg_file


;; STEP 2: input the fluxes into OSPEX
algo = 'fwdfit'
stx_img_spectra_sav2ospex, path_new_folder, algo



end ; End of the script