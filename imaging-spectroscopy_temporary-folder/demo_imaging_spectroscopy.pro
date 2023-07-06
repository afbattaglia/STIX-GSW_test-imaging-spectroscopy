; **************************************************************
; Demo for the imaging spectroscopy software.
;
; ---------
; Please, read the following document for more information:
;     User-Guide-to-STIX-Imaging-Spectroscopy_Status-July2023.pdf
; ---------
; **************************************************************


; ******************* Example with the 2022 Sep 30th flare ******************************

;;;;; Path where to store the save and png files that will be generated
path_sav_folder = '/home/afbattaglia/Software/idl_lib/STIX-GSW_test-imaging-spectroscopy/'

;;;;; Path to the science, bkg and auxiliary files
;; It only works by giving the ABSOLUTE PATH!!!
;; You can find the data in the following folder:
;;          imaging-spectroscopy_temporary-folder/data4demo/
path_folder_data = '/home/afbattaglia/Software/idl_lib/STIX-GSW_test-imaging-spectroscopy/STIX-GSW_test-imaging-spectroscopy/imaging-spectroscopy_temporary-folder/data4demo/'
this_path_sci_file = path_folder_data + 'cpd_4s/solo_L1_stix-sci-xray-cpd_20220930T160153-20220930T162949_V01_2209303250-64438.fits'
this_path_bkg_file = path_folder_data + 'bkg/solo_L1_stix-sci-xray-cpd_20221002T092421-20221002T101741_V01_2210022580-50759.fits'
this_aux_fits_file = path_folder_data + 'aux/solo_L2_stix-aux-ephemeris_20220930_V01.fits'

;;;;; Time range
; It can be given in Solar Orbiter UT or Earth UT
; If you want to use Earth UT times, then you MUST set the keyword earth_ut
this_time_range_so = ['30-Sep-2022 16:15:35', '30-Sep-2022 16:16:59']

;;;;; Energy range
;; By setting this, the script will loop all native energy bins within the specified range
energy_range = [9, 32]
;; By setting this, you can bin in energy
;energy_low = [32,40,50]
;energy_high = [40,50,70]

;;;;; Maximum energy to use for the inversion to calculate the regularized visibilities
energy_max_inversion = 63

;;;;; Suffix to append at the end of the newly created folder
;suffix_folder = '_test-new-software_std-visibilities'
suffix_folder = '_test-new-software_reg-visibilities'

;;;;; Set the minimum/maximum size of the source FWHM
min_fwhm = 14.6  ; corresponding to the resolution of sc3
max_fwhm = 178.6  ; corresponding to the resolution of sc10

;;;;; If to calculate the uncertainty on the fwdfit parameters
;; Default: uncertainty = 1 (i.e., calculate the uncertainties)
uncertainty = 1

; ******************************************************************************************







;; STEP 1: run the imaging spectroscopy tool
stx_imaging_spectroscopy, $
  ;; --- Necessary inputs
  this_path_sci_file, $
  this_path_bkg_file, $
  this_aux_fits_file, $
  this_time_range_so, $
  energy_max_inversion, $
  ;; --- Optional inputs and keywords
  ;/select_loc, $
  ;/observed_vis, $
  path_sav_folder = path_sav_folder, $
  ;/stop_here, $
  /earth_ut, $
  ;source_fwhm = source_fwhm
  min_fwhm = min_fwhm, $
  ;max_fwhm = max_fwhm, $
  ;/ellipse_shape, $
  /select_box, $
  ;energy_low = energy_low, $
  ;energy_high = energy_high, $
  ;configuration_fwdfit = ['circle','circle'], $
  ;source_loc = [[1188.896, 404.447], [1188.896, 404.447]], $
  uncertainty = uncertainty, $
  energy_range = energy_range, $
  suffix_folder = suffix_folder, $
  ;; --- Output
  path_new_folder = path_new_folder


;; STEP 2: get the fluxes of all sources
stx_plot_imaging_spectra, $
  ;; -- Necessary inputs
  path_new_folder, $
  this_path_sci_file, $
  this_path_bkg_file, $
  ;; -- Optional output
  flux_str = flux_str


;; STEP 3: fit the imaging spectra in OSPEX
stx_flux2ospex, $
  ;; -- Necessary input
  flux_str, $
  ;; -- Optional input
  ind_sources = [1,2], $
  ;; -- Optional output
  ospex_obj = ospex_obj


end ; End of the script