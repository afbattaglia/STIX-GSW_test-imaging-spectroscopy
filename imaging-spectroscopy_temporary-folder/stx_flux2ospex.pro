;+
;
; PROJECT:
;   STIX
;
; NAME:
;   stx_flux2ospex
;
; PURPOSE:
;   This procedure takes the flux structure obtained from stx_plot_imaging_spectra
;   and puts it in OSPEX, which allows to fit the imaging spectrum.
;
; CALLING SEQUENCE:
;   stx_flux2ospex, flux_str
;
; INPUTS:
;   flux_str     : the flux structure obtained from stx_plot_imaging_spectra
;
; OPTIONAL INPUTS AND KEYWORDS:
;   ind_sources  : the indices of the sources to sum up and pass to OSPEX.
;                  Default: sum up all sources
;
;   stop_here    : if set, the procedure stops before exiting (for debugging)
;
; OPTIONAL OUTPUTS:
;   ospex_obj    : the OSPEX object containing the imaging spectrum
;
; HISTORY:
;   July 2023, Battaglia A. F. (FHNW & ETHZ), initial release
;
; CONTACT:
;   andrea.battaglia@fhnw.ch
;
;-

pro stx_flux2ospex, flux_str, $ ; input
  ;; --- Optional inputs and keywords
  ind_sources = ind_sources, $
  stop_here = stop_here, $
  ;; --- Optional output
  ospex_obj = ospex_obj


  ;;;;; Default parameters
  default, ind_sources, indgen(n_elements(flux_str.FF_SOURCES[*,0]))

  ;;;;; Extract the fluxes
  all_fluxes = flux_str.ff_sources
  all_errors = flux_str.ff_sources_err

  ;;;;; Extract the SRM structure
  this_srm = flux_str.srm

  ;;;;; Time range and duration
  time_range = flux_str.time_range
  duration = anytim(time_range[1]) - anytim(time_range[0])

  ;;;;; Energy axis
  used_ct_energy_edges = [flux_str.e_min, flux_str.e_max]
  ct_energy_edges = flux_str.srm.edges_out[0:29]
  ;ph_energy_edges = this_srm.ph_edges ;flux_str.srm.edges_in;[0:-1]
  ph_energy_edges = flux_str.srm.edges_in;[0:-1]
  used_e_bin_width = used_ct_energy_edges[1,*] - used_ct_energy_edges[0,*]
  e_bin_width = ct_energy_edges[1:-1] - ct_energy_edges[0:-2]
  ct_energy_edges = transpose([[ct_energy_edges], [ct_energy_edges + e_bin_width, 30.]])

  ;;;;; Detector area
  det_area = flux_str.det_area

  ;;;;; Sum up the sources, if needed
  if n_elements(ind_sources) gt 1 then begin
    summed_sources = total(all_fluxes[ind_sources,*],1) * used_e_bin_width * det_area; * duration
    summed_sources_err = sqrt(total(all_errors[ind_sources,*]^2,1)) * used_e_bin_width * det_area; * duration

  endif else begin
    summed_sources = transpose(all_fluxes[ind_sources,*]) * used_e_bin_width * det_area; * duration
    summed_sources_err = transpose(all_errors[ind_sources,*]) * used_e_bin_width * det_area; * duration

  endelse

  ;;;;; Now put the input spectrum in the correct format, wrt the starting energy
  input_spectrum     = fltarr(30) * !values.f_nan
  input_spectrum_err = fltarr(30) * !values.f_nan
  start_en = flux_str.e_min[0]
  ind_start_en = where(abs(start_en - this_srm.edges_out) eq min(abs(start_en - this_srm.edges_out)))
  ind_end_en = ind_start_en + n_elements(summed_sources)
  input_spectrum[ind_start_en:ind_end_en-1] = summed_sources
  input_spectrum_err[ind_start_en:ind_end_en-1] = summed_sources_err


  ;;;;; Adjust the size of the input spectrum
  ;  n_en = n_elements(input_spectrum)
  ;  for i=n_en,29 do input_spectrum = [input_spectrum, !values.f_nan]

  ;  input_spectrum = summed_sources
  ;  input_spectrum_err = summed_sources_err
 
 ;cut out unused energies from OSEPX input
  valid = where(finite(input_spectrum) eq 1)

  input_spectrum = input_spectrum[valid]
  input_spectrum_err = input_spectrum_err[valid]
  ct_energy_edges = ct_energy_edges[*,valid]

  this_srm = rep_tag_value(this_srm, (this_srm.drm)[valid,all],'drm')
  this_srm = rep_tag_value(this_srm, (this_srm.edges_out)[valid],'edges_out')

  ;;;;; Open OSPEX
  ospex_obj = ospex(/no)

  ;;;;; Set the parameters
  ospex_obj -> set, spex_data_source = 'SPEX_USER_DATA'
  ospex_obj -> set, spectrum = input_spectrum, $
    spex_ct_edges = ct_energy_edges, $
    spex_ut_edges = anytim(time_range), $
    ;livetime = livetime, $
    errors = input_spectrum_err
  ospex_obj -> set, spex_drm_ct_edges = ct_energy_edges
  ospex_obj -> set, spex_drm_ph_edges = ph_energy_edges
  ospex_obj -> set, spex_respinfo = this_srm
  ospex_obj -> set, spex_area = this_srm.area
  ospex_obj -> set, spex_detectors = 'STIX (fwdfit)'
  ospex_obj -> gui

  if keyword_set(stop_here) then stop

end ; End of the script
