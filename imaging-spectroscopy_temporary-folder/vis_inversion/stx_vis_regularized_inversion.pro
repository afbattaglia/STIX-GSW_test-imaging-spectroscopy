;+
;
; NAME:
;   stx_vis_regularized_inversion.pro
;
; PURPOSE:
;     for a given (real or imaginary) count visibility spectrum, performs the regularized inversion as follows:
;     1) check the consistency of the input vectors
;     2) rescale the input count visibility spectrum
;     3) invert the count visibility spectrum to provide regularized electron visibility spectrum according to Eq.3
;        in Piana et al. The Astrophysical Journal, 665:846-855, 2007 August 10
;     4) repeat the inversion several times in order to provide a confidence strip around the regularized electron visibility spectrum
;     5) save the inverted spectra in a structure
;
; CALLING SEQUENCE:
;     stx_vis_regularized_inversion, l, part, count_vis_spectra, datastruct, drmini, electron_bin=electron_bin, confidencestrip=confidencestrip
;
; CALLED BY:
;   - stx_visibilities_inversion.pro
;
; CALLS TO:
;   - stx_vis_reg_ge_sampling.pro
;   - stx_vis_reg_ge_cross3bn.pro
;   - stx_vis_reg_ge_regularization.pro
;   - stx_vis_reg_ge_confidence.pro
;
; INPUTS:
;   l:                 label of the spectrum (l=0, .., Total number of spectra to invert -1)
;   part:              'Real' or 'Imaginary'
;   count_vis_spectra: the 7-tag structure, output of the build_count_vis_spectra.pro code, containing
;                        the information on the count visibility spectra
;   drmini:            the drm computed by stx_build_drm.pro (drm is assumed to be the same for all detectors)
;   dist_solo_sun:     SolO - Sun distance
;
; KEYWORDS:
;   confidencestrip    - Number of solution realizations to determine. 
;                        (Range [1,50] - Default = 10)
;
; OUTPUTS:
;   DATASTRUCT: a structure storing 2 arrays
;               - Tag 1 :     STRIP       DOUBLE    Array[NEE, N+2] where NEE = number of electron energy points
;                                                                         N   = number of solution realizations 
;                                                                               (defined by confidencestrip value)
;              
;                     The array contains the electron energies and each realization of the regularized solution
;                     in the following format:              
;                         datastruct.strip[*,0]     = a vector of NEE elements containing the electron energies
;                         datastruct.strip[*,1]     = a vector of NEE elements containing the original (unperturbed) regularized
;                                                     solution at the specified energies
;                         datastruct.strip[*,2:N+1] = array containing N different regularized electron visibilities, at the
;                                                     specified energies, produced by randomly perturbing the input data
;              
;               - Tag 2 : RESIDUALS       DOUBLE    Array[NPH, 6] where NPH = number of count energy points
;              
;                     The data is in six columns, with each row arranged as follows:;              
;                         datastruct.residuals[*,0] = Count energy used
;                         datastruct.residuals[*,1] = Count visibility values (input data)
;                         datastruct.residuals[*,2] = Uncertainty in count visibility values (input data)
;                         datastruct.residuals[*,3] = Count visibilities corresponding to the recovered electron visibility array
;                         datastruct.residuals[*,4] = Residual count visibilities (actual - recovered, i.e., column 2 - column 4)
;                         datastruct.residuals[*,5] = Cumulative residual, defined as
;                                                        C_j = (1/j) sum_[i=1]^j res_i
;
;+
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro stx_vis_regularized_inversion, l, part, opt, count_vis_spectra, datastruct, drmini, dist_solo_sun, $
                                    confidencestrip = confidencestrip ; Number of solution realizations to determine. 


  
  DEFAULT, el_energy_max_factor,  2
  DEFAULT, crosssection,          'Cross3BN'
  DEFAULT, Z,                     1.2
  
  photon_bin_position = 0.5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (part eq 'Real') then vis_array = count_vis_spectra.RealPart[*,l] else vis_array = count_vis_spectra.ImagPart[*,l]
  photon_energy_low_array  = count_vis_spectra.energyl
  photon_energy_high_array = count_vis_spectra.energyh
  
  photon_energy_edges = [photon_energy_low_array, photon_energy_high_array[-1]]
  
  bin      = photon_energy_high_array-photon_energy_low_array
  bin_unif = bin[where(photon_energy_low_array le 16.)]
  
  epsini = photon_energy_low_array + bin*photon_bin_position
  eemax = el_energy_max_factor*epsini[-1]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;; Some cheking;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;  Testing vector sizes and contents

  n = N_elements(count_vis_spectra.energyl)
  
  if ( (n ne N_elements(count_vis_spectra.energyh)) or $
       (n ne N_elements(vis_array)) or $
       (n ne N_elements(count_vis_spectra.error[*,l])) ) then message, 'INPUT VECTOR SIZES NOT COMPATIBLE'

  if ( (total(count_vis_spectra.energyl) eq 0) or $
       (total(count_vis_spectra.energyh) eq 0) or $
       (total(vis_array) eq 0) or $
       (total(count_vis_spectra.error[*,l]) eq 0)) then message, 'INPUT VECTOR CONTENTS ARE ZEROS'

  if ( where(vis_array eq 'Inf') ne -1 )  then message, 'Inf values in the visibility array'
  
  if ( where(count_vis_spectra.error[*,l] eq 'Inf') ne -1 )  then message, 'Inf values in the visibility array'


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;; The code automatically fits the photon spectrum with a single, broken
  ;;;;;;;;;;; power law or sum of power laws. Then it rescales the spectrum.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  stx_vis_reg_ge_rescaling, count_vis_spectra.energyl, count_vis_spectra.energyh, vis_array, $
                            count_vis_spectra.error[*,l], drmini, el_energy_max_factor, dist_solo_sun, $
                            nph,eps,gs,err,wei,fitph,fitel,drm

  if( nph le 3) then begin
    datastruct={strip:0,residuals:0}
    return
  endif
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;,;
  ;;;;;;;;;;; Electron energies sampling
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  e_stx_bin = [indgen(13)+4., indgen(3.)*2.+18., indgen(2.)*3.+25., indgen(3.)*4.+32., indgen(2.)*5.+45., indgen(2.)*7.+56., 70., 76., 84., 100., 120., 150.]
  delta_stx = [fltarr(13)+1., fltarr(3.)+2., fltarr(2.)+3., fltarr(3.)+4., fltarr(2.)+5., indgen(2.)+6., 7., 6., 8.,16., 20., 30.]

  ee_sampl = photon_energy_edges
  delta_e  = bin
  index_min = min(where(e_stx_bin gt photon_energy_edges[-1]))
  index_max = min(where(e_stx_bin gt eemax))
  ee_sampl = [ee_sampl, e_stx_bin[index_min:index_max]]
  delta_e = [delta_e, delta_stx[index_min:index_max+1]]
  
  ee  = ee_sampl + delta_e*photon_bin_position
  nee = n_elements(ee)
  
  ;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;; Cross-section computation
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  eebrem=0
  stx_vis_reg_ge_cross3bn,eebrem,nph,nee,eps,ee,drm,fitph,fitel,Z,dist_solo_sun, w,u,sf

  ;;;;;;;;; Rescaling of the singular vectors and
  ;;;;;;;;; values according to the weights
  for k=0,nph-1 do begin
    w[k]=sqrt(wei[k])*w[k]
    u[*,k]= u[*,k]/sqrt(wei[k])
  end
  ;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;; Regularized solution computation
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;rescaling
  gs=gs/fitph
  err=err/fitph
  ;;;;;;;;;;;;;;;;;;;;

  vis_reg_ge_regularization,nph,nee,wei,eps,gs,err,u,w,sf,fitph,fitel,$
    ee,opt,gopt,regsol,residual_array, cumulativeresidual_array

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;; Confidence strip computation
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  vis_reg_ge_confidence,nph,nee,ee,regsol,$
    err,gs,eps,wei,confidencestrip,u,w,opt,sf,drm,fitph,fitel,Z,datastruct

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;  Fill the fourth tag of datastruct structure and save it ;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  datastruct.residuals[*,0]=eps[*]
  datastruct.residuals[*,1]=gs[*]
  datastruct.residuals[*,2]=err[*]
  datastruct.residuals[*,3]=gopt[*]
  datastruct.residuals[*,4]=residual_array[*]
  datastruct.residuals[*,5]=cumulativeresidual_array[*]

end