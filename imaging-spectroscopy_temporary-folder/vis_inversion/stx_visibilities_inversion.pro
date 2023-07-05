;+
;
; NAME:
;   stx_visibilities_inversion.pro
;
; PURPOSE:
;     this is the main routine of the Visibility Inversion Software (VIS) implementing the
;     method described in Piana et al. The Astrophysical Journal, 665:846-855, 2007 to compute
;     electron visibilities from STIX count visibilities.
;     In particular this routine:
;     1) stx_build_count_vis_spectra, for each (u,v) pair, as a function of the count energy
;     2) invert the count visibility spectra to provide regularized electron visibility spectra
;     3) save the inverted spectra in arrays of visibility structures in the same format of the
;        standard arrays of count visibility structures.
;
; CALLING SEQUENCE:
;       stx_visibilities_inversion, vis, dist_solo_sun, attenuator = attenuator, $
;                                    reg_el_vis, orig_ph_vis, reg_ph_vis
;
;
; CALLS TO:
;
;   - stx_build_count_vis_spectra.pro
;   - stx_include_drm_effects.pro
;   - stx_vis_regularized_inversion.pro
;   - stx_el_spectra_2_el_vis.pro
;
; INPUTS:
;   vis:           array of STIX count visibility structures (counts cm^-2 s^-1 keV^-1).
;                  Each element of the array must be a STIX visibility structure corresponding to a specific energy range              
;   dist_solo_sun: SolO - Sun distance 
;   
; KEYWORDS:
;   attenuator = 0 - no attenuator, 1 - include the attenuator
;   
; OUTPUTS
;   reg_el_vis:  array of electron visibility structures for the basic regularized solution (electrons cm^-2 s^-1 keV^-1)
;   orig_ph_vis: array of the sub-set of the count visibility structures used as input for the inversion (counts cm^-2 s^-1 keV^-1)
;   reg_ph_vis:  array of regularized count visibility structures corresponding to the recovered electron visbilities (counts cm^-2 s^-1 keV^-1)
;
;+

pro stx_visibilities_inversion, vis, dist_solo_sun, attenuator, confidencestrip = confidencestrip, $
                                reg_el_vis, orig_ph_vis, reg_ph_vis
                                
                                
  default, attenuator, 0 ; no attenuator
  default, confidencestrip, 10
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;; Check for the detectors involved
  detOK = intarr(32)
  FOR isc = 1, 32 DO BEGIN
    isthere = where(vis.isc eq isc, cont)
    if (cont GT 0) then detOK[isc-1] = 1
  END
  
  ct_edges_lower  = vis.energy_range[0]
  ct_edges_upper  = vis.energy_range[1]
  ct_edges_lower1 = ct_edges_lower[UNIQ(ct_edges_lower, SORT(ct_edges_lower))]
  ct_edges_upper1 = ct_edges_upper[UNIQ(ct_edges_upper, SORT(ct_edges_upper))]
  ct_edges = [ct_edges_lower1, ct_edges_upper1[-1]]
  
  deltaeps = ct_edges_upper1 - ct_edges_lower1
  neps     = n_elements(ct_edges_lower1)  
  
  if ((ct_edges_upper1[0]+ct_edges_lower1[0])/2 lt 9) then begin
    message, 'Check for low energy values: Please the first energy bin must be higher than 9 keV'
  endif

  ;;;;;;;; Build count visibility spectra
  stx_build_count_vis_spectra, detOK, vis, count_vis_spectra

  ;;;;;;;; DRM computation for the involved detectors (drm is assumed to be the same for all detectors)
  e_stx_bin = [indgen(13)+4., indgen(3.)*2.+18., indgen(2.)*3.+25., indgen(3.)*4.+32., indgen(2.)*5.+45., indgen(2.)*7.+56., 70., 76., 84., 100., 120., 150.]
  ct_edges2 = [ct_edges, e_stx_bin[min(where(e_stx_bin gt ct_edges[-1]))]]
  stx_drm = stx_build_drm(ct_edges2, attenuator=attenuator)
  srm = stx_drm.SMATRIX ; smatrix is DRM (counts/keV/photons - Array[32, 32]). 
  
  ;;;;;;;; Inversion
  Nspectra   = total(detOK, /integer)
  spectraptr = ptrarr(Nspectra)

  confidencestrip = 10
  for l = 0, Nspectra-1 do begin
      print, 'Inversion of spectra n.', l
      ;;;; spectrum 'l' comes from detector 'count_vis_spectra.det[l]'
      drmini = srm;
      
      ;;;;;;;; Real part inversion
      part = 'Real'
      stx_vis_regularized_inversion, l, part, opt, count_vis_spectra, outdataRe, drmini, dist_solo_sun, $
                                     confidencestrip = confidencestrip

      nphRe = N_elements(outdataRe.residuals[*,0])
      neeRe = N_elements(outdataRe.strip[*,0])

      ;;;;;;;; Imaginary part inversion
      part='Imaginary'
      stx_vis_regularized_inversion, l, part, opt, count_vis_spectra, outdataIm, drmini, dist_solo_sun, $
                                     confidencestrip = confidencestrip


      nphIm = N_elements(outdataIm.residuals[*,0])
      neeIm = N_elements(outdataIm.strip[*,0])

      ;;;;;;;; Save inverted spectra

      outdata = { stripRe : dblarr(neeRe,confidencestrip+2), $
                  residualsRe : dblarr(nphRe,6), $
                  stripIm:dblarr(neeIm,confidencestrip+2),$
                  residualsIm:dblarr(nphIm,6) }
                  
      outdata.stripRe     = outdataRe.strip
      outdata.stripIm     = outdataIm.strip
      outdata.residualsRe = outdataRe.residuals
      outdata.residualsIm = outdataIm.residuals
      
      print, ' '
      print, nphRe
      print, neeRe
      print, ' '
      print, nphIm
      print, neeIm

      spectraptr[l] = ptr_new(outdata)
      
  end

  print, ''

  stx_vis_el_spectra_2_el_vis, vis, count_vis_spectra, spectraptr, reg_el_vis, orig_ph_vis, reg_ph_vis

end
