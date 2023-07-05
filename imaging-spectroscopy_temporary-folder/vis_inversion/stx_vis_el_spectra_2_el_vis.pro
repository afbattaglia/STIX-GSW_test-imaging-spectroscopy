;+
;
; NAME:
;   stx_el_spectra_2_el_vis.pro
;
; PURPOSE:
;     take in input the regularized visibility electron spectra and rearrange
;     the electron visibilities in standard array of visibility structures.
;
; CALLING SEQUENCE:
;     stx_el_spectra_2_el_vis, vis, countspectra, spectraptr, reg_el_vis, count_vis, reg_count_vis
;
; CALLED BY:
;
;   - stx_visibilities_inversion.pro
;
; CALLS TO:
;
;   - stx_vis_converter.pro
;
; INPUTS:
;
;   vis:             - Array of count visibility structures
;   countspectra:    - A 7-tag structure containing all the information about the count visibility spectra inverted by the software
;   spectraptr       - A pointer to 2-tag structures containing all information about the inverted regularized electron spectra
;
; OUTPUTS
;
;   reg_el_vis:  array of electron visibility structures (for the basic regularized solution)
;   stripvisptr: array of pointers to arrays of electron visibility structures (for each solution in the strip of confidence)
;   orig_ph_vis: array of the count visibility structures used as input for the inversion
;   reg_ph_vis:  array of regularized count visibility structures corresponding to the recovered electron visbilities
;
;+

pro stx_vis_el_spectra_2_el_vis, vis, countspectra, spectraptr, reg_el_vis, count_vis, reg_count_vis

DEFAULT, el_energy_max_factor,  2
DEFAULT, photon_bin_position,   0.5

  ;;;;;;;; Check for the count energy ranges of the visibility bag
  ct_edges_lower  = vis.energy_range[0]
  ct_edges_upper  = vis.energy_range[1]
  ct_edges_lower1 = ct_edges_lower[UNIQ(ct_edges_lower, SORT(ct_edges_lower))]
  ct_edges_upper1 = ct_edges_upper[UNIQ(ct_edges_upper, SORT(ct_edges_upper))]
  
  epsmin   = min(ct_edges_lower1)
  deltaeps = ct_edges_upper1 - ct_edges_lower1
  neps     = n_elements(ct_edges_lower1)
  bin_unif = deltaeps[where(ct_edges_lower1 le 16.)]
  epsini = ct_edges_lower1 + deltaeps*photon_bin_position

  ;;;;;;;; Count energy sampling
  eps = [ct_edges_lower1, ct_edges_upper1[-1]]

  ;;;;;;;; Electron energy sampling
  eemax = el_energy_max_factor*epsini[-1]
  
  e_stx_bin = [indgen(13)+4., indgen(3.)*2.+18., indgen(2.)*3.+25., indgen(3.)*4.+32., indgen(2.)*5.+45., indgen(2.)*7.+56., 70., 76., 84., 100., 120., 150.]
  delta_stx = [fltarr(13)+1., fltarr(3.)+2., fltarr(2.)+3., fltarr(3.)+4., fltarr(2.)+5., indgen(2.)+6., 7., 6., 8.,16., 20., 30.]
  
  ee_sampl = eps
  delta_e  = deltaeps
  index_min = min(where(e_stx_bin gt eps[-1]))
  index_max = min(where(e_stx_bin gt eemax))
 
  if (index_min eq -1) then begin
    ee_sampl = eps
    delta_e = deltaeps
  endif else begin
    if (index_max eq -1) then begin
      ee_sampl = [ee_sampl, e_stx_bin[index_min:-1]]
      delta_e = [delta_e, delta_stx[index_min:-1]]
    endif else begin
      ee_sampl = [ee_sampl, e_stx_bin[index_min:index_max]]
      delta_e = [delta_e, delta_stx[index_min:index_max]]
    endelse
  endelse
  
  ee  = ee_sampl + delta_e*photon_bin_position
  
  Nspectra = N_elements(countspectra.det)

  det_used = vis[0:n_elements(vis)/neps-1].isc
  
  viscont=1
  reg_el_vis = []
  econt = 0

  for econt =0,n_elements(ee)-1 do begin
    
  viscont=0
    
    for i=0,Nspectra-1 do begin
        outdata =*spectraptr[i]
        ind_e_re = where(outdata.stripRe[*,0] eq ee[econt],countRe)
        ind_e_im = where(outdata.stripIm[*,0] eq ee[econt],countIm)
        if (countRe eq 1 and countIm eq 1) then viscont = viscont+1
        if (viscont eq 1 and countRe eq 1 and countIm eq 1) then begin
          u           = countspectra.uv[0,i]
          v           = countspectra.uv[1,i]
          realvis     = outdata.stripRe[ind_e_re,1]
          imvis       = outdata.stripIm[ind_e_im,1]
          errRe       = 0.5*(max(outdata.stripRe[ind_e_re,2:11]) - min(outdata.stripRe[ind_e_re,2:11]))
          errIm       = 0.5*(max(outdata.stripIm[ind_e_im,2:11]) - min(outdata.stripIm[ind_e_im,2:11]))
          err         = 1./sqrt(2)*sqrt(errRe^2. + errIm^2.)
          det_used    = countspectra.det[i]
          lab_detused = countspectra.lab[i]
        endif
        if (viscont GT 1 and countRe eq 1 and countIm eq 1) then begin
          u           = [u,countspectra.uv[0,i]]
          v           = [v,countspectra.uv[1,i]]
          realvis     = [realvis, outdata.stripRe[ind_e_re,1]]
          imvis       = [imvis, outdata.stripIm[ind_e_im,1]]
          errRe       = 0.5*(max(outdata.stripRe[ind_e_re,2:11]) - min(outdata.stripRe[ind_e_re,2:11]))
          errIm       = 0.5*(max(outdata.stripIm[ind_e_im,2:11]) - min(outdata.stripIm[ind_e_im,2:11]))
          err         = [err, 1./sqrt(2)*sqrt(errRe^2. + errIm^2.)]
          det_used    = [det_used, countspectra.det[i]]
          lab_detused = [lab_detused, countspectra.lab[i]]
        endif
    end

    if (econt eq 0) then begin
      reg_el_vis = stx_vis_converter(ee_sampl[0],ee_sampl[1], u, v, realvis*(ee_sampl[1]-ee_sampl[0]), imvis*(ee_sampl[1]-ee_sampl[0]), err*(ee_sampl[1]-ee_sampl[0]), vis[0].time_range, vis[0].xyoffset, det_used, lab_detused)
    endif else begin
      visE = stx_vis_converter(ee_sampl[econt],ee_sampl[econt+1.], u, v,  realvis*(ee_sampl[econt+1]-ee_sampl[econt]),  imvis*(ee_sampl[econt+1]-ee_sampl[econt]),  err*(ee_sampl[econt+1]-ee_sampl[econt]), vis[0].time_range, vis[0].xyoffset, det_used, lab_detused)
      reg_el_vis = [reg_el_vis, visE]
    endelse
    
  endfor

  ;;;;;;;;;;; Repeat for count visibilities

  for ieps=0,neps-1 do begin
    viscont=0
    for i=0,Nspectra-1 do begin
        outdata =*spectraptr[i]
        ind_eps_re = where(outdata.residualsRe[*,0] eq eps[ieps]+deltaeps[ieps]/2.,countRe)
        ind_eps_im = where(outdata.residualsIm[*,0] eq eps[ieps]+deltaeps[ieps]/2.,countIm)
        if (countRe eq 1 and countIm eq 1) then viscont = viscont+1
        if (viscont eq 1 and countRe eq 1 and countIm eq 1) then begin
          u            = countspectra.uv[0,i]
          v            = countspectra.uv[1,i]
          countrealvis = outdata.residualsRe[ind_eps_re,1]
          countimvis   = outdata.residualsIm[ind_eps_im,1]
          err          = outdata.residualsRe[ind_eps_re,2]
          regcountRe   = outdata.residualsRe[ind_eps_re,3]
          regcountIm   = outdata.residualsIm[ind_eps_im,3]
          det_used     = countspectra.det[i]
          lab_detused  = countspectra.lab[i]
        endif
        if (viscont GT 1 and countRe eq 1 and countIm eq 1) then begin
          u            = [u,countspectra.uv[0,i]]
          v            = [v,countspectra.uv[1,i]]
          countrealvis = [countrealvis, outdata.residualsRe[ind_eps_re,1]]
          countimvis   = [countimvis, outdata.residualsIm[ind_eps_im,1]]
          err          = [err, outdata.residualsRe[ind_eps_re,2]]
          regcountRe   = [regcountRe, outdata.residualsRe[ind_eps_re,3]]
          regcountIm   = [regcountIm, outdata.residualsIm[ind_eps_im,3]]
          det_used     = [det_used, countspectra.det[i]]
          lab_detused  = [lab_detused, countspectra.lab[i]]
        endif
    end
    if (ieps eq 0) then begin
      count_vis     = stx_vis_converter(eps[0], eps[1], u, v, countrealvis*(eps[1]-eps[0]), countimvis*(eps[1]-eps[0]), err*(eps[1]-eps[0]), vis[0].time_range, vis[0].xyoffset, det_used, lab_detused)
      reg_count_vis = stx_vis_converter(eps[0], eps[1], u, v, regcountRe*(eps[1]-eps[0]), regcountIm*(eps[1]-eps[0]), err*(eps[1]-eps[0]), vis[0].time_range, vis[0].xyoffset, det_used, lab_detused)
    endif else begin
      viscount      = stx_vis_converter(eps[ieps], eps[ieps+1], u, v, countrealvis*(eps[ieps+1]-eps[ieps]), countimvis*(eps[ieps+1]-eps[ieps]), err*(eps[ieps+1]-eps[ieps]), vis[0].time_range, vis[0].xyoffset, det_used, lab_detused)
      reg_viscount  = stx_vis_converter(eps[ieps], eps[ieps+1], u, v, regcountRe*(eps[ieps+1]-eps[ieps]), regcountIm*(eps[ieps+1]-eps[ieps]), err*(eps[ieps+1]-eps[ieps]), vis[0].time_range, vis[0].xyoffset, det_used, lab_detused)
      count_vis     = [count_vis,viscount]
      reg_count_vis = [reg_count_vis,reg_viscount]
    endelse
  endfor

end