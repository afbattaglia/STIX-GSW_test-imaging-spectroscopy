;+
;
; NAME:
;   stx_build_count_vis_spectra.pro
;
; PURPOSE:
;     takes in input a standard array of visibility structures and
;     1) for each u,v pair build a count visibility spectrum as a function of the count energy.
;        A spectrum is accepetd only if it is made of 5 points at least.
;     2) save the spectra in the structure "count_vis_spectra"
;
; CALLING SEQUENCE:
;     stx_build_count_vis_spectra, detOK, vis, count_vis_spectra
;
; CALLED BY:
;
;   - stx_visibilities_inversion.pro
;
; INPUTS:
;    vis  : Array of count visibility structures.
;    detOK: Detectors involved in the inversion
;
; OUTPUTS:
;   count_vis_spectra: a 9-tag structure containing
;                            - energyl       dbalrr(neps)            Lower energy of each count bin                  (keV)
;                            - energyh       dbalrr(neps)            Upper energy of each count bin                  (keV)
;                            - RealPart      dblarr(neps,Nspectra)   Visibility spectra (Real Part) vs count energy  (cnt/cm^2/s/keV)
;                            - ImagPart      dblarr(neps,Nspectra)   Visibility spectra (Imag Part) vs count energy  (cnt/cm^2/s/keV)
;                            - error         dblarr(neps,sNspectra)  Errors in visibility values (= for both parts)  (cnt/cm^2/s/keV)
;                            - det           intarr(Nspectra)        Detector associated to each spectrum
;                            - uv            dblarr(2,Nspectra)      u,v pair associated to each spectrum
;
;+


pro stx_build_count_vis_spectra, detOK, vis, count_vis_spectra

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  threshold = 5     ;;;;; minimum number of samples to build a spectrum for the inversion
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;; Check for the count energy ranges of the visibility bag 
  ct_edges_lower  = vis.energy_range[0]
  ct_edges_upper  = vis.energy_range[1]
  ct_edges_lower1 = ct_edges_lower[UNIQ(ct_edges_lower, SORT(ct_edges_lower))]
  ct_edges_upper1 = ct_edges_upper[UNIQ(ct_edges_upper, SORT(ct_edges_upper))]
  deltaeps = ct_edges_upper1 - ct_edges_lower1

  epsmin   = min(vis.energy_range[0])
  deltaeps = ct_edges_upper1 - ct_edges_lower1
  neps     = n_elements(ct_edges_lower1)

  ;;;;;;;; Count energy sampling
  eps = [ct_edges_lower1, ct_edges_upper1[-1]]

  ;;;;;;;;;;;; The number of spectra to be inverted is = to the number of subcollimators used
  Nspectra = total(detOK,/INTEGER)
  det_used = vis[0:n_elements(vis)/neps-1].isc
  label_used =  vis[0:n_elements(vis)/neps-1].label
  
  ;;;;;;;;;;; Construct count_vis_spectra structure
  count_vis_spectra = { energyl:  eps[0:neps-1],$              ; Lower energy of each count bin (keV)
                        energyh:  eps[0:neps-1]+deltaeps,$     ; Upper energy of each count bin (keV)
                        RealPart: dblarr(neps,Nspectra),$      ; Visibility spectra (Real Part) vs count energy (cnt/cm^2/s/keV)
                        ImagPart: dblarr(neps,Nspectra),$      ; Visibility spectra (Imag Part) vs count energy (cnt/cm^2/s/keV)
                        error:    dblarr(neps,Nspectra),$      ; Errors in visibility values (= for both parts)  (cnt/cm^2/s/keV)
                        det:      det_used, $                  ; Detector associated to each spectrum
                        lab:      label_used, $                ; Label of detector associated to each spectrum
                        uv:       transpose([[vis[0:n_elements(vis)/neps-1].u], [vis[0:n_elements(vis)/neps-1].v]]) } ; u,v pair associated to each spectrum
                        
  if (neps lt threshold) then begin
    print, ' '
    message, "Please use at least 5 energy bins."
    print, ' '
  endif else begin
    for i = 0, neps -1 do begin
      index = where(vis.energy_range[0] eq eps[i], loc)
      vis_j = vis[index]
      count_vis_spectra.RealPart[i,*] = float(vis_j.obsvis)
      count_vis_spectra.ImagPart[i,*] = imaginary(vis_j.obsvis)
      count_vis_spectra.error[i,*]    = vis_j.sigamp
    endfor
  endelse

end
