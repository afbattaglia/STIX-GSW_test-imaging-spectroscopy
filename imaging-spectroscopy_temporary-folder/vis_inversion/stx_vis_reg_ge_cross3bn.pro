;+
;
; NAME:
;   stx_vis_reg_ge_cross3bn.pro
;
; PURPOSE:
;       computes the Singular Value Decomposition of the integral operator associated
;       to Eq.3 in Piana et al. The Astrophysical Journal, 665:846-855, 2007 August 10.
;       This operator is derived from the fully relativistic, solid-angle-averaged,
;       cross section (Koch & Motz 1959, Rev. Mod. Phys.,31, 920-955) including the
;       Elwert collisional correction factor.
;       It optionally includes the term due to electron-electron bremsstrahlung
;
;
; CALLING SEQUENCE:
;
; stx_vis_reg_ge_cross3bn,eebrem,nph,nee,eps,ee,drm,visfit,visfitel,Z,w,u,sf
;
; CALLED BY:
;
;   - stx_vis_regularized_inversion.pro
;   - stx_vis_reg_ge_rescaling.pro
;
; CALLS TO:
;
;   - vis_reg_ge_bremss_cross.pro
;
; INPUTS:
;  Eebrem:   a number (0 or 1) corresponding to:
;                 0 = only cross3bn
;                 1 = cross3bn plus the term due to electron-electron bremsstrahlung
;  nph:      number of count energy bins used
;  nee:      number of electron energy bins used
;  eps:      count energy vector
;  ee:       electron energy vector
;  drm:      detector response matrix
;  visfit:   fit of the count visibility spectrum to rescale it
;  visfitel: rescaling of the electron visibility spectrum
;  Z:        value of the root-mean-square atomic number of the target. Default = 1.2
;
;
; OUTPUTS:
;
;    W: singular values of the SVD of the Bremsstrahlung integral operator
;
;    U: singular vectors of the SVD of the Bremsstrahlung integral operator
;
;   SF: singular functions of the SVD of the Bremsstrahlung integral operator
;
;
; RESTRICTIONS:
;
;      We caution the user that addition of the electron-electron bremsstrahlung term
;      increases the computational time considerably.  Unless photon energies significantly
;      above 100 keV are used, the results for 'cross3bn' and 'cross3bnee' are very similar.
;
;-


pro stx_vis_reg_ge_cross3bn,eebrem,nph,nee,eps,ee,drm,visfit,visfitel,Z,dist_solo_sun, w,u,sf

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;; Electron energies sampled starting from the same initial photon energy
  ;;;;;;;;;; and with the same bin


  eebin = fltarr(n_elements(ee))
  eebin[0] = ee[1]-ee[0]
  for i = 0, n_elements(ee)-2 do eebin[i+1] = ee[i+1]-ee[i]
  epsbin = fltarr(n_elements(eps))
  epsbin[0] = eps[1]-eps[0]
  for i = 0, n_elements(eps)-2 do epsbin[i+1] = eps[i+1]-eps[i]
  
  q=dblarr(nee,nph)

  for j=0,nee-1 do begin
    E=ee[j]
    for i=0,nph-1 do begin
      eph=eps[i]
      if (eebrem eq 0) then Qnew=vis_reg_ge_bremss_cross(E,eph,Z,/Noelec)
      if (eebrem eq 1) then Qnew=vis_reg_ge_bremss_cross(E,eph,Z)
      q[j,i]=Qnew
    end
  end


  counts=eps
  Ker=dblarr(max(nee),nph)

  R = dist_solo_sun  ; distance in cm
  R2pi4=4.*!PI*R^2

  for iq=0,nph-1 do begin
    for j=0,max(nee)-1 do begin
      Ker(j,iq)=0.
      for i=0,nph-1 do begin
        if( (eps[i] ge counts[iq]) and (eps[i] le ee[j])) then Ker[j,iq]= Ker[j,iq]+drm[iq,i]*q[j,i]*epsbin[i]
      end
      Ker[j,iq]=1.d50*eebin[j]*Ker[j,iq]/R2pi4
    end
  end

  ;;;;;;;;;; rescaling


  for iq=0,nph-1 do begin
    for j=0,max(nee)-1 do begin
      Ker[j,iq]=Ker[j,iq]*visfitel[j]/visfit[iq]
    endfor
  endfor


  ;;;;;;;;;;;;;;;;;;;;

  svdc,Ker,w,u,v,/double

  sf=v    ;;;; singular functions

end
