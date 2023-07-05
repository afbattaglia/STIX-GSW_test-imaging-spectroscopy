;+
;
; NAME:
;   stx_vis_reg_ge_rescaling.pro
;
; PURPOSE:
;     for the input count visibility spectrum:
;     1) check the sampling, check the high and low energy values
;     2) rescale the input count visibility spectrum by using a single PL, a broken PL or a sum of PLs.
;     3) reject samples at low energies if their error bars are too big otherwise the regularization parameter
;        by using the discrepancy principle cannot be found and the spectrum cannot be inverted.
;
; CALLING SEQUENCE:
;     stx_vis_reg_ge_rescaling, count_vis_spectra.energyl, count_vis_spectra.energyh, vis_array, $
;                               count_vis_spectra.error[*,l], drmini, el_energy_max_factor, $
;                               nph, eps, gs, err, wei, fitph, fitel, drm
;
; CALLED BY:
;   - stx_vis_regularized_inversion.pro
;
; CALLS TO:
;   - stx_vis_reg_ge_cross3bn.pro
;
; INPUTS:
;   count_vis_spectra.energyl:    array of the lower energies (keV) of count bins
;   count_vis_spectra.energyh:    array of the higher energies (keV) of count bins
;   vis_array:                    for a fixed (u,v) pair, array of the visibility values (real or imaginary) as a
;                                 function of count energy
;   count_vis_spectra.error[*,l]: corresponding error bars.
;   drmini:                       the drm computed by stx_build_drm.pro
;   el_energy_max_factor:         the maximum electron energy (in units of the maximum count energy provided)
;                                 in the returned electron energy spectrum array (i.e., Emax=el_energy_max_factor*eps_max).
;                                  The default is 2.
;  dist_solo_sun:                 SolO - Sun distance
;
;
; OUTPUTS
;   nph:         number of count energy points in the spectrum to be inverted
;   eps:         count energies points (keV)
;   gs:          visibility spectrum (real or imaginary part)  (cnt/cm^2/s/keV)
;   err:         errors in visibility values                   (cnt/cm^2/s/keV)
;   wei:         weight vector
;   fitph:       fit of the count spectrum to rescale it
;   fitel:       rescaling of the electron spectrum according to the rescling performed in the count space
;   drm:         sub-matrix of the initial drmini selected according to the selected count energy points.
;
;+


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fit_function_PL1, X, Apar, F, pder

  F = Apar[0] * x^(-Apar[1])

  IF N_PARAMS() GE 4 THEN $
    pder = [[x^(-Apar[1])], [-Apar[0] * x^(-Apar[1]) * alog(X)] ]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fit_function_PL2, X, Apar, F, pder

  F = -Apar[0] * x^(-Apar[1])

  IF N_PARAMS() GE 4 THEN $
    pder = [[-x^(-Apar[1])], [Apar[0] * x^(-Apar[1]) * alog(X)] ]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fit_function_SPL, X, Apar3, F, pder

  F = Apar3[0] * x^(-Apar3[1])-Apar3[2] * x^(-Apar3[3])

  IF N_PARAMS() GE 4 THEN $
    pder = [[x^(-Apar3[1])], [-Apar3[0] * x^(-Apar3[1]) * alog(X)],[-x^(-Apar3[3])], [Apar3[2] * x^(-Apar3[3]) * alog(X)] ]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO fit_function_BPL, X, AparBPL, F, pder

  Eth=x[0]+(x[N_ELEMENTS(x)-1]-x[0])/5.
  F=dblarr(N_ELEMENTS(x))

  for i=0,N_ELEMENTS(x)-1 do begin
    if(x(i) lE Eth) then F(i) = AparBPL[0] * x(i)^(-AparBPL[1]) else F(i) = AparBPL[0] * Eth^(AparBPL[2]-AparBPL[1])*x(i)^(-AparBPL[2])
  end

  IF N_PARAMS() GE 4 THEN begin

    pder=dblarr(N_ELEMENTS(x),N_ELEMENTS(AparBPL))
    for i=0,N_ELEMENTS(x)-1 do begin
      if(x[i] lE Eth ) then begin
        pder[i,0] = x[i]^(-AparBPL[1])
        pder[i,1] = -AparBPL[0] * x[i]^(-AparBPL[1]) * alog(X[i])
        pder[i,2] = 0.
      endif else begin
        pder[i,0] = Eth^(AparBPL[2]-AparBPL[1])*x[i]^(-AparBPL[2])
        pder[i,1] = -AparBPL[0] * Eth^(AparBPL[2]-AparBPL[1])*x[i]^(-AparBPL[2]) * alog(Eth)
        pder[i,2] = AparBPL[0] * Eth^(AparBPL[2]-AparBPL[1])*x[i]^(-AparBPL[2])*(alog(Eth)-alog(X[i]))
      endelse
    end

  end

END


function rescaling_lookup_table,gamma

  ptab=fltarr(17)
  ptab[0]=0
  ptab[1]=0
  ptab[2]=0.1
  ptab[3]=0.2
  ptab[4]=0.35
  ptab[5]=0.5
  ptab[6]=0.65
  ptab[7]=0.75
  ptab[8]=0.85
  ptab[9]=0.95
  ptab[10]=1.1
  ptab[11]=1.2
  ptab[12]=1.4
  ptab[13]=1.55
  ptab[14]=1.6
  ptab[15]=1.7
  ptab[16]=1.75

  gamma1=fltarr(17)
  gamma1[0]=2.20
  gamma1[1]=2.38
  gamma1[2]=2.58
  gamma1[3]=2.80
  gamma1[4]=3.03
  gamma1[5]=3.26
  gamma1[6]=3.53
  gamma1[7]=3.79
  gamma1[8]=4.02
  gamma1[9]=4.29
  gamma1[10]=4.54
  gamma1[11]=4.79
  gamma1[12]=5.06
  gamma1[13]=5.32
  gamma1[14]=5.55
  gamma1[15]=5.82
  gamma1[16]=6.09

  gammainterp = FINDGEN(1000)/250.+2.20
  pinterp=spline(gamma1,ptab,gammainterp)

  distmin=abs(gamma-gammainterp(0))
  for i=1,999 do begin
    dist=abs(gamma-gammainterp(i))
    if (dist le distmin) then distmin=dist else break
  end

  if(pinterp[i-1] lt 0) then pinterp[i-1]=0

  return, pinterp[i-1]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



pro stx_vis_reg_ge_rescaling,photon_energy_low_array, photon_energy_high_array, vis_array, $
                             vis_uncertainty_array, drmini, el_energy_max_factor, dist_solo_sun, $
                             nph, eps, gs, err, wei, fitph, fitel, drm

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  N = N_elements(photon_energy_low_array)
  epsini = dblarr(N)

  ;;;;;;;; Energy bin size 
  photon_bin_position = 0.5
  bin = photon_energy_high_array - photon_energy_low_array
  bin_unif = bin[where(photon_energy_low_array le 16.)]
  
  photon_energy_edges = [photon_energy_low_array, photon_energy_high_array[-1]]
  
  epsini = photon_energy_low_array + bin*photon_bin_position

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; Input visibilities are in units of (cm^-2 s^-1). Dividing by the
  ;;; binning makes them in the correct units (cm^-2 s^-1 keV^-1) for inversion.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  vis_array             = vis_array/bin
  vis_uncertainty_array = vis_uncertainty_array/bin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Check for high energy values: cancel values for which the uncertainty = 0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Nend = max(where(vis_uncertainty_array ne 0))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Check for low energy values: eps(0) must be around 9
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minphotonenergy = 9
  Nstart = min(where(epsini ge minphotonenergy))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;; Rescaling
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ;;;;;;;;;;; Fit the visibility spectrum with a single power law

  if (vis_array[0] GT vis_array[Nend]) then begin
    AparPL = [vis_array[0]*epsini[0]^4,4]
    zPL    = curvefit(epsini,vis_array,vis_uncertainty_array,AparPL,SIGMA, CHISQ=chiquadroPL, FUNCTION_NAME='fit_function_PL1', ITMAX=1000, status=status)
    ;print, 'Chiquadro=', chiquadroPL , '   A=', AparPL
    flagel=1
  endif else begin
    AparPL = [-vis_array[0]*epsini[0]^4,4]
    zPL    = curvefit(epsini,vis_array, vis_uncertainty_array, AparPL, SIGMA, CHISQ=chiquadroPL, FUNCTION_NAME='fit_function_PL2', ITMAX=1000, status=status)
    ;print,'Chiquadro=', chiquadroPL , '   A=', AparPL
    flagel=-1
  endelse

  ;;;;;;;;;;;;;;;;;;;;;;;; Fit with sum of Power Laws
  gamma1=4.
  gamma2=2.
  AparSPL = [vis_array[0]*epsini[0]^gamma1,gamma1,-vis_array[N-1]*epsini[N-1]^gamma2,gamma2]
  zSPL    = curvefit(epsini,vis_array,vis_uncertainty_array,AparSPL,SIGMA, CHISQ=chiquadroSPL, FUNCTION_NAME='fit_function_SPL', ITMAX=1000, status=status)
  ;print,'Chiquadro=', chiquadroSPL , '   A=', AparSPL
  pSPL1   = rescaling_lookup_table(AparSPL[1])
  pSPL2   = rescaling_lookup_table(AparSPL[3])

  ;;;;;;;;;;;;;;;;;;;;;;; Fit with a Broken Power Law
  Eth = epsini[0]+(epsini[Nend]-epsini[0])/5.
  ith = where( abs(epsini-Eth) eq min(abs(epsini-Eth)) )
  if (vis_array[0] GT vis_array[Nend]) then begin
    AparBPL = [vis_array[0]*epsini[0]^4,4,5]
    zBPL    = curvefit(epsini,vis_array,vis_uncertainty_array,AparBPL,SIGMA, CHISQ=chiquadroBPL, FUNCTION_NAME='fit_function_BPL', ITMAX=1000,  status=status)
    ;print, 'Chiquadro BPL+   =', chiquadroBPL , '   A=', AparBPL
    pBPL1   = rescaling_lookup_table(AparBPL[1])
    pBPL2   = rescaling_lookup_table(AparBPL[2])
  endif else begin
    AparBPL = [-vis_array[0]*epsini[0]^4,4,5]
    zBPL    = curvefit(epsini,vis_array,vis_uncertainty_array,AparBPL,SIGMA, CHISQ=chiquadroBPL, FUNCTION_NAME='fit_function_BPL', ITMAX=1000, status=status)
    ;print, 'Chiquadro BPL-   =', chiquadroBPL , '   A=', AparBPL
    pBPL1   = rescaling_lookup_table(AparBPL[1])
    pBPL2   = rescaling_lookup_table(AparBPL[2])
  endelse
  pBPL=rescaling_lookup_table((AparBPL[1]+AparBPL[2])/2.)
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;; Look for the minimum Chisquare
  chiquadro = [chiquadroPL,chiquadroSPL,chiquadroBPL]
  whichfit  = where(chiquadro EQ min(chiquadro))
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  visfit    = dblarr(N)
  
  print, ' '
  CASE whichfit OF
     0: print, 'Fit with Single PL'
     1: print, 'Fit with sum of PL'
     2: print, 'Fit with Broken PL'
  ENDCASE
  print, ' '

  CASE whichfit OF
    0: if (flagel eq 1) then begin
      p = rescaling_lookup_table(AparPL[1])
      visfit=AparPL[0]*epsini^(-p)
    endif else begin
      p = rescaling_lookup_table(AparPL[1])
      visfit=-AparPL[0]*epsini^(-p)
    endelse
    1: visfit = AparSPL[0] * epsini^(-pSPL1)-AparSPL[2] * epsini^(-pSPL2)
    2: visfit = AparBPL[0] * epsini^(-pBPL)
  ENDCASE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Visibilities at low energy have big error bars. This can be a problem in order to
  ; find the regularization parameter. In fact the reg parameter in a first step is
  ; found using the discrepancy principle, but if the error bars are too big
  ; it may happen that for each lambda ||Af-g||<||err|| ==> optimal lambda does
  ; not exist. In the following, we fix a very high value for lambda (1.e15) and
  ; we look for the minimum energy verifing the condition ||Af-g||>||err|| ==>
  ; i.e. it is possible to find lambda.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ;;;;; select in the drm matrix the submatrix corresponding to energy between eps(Nstart)-eps(Nend)

  drm0  = drmini[Nstart:Nend,Nstart:Nend]

  mu    = 1.e15
  eemax = el_energy_max_factor*epsini[Nend]

  Nstart0 = Nstart

  ;;;;;;

  nph = Nend-Nstart+1

  arg0 = dblarr(nph)

  eps0      = epsini[Nstart:Nend]
  vis0      = vis_array[Nstart:Nend]
  deltavis0 = vis_uncertainty_array[Nstart:Nend]
  visfit0   = visfit[Nstart:Nend]

  ;;;;;;;;;;rescaling

  vis0      = vis0/visfit0
  deltavis0 = deltavis0/visfit0

  ;;;;;;;;;;;;;;;;;;;;
  wei0 = fltarr(nph)+bin

  ;;;;;;;;;;;;;; energy sampling 
  e_stx_bin = [indgen(13)+4., indgen(3.)*2.+18., indgen(2.)*3.+25., indgen(3.)*4.+32., indgen(2.)*5.+45., indgen(2.)*7.+56., 70., 76., 84., 100., 120., 150.]
  delta_stx = [fltarr(13)+1., fltarr(3.)+2., fltarr(2.)+3., fltarr(3.)+4., fltarr(2.)+5., indgen(2.)+6., 7., 6., 8.,16., 20., 30.]

  ee_sampl = photon_energy_edges
  delta_e  = bin
  index_min = min(where(e_stx_bin gt photon_energy_edges[-1]))
  index_max = min(where(e_stx_bin gt eemax))
  if (index_min eq -1) then begin
    ee_sampl = photon_energy_edges
    delta_e = bin
  endif else begin
    if (index_max eq -1) then begin
      ee_sampl = [ee_sampl, e_stx_bin[index_min:-1]]
      delta_e = [delta_e, delta_stx[index_min:-1]]
    endif else begin
      ee_sampl = [ee_sampl, e_stx_bin[index_min:index_max-1]]
      delta_e = [delta_e, delta_stx[index_min:index_max]]      
    endelse
  endelse

  ee  = ee_sampl + delta_e*photon_bin_position
  nee = n_elements(ee)
  visfitel = dblarr(nee)

  CASE whichfit OF
    0: if (flagel eq 1) then visfitel = AparPL[0]*ee^(-p) else visfitel = -AparPL[0]*ee^(-p)
    1: visfitel = AparSPL[0] * ee^(-pSPL1)-AparSPL[2] * ee^(-pSPL2)
    2: visfitel = AparBPL[0] * ee^(-pBPL)
  ENDCASE

  eebrem = 0
  Z = 1.2
  stx_vis_reg_ge_cross3bn,eebrem,nph,nee,eps0,ee,drm0,visfit0,visfitel,Z,dist_solo_sun, w,u,sf

  ;;;;;;;;; Rescaling of the singular vectors and
  ;;;;;;;;; values according to the weights
  for k=0,nph-1 do begin
    w[k]   = sqrt(wei0[k])*w[k]
    u[*,k] = u[*,k]/sqrt(wei0[k])
  end
  ;;;;;;;;;;;;;
  for k=0,nph-1 do begin
    scal    = abs((wei0*vis0)##u[k,*])
    arg0[k] = (mu*scal/(w[k]*w[k]+mu))*(mu*scal/(w[k]*w[k]+mu))
  end

  rerr = total(deltavis0*deltavis0*wei0)

  while (rerr gt total(arg0)) do begin

    Nstart+=1
    nph    = Nend-Nstart+1

    drm1 = drmini[Nstart:Nend,Nstart:Nend]

    if( nph le 3) then break

    arg0 = dblarr(nph)

    eps1      = epsini[Nstart:Nend]
    vis1      = vis_array[Nstart:Nend]
    deltavis1 = vis_uncertainty_array[Nstart:Nend]
    wei1      = wei0[Nstart-Nstart0:Nend-Nstart0]
    visfit1   = visfit[Nstart:Nend]

    ;;;;;;;;;;rescaling

    vis1      = vis1/visfit1
    deltavis1 = deltavis1/visfit1

    ;;;;;;;;;;;;;;;;;;;;

    rerr = total(deltavis1*deltavis1*wei1)

    ee = ee[1:-1]
    nee = n_elements(ee)
    visfitel = dblarr(nee)


    CASE whichfit OF
      0: if (flagel eq 1) then visfitel1 = AparPL[0]*ee^(-p) else visfitel1 = -AparPL[0]*ee^(-p)
      1: visfitel1 = AparSPL[0] * ee^(-pSPL1)-AparSPL[2] * ee^(-pSPL2)
      2: visfitel1 = AparBPL[0] * ee^(-pBPL)
    ENDCASE

    eebrem=0
    stx_vis_reg_ge_cross3bn,eebrem,nph,nee,eps1,ee,drm1,visfit1,visfitel1,Z,dist_solo_sun, w,u,sf

    ;;;;;;;;; Rescaling of the singular vectors and
    ;;;;;;;;; values according to the weights
    for k=0,nph-1 do begin
      w[k]   = sqrt(wei1[k])*w[k]
      u[*,k] = u[*,k]/sqrt(wei1[k])
    end
    ;;;;;;;;;;;;;
    for k=0,nph-1 do begin
      scal    = abs((wei1*vis1)##u[k,*])
      arg0[k] = (mu*scal/(w[k]*w[k]+mu))*(mu*scal/(w[k]*w[k]+mu))
    end

  endwhile

  if( nph le 3) then  return

  ;;;;;;;;;;;;;;;;;;;;;
  nph = Nend-Nstart+1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;   Building the weights vector according to the sampling
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wei = fltarr(nph)+bin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;   Saving values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  eps   = fltarr(nph)
  gs    = fltarr(nph)
  err   = fltarr(nph)
  fitph = fltarr(nph)

  for i=Nstart,Nend do begin
    eps[i-Nstart] = epsini[i]
    gs[i-Nstart]  = vis_array[i]
    err[i-Nstart] = vis_uncertainty_array[i]
    fitph[i-Nstart] = visfit[i]
  end
  drm = drmini[Nstart:Nend,Nstart:Nend]

  eemax = el_energy_max_factor*eps[nph-1]

  ee_sampl = photon_energy_edges
  delta_e  = bin
  index_min = min(where(e_stx_bin gt photon_energy_edges[-1]))
  index_max = min(where(e_stx_bin gt eemax))
  if (index_min eq -1) then begin
    ee_sampl = photon_energy_edges
    delta_e = bin
  endif else begin
    if (index_max eq -1) then begin
      ee_sampl = [ee_sampl, e_stx_bin[index_min:-1]]
      delta_e = [delta_e, delta_stx[index_min:-1]]
    endif else begin
      ee_sampl = [ee_sampl, e_stx_bin[index_min:index_max-1]]
      delta_e = [delta_e, delta_stx[index_min:index_max]]      
    endelse
  endelse
   
  ee  = ee_sampl + delta_e*photon_bin_position

  nee = n_elements(ee)
  fitel = dblarr(nee)

  CASE whichfit OF
    0: if (flagel eq 1) then fitel = AparPL[0]*ee^(-p) else fitel =-AparPL[0]*ee^(-p)
    1: fitel = AparSPL[0] * ee^(-pSPL1)-AparSPL[2] * ee^(-pSPL2)
    2: fitel = AparBPL[0] * ee^(-pBPL)
  ENDCASE

end

