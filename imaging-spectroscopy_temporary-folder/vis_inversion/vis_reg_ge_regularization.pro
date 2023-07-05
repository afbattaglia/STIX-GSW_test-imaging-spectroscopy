;+
;
; NAME:
;   vis_reg_ge_regularization.pro
;
; PURPOSE:
;       computes the regularized solution. It is divided in two parts:
;           the first one calculates the regularization parameter while
;           the second one the solution.
;
; CALLING SEQUENCE:
;
; vis_reg_ge_regularization,nph,nee,wei,eps,gs,err,u,w,sf,fitph,fitel,$
;    ee,opt,gopt,regsol,residual_array, cumulativeresidual_array
;
; CALLED BY:
;
;   - stx_vis_regularized_inversion.pro
;
; CALLS TO:
;
;   none
;
; INPUTS:
;
;   nph:    number of count energy bins used
;   nee:    number of electron energy bins used
;   wei:    weight vector
;   eps:    count energy vector
;   gs:     count visibility spectrum (cnts/cm^2/s/kev at detector)
;   err:    uncertainty in visibility value
;   u:      singular vectors of the svd of the bremsstrahlung integral operator
;   w:      singular values of the svd of the bremsstrahlung integral operator
;   sf:     singular functions of the svd of the bremsstrahlung integral operator
; fitph:    rescaling of the count visibility spectrum
; fitel:    rescaling of the electron visibility spectrum
;   ee:     electron energy vector
;
;
; OUTPUTS:
;
; opt:                      value of the optimal regularization parameter
; gopt:                     count visibility spectrum corresponding to the recovered electron visibility array
; regsol:                   the recovered electron visibility array
; residual_array:           residual count visibility spectrum (actual - recovered)
; cumulativeresidual_array: cumulative residual, defined as : c_j = (1/j) sum_[i=1]^j res_i
;
; RESTRICTIONS:
;  - CHOICE OF THE REGULARIZATION PARAMETER:
;     The code uses the cumulative residuals to fix the optimal regularization parameter.
;     The cumulative residuals can be used as a measure of residual 'clustering'; for an
;     acceptable fit to the data, abs(C_j) should be less than n/sqrt(j), where n is the
;     number of allowed standard deviations (n=1 is used in the code to drive the solution
;     to an 'acceptable' solution).
;
;-




pro vis_reg_ge_regularization,nph,nee,wei,eps,gs,err,u,w,sf,fitph,fitel,$
    ee,opt,gopt,regsol, residual_array, cumulativeresidual_array

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;; Regularization parameter choice: the program starts with the
;;;;;;;;;;;;;;;;;; discrepancy principle and in sequence with the 3-sigma criterion.
;;;;;;;;;;;;;;;;;; If residuals corresponding to the lambda chosen with the discrepancy
;;;;;;;;;;;;;;;;;; principle are not inside 3-sigma the parameter is further decreased.


;;;;;;;;;;;;;;;;;; DISCREPANCY PRINCIPLE


Nmu=1000.

arg=dblarr(nph,Nmu)
mu=dblarr(Nmu)        ;;;;;;; regularization parameter values vector
discr=dblarr(Nmu)     ;;;;;;; discrepancy values vector

mu_min=1.e-4    ;Initial minimum value for the regularization parameter

discrep=0
iterazione=0
left=0
right=0

rerr=total(err*err*wei)

while discrep EQ 0 do begin
    mu_max=1.e2*mu_min

    step=(mu_max-mu_min)/(Nmu-1)

    for i=0,Nmu-1 do mu[i]=mu_min+i*step

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;; for all cross sections use WEI in the discrepancy
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    for k=0,nph-1 do begin
       scal=abs((wei*gs)##u[k,*])
       for i=0,Nmu-1 do arg[k,i]=(mu[i]*scal/(w[k]*w[k]+mu[i]))*(mu[i]*scal/(w[k]*w[k]+mu[i]))
    end

    discr=total(arg,1)-rerr


    opt=mu[where(abs(discr) eq min(abs(discr)))]
    ;opt = opt[0]
    lambda=opt[0]

    if (lambda eq mu[0]) then begin
        mu_min=1.e-2*mu_min
        left=left+1
    endif
    if (lambda eq mu[Nmu-1]) then begin
        mu_min=1.e2*mu_min
        right=right+1
    endif
    if (lambda GT mu[0] and lambda LT mu[Nmu-1]) then discrep=1

    if (left ne 0 and right ne 0) then discrep=1

endwhile

;;;;;;;;;;;;;;;;;;;;; END OF THE DISCREPANCY PRINCIPLE'

;;;;;;;;;;;;;;;;;;;;; 3-SIGMA CRITERION

optmin=1.e-4*lambda
iterations=100.
deltaopt=(optmin/lambda)^(1/iterations)

arg1=dblarr(nph,max(nee))
regsol=dblarr(max(nee))
argopt=dblarr(nph,nph)
gopt=dblarr(nph)

n3sigma=floor(nph*0.01)

for l=0,iterations-1 do begin

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;; for all cross sections use WEI in the regularized solution formula
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    for k=0,nph-1 do begin
       scal=(wei*gs)##u[k,*]
       for i=0,nph-1 do argopt[k,i]=scal*w[k]*w[k]*u[k,i]/(w[k]*w[k]+opt)
    end
    gopt=total(argopt,1)

    ;;;;;;;;;;rescaling
    gs=gs*fitph
    err=err*fitph
    gopt=gopt*fitph
    ;;;;;;;;;;;;;;;;;;;;

    residual = gs - gopt

    ressum=fltarr(nph)
    for i=0,nph-1 do begin
       for j=0,i do ressum[i]=ressum[i]+(residual[j]/err[j])/float(i+1)
    endfor
    ressum3sig=fltarr(nph-1)
    for i=0,nph-2 do ressum3sig[i]=ressum(i+1)

    expected=1./sqrt(findgen(nph)+1)
    expected3sig=fltarr(nph-1)
    for i=0,nph-2 do expected3sig[i]=expected[i+1]

    ;;;;;;;;;;rescaling
    gs=gs/fitph
    err=err/fitph
    ;;;;;;;;;;;;;;;;;;;;

    count=0
    for i=0,nph-2 do begin
       if (expected3sig[i]-ressum3sig[i])<0 or (expected3sig[i]+ressum3sig[i])<0 then count=count+1
    end

    if (count gt n3sigma) then opt=opt*deltaopt else break

    if (min(opt) le optmin) then opt=lambda ;;;;; if after 100 iterations the 3-sigma criterion is not
                                       ;;;;; satisfied, chose the value obtained with the discrepancy
                                       ;total(abs(opt))

end ;;; End of loop on the iterations

;;;;;;;;;;;;;;;;;;;;; END OF THE 3-SIGMA CRITERION
;print, 'discr   =  ', lambda
;print, '3 sigma =', opt
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;; Compute the regularised solution and the residuals
;;;;;;;;;;;;;; corresponding to the selected optimal value of lambda.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for k=0,nph-1 do begin
    scal=(wei*gs)##u[k,*]
    for j=0,max(nee)-1 do arg1[k,j]=scal*w[k]*sf[k,j]/(w[k]*w[k]+opt)
end

regsol=total(arg1,1)
regsol=fitel*regsol ;;;; rescaling

;;;;;;;;;;;;;;;;;;;;; Residuals

gs=gs*fitph
err=err*fitph

normres=residual/err
residual_array=normres
cumulativeresidual_array=ressum

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end
