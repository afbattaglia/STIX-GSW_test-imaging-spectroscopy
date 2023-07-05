;+
;
; NAME:
;   vis_reg_ge_confidence.pro
;
; PURPOSE:
;      determines a 'set' of solutions in order to estimate the uncertainty in the recovered
;      solution.  Each realization is produced by randomly perturbing the input data
;      using the noise levels specified therein.
;
; CALLING SEQUENCE:
;
; vis_reg_ge_confidence,nph,nee,ee,regsol,$
;    err,gs,eps,wei,confidencestrip,u,w,opt,sf,drm,fitph,fitel,Z,datastruct
;
; CALLED BY:
;
;   - vis_regularized_inversion.pro
;
; CALLS TO:
;
;   none
;
; INPUTS:
;
;   nph:      number of count energy bins used
;   nee:      number of electron energy bins used
;   ee:       electron energy vector
;   regsol:   recovered electron visibility array
;   err:      uncertainty in visibility value
;   gs:       count visibility spectrum (cnts/cm^2/s/kev at detector)
;   eps:      count energy vector
;   wei:      weight vector
;   confidencestrip:    number of solution realizations to determine. (default = 10)
;   u:        singular vectors of the svd of the bremsstrahlung integral operator
;   w:        singular values of the svd of the bremsstrahlung integral operator
;   opt:      regularization parameter
;   sf:       singular functions of the svd of the bremsstrahlung integral operator
;   fitph:    rescaling of the count visibility spectrum
;   fitel:    rescaling of the electron visibility spectrum
;   Z:        value of the root-mean-square atomic number of the target. Default = 1.2
;
; OUTPUTS
;
;   DATASTRUCT:  a structure storing 2 arrays
;                   - Tag 1 :     STRIP       DOUBLE    Array[NEE, N+2]   where NEE = number of electron energy points
;                                                                                 N = number of solution realizations 
;                                                                                     (defined by confidencestrip value)
;              
;                     The array contains the electron energies and each realization of the regularized solution
;                     in the following format:;              
;                        datastruct.strip[*,0]     =  a vector of NEE elements containing the electron energies
;                        datastruct.strip[*,1]     =  a vector of NEE elements containing the original (unperturbed) regularized
;                                                     solution at the specified energies
;                        datastruct.strip[*,2:N+1] =  array containing N different regularized electron visibilities, at the
;                                                     specified energies, produced by randomly perturbing the input data
;              
;              
;                   - Tag 2 : RESIDUALS       DOUBLE    Array[NPH, 6]   where NPH = number of count energy points
;              
;                     The data is in six columns, with each row arranged as follows:;              
;                        datastruct.residuals[*,0] = Count energy used
;                        datastruct.residuals[*,1] = Count visibility values (input data)
;                        datastruct.residuals[*,2] = Uncertainty in count visibility values (input data)
;                        datastruct.residuals[*,3] = Count visibilities corresponding to the recovered electron visibility array
;                        datastruct.residuals[*,4] = Residual count visibilities (actual - recovered, i.e., column 2 - column 4)
;                        datastruct.residuals[*,5] = Cumulative residual, defined as
;                                                     C_j = (1/j) sum_[i=1]^j res_i
;
;
; RESTRICTIONS:
;     - CHOICE OF THE REGULARIZATION PARAMETER:
;         The code uses the cumulative residuals to fix the optimal regularization parameter.
;         The cumulative residuals can be used as a measure of residual 'clustering'; for an
;         acceptable fit to the data, abs(C_j) should be less than n/sqrt(j), where n is the
;         number of allowed standard deviations (n=1 is used in the code to drive the solution
;         to an 'acceptable' solution).
;
;-


pro vis_reg_ge_confidence,nph,nee,ee,regsol,$
    err,gs,eps,wei,conf,u,w,opt,sf,drm,fitph,fitel,Z,datastruct

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

num=conf                ;;;;;;; Number of iterations for the strip of confidence
if (conf eq 1 ) then num=num-1

regsoliter=dblarr(max(nee),num+1)

regsoliter[*,0]=regsol  ;;;;;;; save original values for regsol and gs
gsorig=gs

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; Confidence strip
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


maxregsol=regsol
minregsol=regsol

if (conf ne 1) then begin

    S=SYSTIME(1)
    S=5.
    arg1=dblarr(nph,max(nee))

    for iter=0,num-1 do begin

        a=err*randomn(S, nph) ; Normally-distributed pseudo-random numbers with a mean of zero and
                              ; a standard deviation of sigma (err) ==> between -3*err and 3*err
                              ; with probability 99%

        ;;; Recalculate values out of 3*sigma;;

        for i=0,nph-1 do begin
          if (abs(a[i]) gt 3*err[i]) then begin
               while (abs(a[i]) gt 3*err[i]) do a[i]=err[i]*randomn(S)
          endif
        end

        gs=gsorig+a

        ;;;;;;; Compute the corresponding regularized solution by using the regularization
        ;;;;;;; parameter computed during the first iteration.

        gs=gs/fitph

        for k=0,nph-1 do begin
            scal=(wei*gs)##u[k,*]
            for j=0,max(nee)-1 do arg1[k,j]=scal*w[k]*sf[k,j]/(w[k]*w[k]+opt)
        end
        regsol=total(arg1,1)                    ;;;;;;; regularised solutions of the strip
        regsol=fitel*regsol  ;;;;;;; rescaling

        gs=gs*fitph

        ; Update maxregsol and minregsol values

        for j=0,max(nee)-1 do begin
            if(regsol[j] le minregsol[j]) then minregsol[j]=regsol[j]
            if(regsol[j] ge maxregsol[j]) then maxregsol[j]=regsol[j]
        end

        regsoliter[*,iter+1]=regsol ;;;;; At each iteration save the current regularised solution

    end

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; End of the confidence strip
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

regsol=regsoliter[*,0] ;;;;; initial solution
gs=gsorig

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;  Define datastruct structure and fill 1 (over 2) tags ;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

datastruct={strip:dblarr(max(nee),num+2),residuals:dblarr(nph,6)}

datastruct.strip[*,0]=ee[*]
datastruct.strip[*,1:num+1]=regsoliter[*,0:num]

end