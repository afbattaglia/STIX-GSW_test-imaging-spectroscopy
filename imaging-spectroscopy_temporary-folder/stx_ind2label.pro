;+
;
; NAME:
;   stx_ind2label
;
; PURPOSE:
;   Return an array of detector labels from the corresponding array of indices
;   Inverse function of stx_label2ind
;
; INPUTS:
;   input_indices: array containing detector indices to be converted to corresponding labels (number 10-1 (resolution) and letter a-c (orientation)))
;
; HISTORY: December 2023, Stiefel M.Z., initial release
;
;-

FUNCTION stx_ind2label, input_indices

  ref_label = ['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c','6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c','2a','2b','2c','1a','1b','1c']
  ref_ind = stx_label2ind(ref_label)
  
  label = strarr(n_elements(input_indices))
  for i=0,n_elements(input_indices)-1 do begin
    
    label[i] = ref_label[where(ref_ind EQ input_indices[i])]
    
  endfor

  return, label

END
