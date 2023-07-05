; **************************************************************
; Procedure that is meant to be used in stx_imaging_spectroscopy
; **************************************************************

pro plot_boxes, these_boxes, box_or_fix, color=this_color
  default,this_color,'slate gray'
  for i=0,n_elements(these_boxes[0,0,*])-1 do begin
    if box_or_fix[i] eq 0 then begin
      plots,[these_boxes[0,0,i],these_boxes[1,0,i]],[these_boxes[0,1,i],these_boxes[0,1,i]],/data,linest=2,color=cgcolor(this_color),thick=2
      plots,[these_boxes[0,0,i],these_boxes[1,0,i]],[these_boxes[1,1,i],these_boxes[1,1,i]],/data,linest=2,color=cgcolor(this_color),thick=2
      plots,[these_boxes[0,0,i],these_boxes[0,0,i]],[these_boxes[0,1,i],these_boxes[1,1,i]],/data,linest=2,color=cgcolor(this_color),thick=2
      plots,[these_boxes[1,0,i],these_boxes[1,0,i]],[these_boxes[0,1,i],these_boxes[1,1,i]],/data,linest=2,color=cgcolor(this_color),thick=2
    endif else begin
      plots,these_boxes[0,0,i],these_boxes[0,1,i],/data,psym=1,syms=3,color=cgcolor(this_color),thick=3
    endelse
  endfor
end    ; End of the plot_boxes procedure