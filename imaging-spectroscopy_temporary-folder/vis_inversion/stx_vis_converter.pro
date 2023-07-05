FUNCTION stx_vis_converter, energymin, energymax, u, v, realvis, imvis, err, time_interval, xyoffset, det_used, lab_detused

  n_vis       = n_elements(u)
  vis         = replicate(stx_visibility(),  n_vis)
  vis.isc     = det_used
  vis.label   = lab_detused
  ;live_time
  vis.energy_range  = [energymin, energymax]
  vis.time_range    = time_interval
  vis.obsvis  = COMPLEX(realvis,imvis)
  ;tot_counts
  ;tot_counts_bkg
  ;totflux
  vis.sigamp  = err
  vis.u       = u
  vis.v       = v
  ;phase_sense
  vis.xyoffset = xyoffset
  ;xy_flare  
 
  RETURN, vis
 
END