;+
;
; NAME:
;
;   stx_imaspec_adapt_map
;
; PURPOSE:
;
;   Adapt a already existing stix map strucutre by changing the amount of saved visibilities
;   Routine developed for stx_imaging_spectroscopy.pro
;
; CALLING SEQUENCE:
;
;   stx_map = stx_imaspec_adapt_map(im_map, vis, pixel)
;
; INPUTS:
;
;   im_map: existing stix map
;
;   vis: new visibility structure, which should be used for the map structure
;   
;   pixel: bi-dimensional array containing the pixel size in arcsec
;
; OUTPUTS:
;
;   Map structure, the same as when using stx_make_map.pro routine
;
; HISTORY: December 2023, Stiefel M.Z., created
;
; CONTACT:
;   muriel.stiefel@fhnw.ch
;-

function stx_imaspec_adapt_map, im_map, vis, pixel

  imsize = size(im_map.data, /dim)
  F = vis_map2vis_matrix(vis.u, vis.v, imsize, pixel)
  mapvis = F # im_map.data[*]

  pred_vis = {mapvis: mapvis, $     ; Complex visibility values predicted from the map
    isc:    vis.ISC, $    ; Subcollimator index
    label:  vis.LABEL, $  ; Subcollimator label
    u: vis.u, $           ; U coordinates of the frequencies sampled by the subcollimators
    v: vis.v  $           ; V coordinates of the frequencies sampled by the subcollimators
  }


  out_map = make_map(im_map.data)

  out_map.ID = im_map.id
  out_map.dx = im_map.dx
  out_map.dy = im_map.dy
  out_map.time = im_map.time
  out_map.dur  = im_map.dur
  out_map.xc = im_map.xc
  out_map.yc = im_map.yc
  out_map.roll_angle = im_map.roll_angle

  ;; Add properties
  energy_range   = im_map.ENERGY_RANGE

  add_prop, out_map, energy_range   = im_map.ENERGY_RANGE 
  add_prop, out_map, obs_vis        = vis                   
  add_prop, out_map, pred_vis       = pred_vis         
  add_prop, out_map, aux_data       = im_map.aux_data              
  add_prop, out_map, time_range     = im_map.time_range            
  add_prop, out_map, rsun = im_map.RSUN
  add_prop, out_map, b0   = im_map.B0
  add_prop, out_map, l0   = im_map.L0
  add_prop, out_map, coord_frame = 'SOLO HPC'
  add_prop, out_map, units = 'counts cm^-2 asec^-2 s^-1 keV^-1'

  return, out_map


end