
ADD_PATH, '/home/afbattaglia/Documents/ETHZ/PhD/Codes/IDL/electron_maps/STIXvisibilityInversion\'
setenv,'STX_DEMO_DATA=/home/afbattaglia/Software/idl_lib/STIX-GSW-grid-2-imaging/STIX-GSW/stix/idl/demo'

;;************************************************ 7-May-2021 18:42:38 ************************************************
; Folder in which the files downloaded for this demonstration are stored
out_dir = concat_dir( getenv('STX_DEMO_DATA'),'imaging', /d)
; URL of the STIX data center
website_url = 'https://datacenter.stix.i4ds.net/download/fits/bsd/'
; UID of the science fits file to be dowloaded from the website
uid_sci_file = '2105070034'
; Download the science fits file (if not already stored in out_dir)
;sock_copy, website_url + uid_sci_file, out_name, status = status, out_dir = out_dir, $
;  local_file=path_sci_file, clobber=0
time_range = ['7-May-2021 18:51:00', '7-May-2021 18:53:40']
flare_loc   = [425., 324.] ; Adapted from imaging
website_url = 'http://dataarchive.stix.i4ds.net/fits/L2/2021/05/07/AUX/'
file_name    = 'solo_L2_stix-aux-ephemeris_20210507_V01.fits'
sock_copy, website_url + file_name, out_name, status = status, out_dir = out_dir, $
  local_file=aux_fits_file, clobber=0
aux_data = stx_create_auxiliary_data(aux_fits_file, time_range)
; SolarOrbiter - Sun distance
stx_read_aux_fits,  aux_fits_file, aux_data=aux_data_str
dist_solo_sun =  average(aux_data_str.SOLO_LOC_CARRINGTON_DIST)*10.^(5)


path_sci_file = '/home/afbattaglia/Software/idl_lib/STIX-GSW-grid-2-imaging/STIX-GSW/stix/idl/demo/imaging/solo_L1A_stix-sci-xray-l1-2105070039_20210507T185059-20210507T185359_031719_V01.fits'


subc_index = stx_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c',$
                             '6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c', $
                             '2a','2b','2c'])

lower_energy_edge = [8.,10.,12.,14.,16.,18.,20.,22.]
upper_energy_edge = [10.,12.,14.,16.,18.,20.,22.,25.]
mapcenter_stix = stx_hpc2stx_coord(flare_loc, aux_data)
xy_flare_stix  = mapcenter_stix

vis_tot = [] ; array of visibilities

for i=0,n_elements(lower_energy_edge)-1 do begin

  energy_range = [lower_energy_edge[i],upper_energy_edge[i]]

  this_estring_0=strtrim(fix(energy_range[0]),2)
  this_estring_1=strtrim(fix(energy_range[1]),2)

  vis = stx_construct_calibrated_visibility(path_sci_file, time_range, energy_range, mapcenter_stix, subc_index=subc_index, $
                                            path_bkg_file=path_bkg_file, xy_flare=xy_flare_stix)

  vis_tot = [vis_tot, vis]
  
end
stop
e_lower  = vis_tot.energy_range[0]
e_high = vis_tot.energy_range[1]
e_lower1 = e_lower[UNIQ(e_lower, SORT(e_lower))]
e_high1 =  e_high[UNIQ( e_high, SORT( e_high))]
print, e_lower1
print, e_high1

;******************************** INVERSION *************************************
attenuator = 0.
stx_visibilities_inversion,vis_tot, dist_solo_sun, attenuator, reg_el_vis, orig_ph_vis, reg_ph_vis

e_lower_e  = reg_el_vis.energy_range[0]
e_high_e = reg_el_vis.energy_range[1]
e_lower_e1 = e_lower_e[UNIQ(e_lower_e, SORT(e_lower_e))]
e_high_e1 =  e_high_e[UNIQ( e_high_e, SORT( e_high_e))]
print, e_lower_e1
print, e_high_e1


; ******************** Reconstruction from photon visibilities **************************************
loadct,5
window, 5, xsize =2500., ysize =1500.
!p.multi = [0,4,2]

for jj = 0, n_elements(lower_energy_edge)-1 do begin
  loc     = where(vis_tot.energy_range[0] eq lower_energy_edge[jj], count)
  vis_ph = replicate(stx_visibility(),  count)
  vis_ph = vis_tot[loc]
  mem_ph = stx_mem_ge(vis_ph,[129, 129], [1.,1.],aux_data, total_flux = max(abs(vis_ph.obsvis)))
  this_estring_0=strtrim(fix(lower_energy_edge[jj]),2)
  this_estring_1=strtrim(fix(upper_energy_edge[jj]),2)
  plot_map, mem_ph, /limb, title =  this_estring_0+'-'+this_estring_1+'keV', charsize=1.5
endfor
stop

; ************************* reg_ph **************************************
loadct,5
window, 6, xsize =2500., ysize =1500.
!p.multi = [0,4,2]

for jj = 0, n_elements(lower_energy_edge)-1 do begin
  loc     = where(reg_ph_vis.energy_range[0] eq lower_energy_edge[jj], count)
  vis_reg_p = replicate(stx_visibility(),  count)
  vis_reg_p = reg_ph_vis[loc]
  mem_reg_p = stx_mem_ge(vis_reg_p,[129, 129], [1.,1.],aux_data, total_flux = max(abs(vis_reg_p.obsvis)))
  this_estring_0=strtrim(fix(lower_energy_edge[jj]),2)
  this_estring_1=strtrim(fix(upper_energy_edge[jj]),2)
  plot_map, mem_reg_p, /limb, title =  this_estring_0+'-'+this_estring_1+'keV', charsize=1.5
endfor
stop

; ************************* reg_e **************************************
loadct,5
window, 7, xsize =2500., ysize =1500.
!p.multi = [0,5,3]

for jj = 0, n_elements(e_lower_e1)-1 do begin
  loc     = where(reg_el_vis.energy_range[0] eq e_lower_e1[jj], count)
  vis_e = replicate(stx_visibility(),  count)
  vis_e = reg_el_vis[loc]
  mem_e = stx_mem_ge(vis_e,[129, 129], [1.,1.],aux_data, total_flux = max(abs(vis_e.obsvis)))
  this_estring_0=strtrim(fix(e_lower_e1[jj]),2)
  this_estring_1=strtrim(fix(e_high_e1[jj]),2)
  plot_map, mem_e, /limb, title =  this_estring_0+'-'+this_estring_1+'keV', charsize=1.5
endfor

end
