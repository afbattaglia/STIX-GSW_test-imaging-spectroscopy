;+
; :description:
;    This procedure reads in all of the imaging spectroscopy map save files produced by stx_imaging_spectroscopy.pro in a given folder. 
;    These are then written to a FITS image cube FITS file. This file is then passed to an OSPEX object and the region of interest selection 
;    widget opened. 
;     
;
; :categories:
;    imaging spectroscopy
;
; :params:
;    folder : in, required, type="string"
;             the path to the IDL save files containing the maps to be passed to OSEPX
;             filenames should be of the format 'stix-imaging-spectroscopy*.sav'
;
;    algo :   in, type="string", default: 'clean'
;             Extract the maps from the given imaging algorithm
;             currently supported are: 'clean', 'fwdfit' and 'memge'
;
; :examples:
;    folder = /STIX-GSW_test-imaging-spectroscopy/2023-05-16_172029-172149_SolarOrbiter-UT_test6_obs-vis_two-sources_min-FWHM
;    stx_img_spectra_sav2ospex, folder,'clean'
;
; :history:
;    30-Aug-2023 - ECMD (Graz), initial release
;
;-
pro stx_img_spectra_sav2ospex, folder, algo

  default, algo, 'clean'

  cd, folder
  flspec = findfile('stix-imaging-spectroscopy*.sav')
  l = list()

  for i = 0, n_elements(flspec)-1 do begin
    restore, flspec[i]
    case algo of
      'clean': image = { $
        type        : "stx_image", $
        algo        : 'clean', $
        time_range  : CLEAN_MAP[0].time_range, $
        energy_range: CLEAN_MAP[0].energy_range, $
        map : CLEAN_MAP[0] $
      }

      'fwdfit': image = { $
        type        : "stx_image", $
        algo        : 'fwd', $
        time_range  : FWDFIT_MAP.time_range, $
        energy_range: FWDFIT_MAP.energy_range, $
        map : FWDFIT_MAP $
      }

      'memge' : image = { $
        type        : "stx_image", $
        algo        : 'memge', $
        time_range  : MEMGE_MAP.time_range, $
        energy_range: MEMGE_MAP.energy_range, $
        map : MEMGE_MAP $
      }
      else: message, 'Algorithm '+algo+' not found. Options are: clean, fwdfit or memge'
    endcase
    l.add,image

  endfor
  s = l.toarray()

  filename = 'stx_image_cube_' +algo+'_'+time2fid(atime(TIME_RANGE_SO[0]),/time,/sec,/full)+'-'+ strmid(time2fid(atime(TIME_RANGE_SO[1]),/time,/sec),7,6)+'.fits'
  
  stx_map2fits_cube, s, filename, maps, energy_axis = energy_axis, time_axis = time_axis

  ospex_obj = ospex()
  ospex_obj->set, spex_specfile = filename
  ospex_obj->set, spex_roi_use = 0
  ospex_obj-> roi


end