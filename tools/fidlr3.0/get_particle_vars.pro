function get_particle_vars, filename

; given a file, determine whether it is HDF5 or pnetcdf, and
; read in the list of variables in its particle data structure.

if (n_elements(filename) EQ 0) then begin
  print, 'ERROR: no filename specified to get_var_list'
  return, -1
endif

; determine what type of file it is (hdf5 or pnetcdf)
file_type = determine_file_type(filename)
if (file_type EQ -1) then begin
  print, 'ERROR: file format of unknown type'
  return, -1
endif

; determine version of flash (2 or 3)
flash_version = determine_flash_version(filename)

; determine precision (single or double) of the file
file_precision = determine_file_precision(filename)

fail = 0

;------------------------------- from read_amr.pro
if (file_type EQ 1) then  begin ; HDF5
  file_identifier = H5F_OPEN(filename)

  dataset = H5D_OPEN(file_identifier, "particle names")
  particle_names =  H5D_READ(dataset)
  H5D_CLOSE, dataset
  H5F_CLOSE,file_identifier
endif else begin  ; NetCDF
  file_identifier = NCDF_OPEN(filename)

  gdata = ncdf_inquire(file_identifier)
  if (flash_version EQ 3) then begin 
    ; plot files and checkpoint files very
    ; similarly configured in flash3
    dimid = ncdf_dimid (file_identifier, 'dim_npart_props')
    ncdf_diminq, file_identifier, dimid, name, size
    particle_names = strarr(1,size)
    for i=0, size-1 do begin
      attname = 'particle_props_name_' + strtrim(string(i), 2)
      ncdf_attget, file_identifier, /global, attname, name
      particle_names[i] = string(name)
    endfor
  endif else begin ; end of flash 3 start of flash2
    if (file_precision eq 1) then begin ; plot files, single precision
      particle_vars = gdata.nvars - 6
      particle_names = strarr(1,particle_vars)
      for i=6, gdata.nvars-1 do begin
        particle_name = ncdf_varinq(file_identifier, i)
        particle_names[i-6] = particle_name.name
      endfor
    endif else begin        ; checkpoint files, double precision (same as flash3)
      ncdf_attget, file_identifier, /global, 'particle_names', particle_vars
      particle_names = strarr(1,particle_vars)
      for i=0, particle_vars-1 do begin
        attname = 'particle_' + strtrim(string(i), 2)
        ncdf_attget, file_identifier, /global, attname, name
        particle_names[i] = string(name)
      endfor
    endelse 
  endelse ; end of flash2
  NCDF_CLOSE,file_identifier
endelse ; end of hdf/cdf split for file open and particle_names read

; make it into a single vector
particle_names = reform(temporary(particle_names))

; get rid of whitespaces in variable names
; DEV no, dont do this now.... will have to be fixed when four space
; limitation is removed
;particle_names = strcompress(particle_names,/REMOVE_ALL)

return, particle_names

end
