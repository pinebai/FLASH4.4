function determine_geometry, filename

; returns the geometry type of a file
; the value is returned as an upper case string e.g. CARTESIAN as
; opposed to "cartesian" or "Cartesian"
; given a file, determine whether it is HDF or HDF5, and read in the 
; list of variables it stores appropriately

if (n_elements(filename) EQ 0) then begin
  print, 'ERROR: no filename specified to get_geometry'
  return, -1
endif

; get the xflash directory from the xflash_dir environmental variable
xflash_dir = get_xflash_path()

; determine what type of file it is
file_type = determine_file_type(filename)

; determine flash version
flash_version = determine_flash_version(filename)

file_version = determine_file_version(filename)

; determine precision
precision = determine_file_precision(filename)

if (file_type EQ -1) then begin
  print, 'ERROR: file is of unknown format'
  return, -1
endif else if (file_type eq 1) then begin ; hdf5
  file_identifier = H5F_OPEN(filename)

  if (file_version LE 7) then begin
    group = H5G_OPEN(file_identifier, "/")
    attribute = H5A_OPEN_NAME(group, "geometry name")

    geometry = H5A_READ(attribute)
    
    H5A_CLOSE, attribute
    H5G_CLOSE, group

  endif else begin
    dataset = H5D_OPEN(file_identifier, "string scalars")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=string_scalars)

    string_scalars = H5D_READ(dataset)
    H5D_CLOSE, dataset
    H5T_CLOSE, datatype

    for i=1, (size(string_scalars))[3] do begin
      if (stregex(string_scalars[i-1].name, '^geometry', /BOOLEAN)) then geometry = string_scalars[i-1].value
    endfor
    endelse
        
    H5F_CLOSE, file_identifier
    
  endif else if (file_type eq 2) then begin ; netCDF
    file_identifier = NCDF_OPEN(filename)

    if (flash_version EQ 2) then begin ; flash version 2
      ncdf_attget,file_identifier,/global,"geometry_name",geom
    endif else begin ; flash version 3
      ncdf_attget,file_identifier,1,"geometry",geom
    endelse 

    geometry = string(geom)
    NCDF_CLOSE, file_identifier
    
endif

; set uppercase to eliminate confusion between various versions
geometry = STRUPCASE(geometry)

return, geometry
end 
