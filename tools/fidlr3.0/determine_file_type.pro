function determine_file_type, filename
;
; determine whether a file is an HDF 5 file, a netCDF
; file, neither, or does not exist
; now uses IDL routines and does not depend upon file naming convention
; now uses IDL routines and does not depend upon file naming convention
;
; return: -1 if file does not exist or is an unknown format
;          1 if HDF 5
;          2 if NetCDF
;

; first determine if the file exists!
file_exists = file_test(filename)
IF (NOT file_exists) THEN BEGIN
  print,'ERROR: file ',filename,' does not exist.'
  return, -1
ENDIF 

; initialize the file type to unknown
file_type = -1

is_HDF5 = H5F_IS_HDF5(filename)
if (is_HDF5 EQ 1) then begin
    file_type = 1
    return, file_type
endif

is_NCDF = NCDF_IS_NETCDF(filename)
if (is_NCDF GT 0) then begin
    file_type = 2
    return, file_type
endif

return, file_type

END ; determine_file_type.pro 

