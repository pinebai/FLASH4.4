FUNCTION determine_file_version, filename
; reads a file and determines the internal format version
;  RETURNS 
;       an longword integer indicating the format version
;       -1 if file not found or not a HDF5 file
;  NOTE
;       FLASH2.5 used format 7
;       FLASH3 uses format   8 (for later files)
;   

format_version = -1

; first, make sure the file is found
file_exists = file_test(filename)
IF (NOT file_exists) THEN BEGIN
  print,'ERROR:  in determine_file_version, file ',filename,' does not exist.'
  return, -1
ENDIF 

; second, make sure the file is HDF5 or NetCDF
file_type = determine_file_type(filename)

IF (file_type LE 0) THEN  BEGIN
  print,'ERROR: in determine_file_version, file ',filename,' is not an known format.' 
  return, -1
ENDIF 

; now open and read the formatVersion
IF (file_type EQ 1) THEN BEGIN  ; HDF5
    file_identifier = H5F_OPEN(filename)
    IF (file_identifier LE 0) THEN BEGIN
        print,'ERROR: in determine_file_version, file ',filename,' could not be opened.'
        return, -1
    ENDIF 

    ; in flash3, we've moved the file format version into the sim info dataset
    ; in flash2, it was it's own dataset
    file_contents = h5_parse(filename)
    tag_contents = tag_names(file_contents)
    tag_where = where(tag_contents EQ "FILE_FORMAT_VERSION")
    ; NOTE that tag_names removes white space and capitalizes everything!
    IF ( tag_where EQ -1) then BEGIN
       dataset = H5D_OPEN(file_identifier,"sim info")
       datatype = H5D_GET_TYPE(dataset)
       idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_info) ;; this line is not strictly necessary
       sim_info = H5D_READ(dataset)
       format_version = strtrim(sim_info.file_format_version, 2)
       H5T_CLOSE, datatype
       H5D_CLOSE, dataset
    ENDIF ELSE BEGIN  ; previous to flash3
       dataset = H5D_OPEN(file_identifier,"file format version")
       format_version = H5D_READ(dataset)
       H5D_CLOSE,dataset
    ENDELSE
    H5F_CLOSE,file_identifier
ENDIF ; end of HDF5

IF (file_type EQ 2) THEN BEGIN ; netCDF
    file_identifier = NCDF_OPEN(filename)
    IF (file_identifier LT 0) THEN BEGIN
        print,'ERROR: in determine_file_version, file ',filename,' could not be opened.'
        return, -1
    ENDIF 
;   see if you can find "file format version" as any of the variables,
;     or attributes.  Draws heavily on Gumley routines, p. 177
    attnames = ncdf_attdir(file_identifier,'');
    where_format = where(attnames EQ "file_format_version")

    IF (where_format GE 0) THEN BEGIN ; works for flash3
        ncdf_attget,file_identifier,/GLOBAL,"file_format_version",format_version
    ENDIF ELSE BEGIN ; suspected flash2
; DEV hack here -- need to check against a file
        print, '  NOTICE could not find file_format_version in file ',filename
        print, '  I will guess that your file is FLASH2 format'
        format_version = 7
    ENDELSE 
    
    NCDF_CLOSE,file_identifier

ENDIF ; end of NCDF

return, format_version

END 
