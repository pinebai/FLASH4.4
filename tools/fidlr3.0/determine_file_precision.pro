function determine_file_precision, filename
; determines the file precision
;  RETURNS precision = 1 for single precision
;          precision = 2 for double precision
;          precision = -1 for an error in reading
;   
if (n_elements(verbose) eq 0) then verbose = 0

; first, make sure the file is found
file_exists = file_test(filename)
if (not file_exists) then begin
  print,'ERROR:  in determine_file_precision, file ',filename,' does not exist.'
  return, -1
endif 

; second, make sure the file is HDF5 or NetCDF
file_type = determine_file_type(filename)
case file_type of 
  1: HDF5 = 1
  2: HDF5 = 0
  -1: begin 
    print, 'ERROR:  Problem determining filetype of file ',filename
    return, -1
  end 
endcase 
if (verbose) then print, ' File Type is ',file_type


; find out what software version we're running
flashVersion = determine_flash_version(filename)
if (verbose) then print,' FLASH Version is ', flashVersion


if (flashVersion eq 3) then begin          ; written into files, can determine
  if (HDF5) then begin 
    fileid = H5F_OPEN(filename)
    dataset = H5D_OPEN(fileid,"logical scalars")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=logical_scalars)
    logical_scalars = H5D_READ(dataset)
    H5D_CLOSE, dataset
    H5T_CLOSE, datatype
    H5F_CLOSE, fileid
    for i=1, (size(logical_scalars))[3] do begin 
      if (stregex(logical_scalars[i-1].name, '^corners', /BOOLEAN)) then $
        corners_bool = logical_scalars[i-1].value
      if (stregex(logical_scalars[i-1].name, '^double_precision', /BOOLEAN)) then $
        double_bool = logical_scalars[i-1].value

    endfor ; end of loop over size of logical_scalars

  endif else begin ; end of HDF, begin of NCDF

    fileid = NCDF_OPEN(filename)
    ; get the right variable id, note that corners is stored in BOTH
    ;  "runtime_parameters" and "scalars"
    varid = ncdf_varid(fileid,"scalars")
    ncdf_attget,fileid,varid,"corners",corners_bool
    ncdf_attget,fileid,varid,"double_precision",double_bool
    NCDF_CLOSE,fileid
    ; flash 2 syntax
    ; ncdf_attget,fileid,/GLOBAL,"corners",corners_bool
    ; ncdf_attget, fileid, /GLOBAL, "doubleprecision",double_bool
  endelse                         ; end of pnetCDF FLASH3 

  ; info is found, now set parameters
  if (corners_bool) then begin
    corners = "yes"
  endif else begin
    corners = "no"
  endelse
  if (double_bool) then begin
    precision = 2
  endif else begin
    precision = 1
  endelse        

endif else begin ; end of flash3, begin flash 2

  ; information not stored in flash2 files, 
  usedHack = 0
  if (strpos(filename, 'plt') GT 0) then begin ; plot files, single precision
    precision = 1
    usedHack = 1
  endif else if (strpos(filename,'chk') GT 0) then begin ; checkpoint files, double precision
    precision = 2
    usedHack = 1
  endif
  ; if (usedHack ) then print, "Used a filename hack to determine that precision is ",precision
  if (not usedHack) then begin ; no clues from the filename, try to determine
    print, "Sorry, no clues in the FLASH2 filename; assuming double precision"
    precision =  2
  endif  

endelse  ; end of FLASH2

return, precision

end   ; of determine_file_precision.pro    
