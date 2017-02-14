function determine_file_dimensionality, filename

; given a filename, return the number of dimensions of the FLASH
; dataset.  Right now, since there is no explicit ndim record in the
; files, do this based on nxb, nyb, and nzb

if (n_elements(filename) EQ 0) then begin
  print, 'ERROR: no filename specified to get_var_list'
  return, -1
endif

; get the xflash directory from the xflash_dir environmental variable
xflash_dir = get_xflash_path()


; determine what type of file it is
;
; itype = 1 is HDF5
;       = 2 is NETCDF
;       = -1 is unknown or file does not exist
itype = determine_file_type(filename)

if (itype EQ -1) then begin
  print, 'ERROR: file is of unrecognized type or does not exist'
  return, -1
endif

if (itype EQ 1) then begin
  file_identifier = H5F_OPEN(filename)
 
  ;Old way of doing things.  Now base this off of GID
  ;dataset = H5D_OPEN(file_identifier, "coordinates")
  ;dataspace = H5D_GET_SPACE(dataset)
  ;dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)

  ;ndim = dims[0]
  
  dataset = H5D_OPEN(file_identifier, "gid")
  dataspace = H5D_GET_SPACE(dataset)
  dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)
  
  switch dims[0] of
    5: begin
       ndim = 1
       break
    end
    9: begin
       ndim = 2
       break
    end
    15: begin
        ndim  = 3
        break
    end
    else: message, 'GID size per block does not match 1, 2 or 3d values. Read will Fail.'
    endswitch

  H5D_CLOSE,dataset
  H5S_CLOSE, dataspace

  H5F_CLOSE, file_identifier
  return, ndim
endif

if (itype EQ 2) then begin
  file_identifier = NCDF_OPEN(filename, /NOWRITE)
  dimid = ncdf_dimid(file_identifier, 'dim_NDIM')
  ncdf_diminq, file_identifier, dimid, nndim, ndim
  NCDF_CLOSE, file_identifier

  return, ndim
endif

end
