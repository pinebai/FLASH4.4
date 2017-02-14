function get_particle_number, filename

; check for the presence of particle data in the file and return the 
; number of particles stored in the file if found.

if (n_elements(filename) EQ 0) then begin
  print, 'ERROR: no filename specified to get_var_list'
  return, -1
endif

; get the xflash directory from the xflash_dir environmental variable
xflash_dir = get_xflash_path()


; determine what type of file it is
itype = determine_file_type(filename)

if (itype EQ -1) then begin
  print, 'ERROR: file does not exist'
  return, -1
endif

; determine flash2 or flash3
iflash = determine_flash_version(filename)

numParticles = 0l

; if hdf5 file
if(itype EQ 1) then begin ; hdf5

  ; sometimes we don't store the particle information -- check first
  file_contents = h5_parse(filename)
  ; flash 2 naming?  or does h5_parse make all upcase and single word?
  if ((where(tag_names(file_contents) EQ "TRACER_PARTICLES"))[0] NE -1) then begin

    file_identifier = H5F_OPEN(filename)
            
    dataset = H5D_OPEN(file_identifier, "tracer particles")
    dataspace = H5D_GET_SPACE(dataset)

    dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)

    ; DEV lbr suspicious change here, but who knows how katy reorganized....
    if (iflash EQ 2) then  numParticles = dims[0]
    if (iflash EQ 3) then  numParticles = dims[1]
            
    H5D_CLOSE, dataset
    H5S_CLOSE, dataspace
    H5F_CLOSE, file_identifier

  endif else begin

    numParticles = 0

  endelse
endif ; end of hdf5 file

; if pnetcdf file
if(itype EQ 2) then begin
  file_identifier = NCDF_OPEN(filename)
  dimid = ncdf_dimid(file_identifier, 'dim_particles')
  if (dimid GE 0) then begin
    ncdf_diminq, file_identifier, dimid, name, size
    numParticles = size
  endif else begin
    numParticles = 0l
  endelse 
  NCDF_CLOSE, file_identifier
endif
      
; print, 'In get_particle_number, numParticles is ',numParticles      
return, numParticles
end
