;==============================================================================
;
; read a FLASH HDF5 file into IDL.  
;
; This routine reads in FLASH data in the HDF 5 format.  
;
; This version is for IDL >= 5.6, as it uses the built in HDF5
; support.
;
; Arguments:  filename -- the name of the file `to read in
;
;             var_name -- a 4-character string specifying the variable
;                         name to read in.  If this keyword is not
;                         present, all variables are read in.
;
;             /VERBOSE  -- output diagnostics while running
;
;==============================================================================

pro read_amr, filename, $
              VAR_NAME=var_name, $
              VERBOSE=verbose, $ 
              TREE=tree, $
              DATA=unk, $
              PARAMETERS=params, $
              STORED_VARS=unk_names, $
              NUM_PARTICLES=numParticle, $
              PARTICLES=particles, $
              INT_PROP_NAMES=IntPropNames, $
              REAL_PROP_NAMES=RealPropNames, $
              GEOMETRY=geometry

; clean up any pre-existing arrays
if n_elements(unk) GT 0 then undefine, unk
if n_elements(particles) gt 0 then undefine, particles
if (not Keyword_set(verbose)) then verbose = 0

; find out what software version we're running
flash_version = determine_flash_version(filename)

; There are no "logical" variables in IDL -- zero = false
if (flash_version EQ 2) then begin 
  FLASH3 = 0
  FLASH2 = 1
endif else begin 
  FLASH3 = 1
  FLASH2 = 0
endelse 

if Keyword_Set(verbose) then print,' FLASH Version is ', flash_version

;  determine the filetype of HDF5(1) or NetCDF(2)
file_type = determine_file_type(filename)

if (file_type EQ -1) then begin
  print, 'ERROR: file does not exist'
  return
endif

if Keyword_Set(verbose) then print, ' File Type is ',file_type

; DOUBLE is a flag to read data in double precision (checkpoint files)
; try to determine the precision from the file
file_precision = determine_file_precision(filename)
if (file_precision eq 1) then double = 0
if (file_precision eq 2) then double = 1

; check to see which variable to read; if var is not defined, read them all
if  n_elements(var_name) EQ 0 then begin
  var_name = 'none'
endif else begin 
  ; check that var name has four characters -- append with blanks if
  ; necessary.  format string pads with blanks on right if too short
  var_name = string(var_name, FORMAT='(A-4)')
endelse 

if Keyword_Set(verbose) then print, ' Variable var_name requested = ', var_name

;-----------------------------------------------------
;  begin reading
;---------------------------------------------------

;--------------------------------------------------------------------
; get the list of variables and make sure that our requested variable
; (if any) exists in the dataset
;-------------------------------------------------------------------
unk_names = get_var_list(filename)

;------------------------------------------------------------------
;  open the file
;------------------------------------------------------------------
if (file_type EQ 1) then file_identifier = H5F_OPEN(filename)
if (file_type EQ 2) then file_identifier = NCDF_OPEN(filename)

;------------------------------------------------------------------
; we can also now set the number of variables
;------------------------------------------------------------------
nvar = (size(unk_names))[1]
if (var_name NE 'none') then begin
  var = (where(unk_names EQ var_name))[0]
  if (var EQ -1) then begin
    print, 'ERROR: requested variable ',var_name,' not found in dataset'
    print, '       therefore reading in all variables'
  endif else begin 
    print, ' Reading in only the variable ', var_name
  endelse 
endif else begin
  var = -1
endelse

;-------------------------------------------------------------------
; read in the simulation parameters
;-------------------------------------------------------------------
if (file_type EQ 1) then begin 
  if(FLASH2) then begin

    dataset = H5D_OPEN(file_identifier, "simulation parameters")
    datatype = H5D_GET_TYPE(dataset)

    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_params)
    sim_params = H5D_READ(dataset)
    H5D_CLOSE, dataset
    H5T_CLOSE, datatype
  endif

  ; read in the integer scalars
  if (FLASH3) then begin
    dataset = H5D_OPEN(file_identifier, "integer scalars")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=int_scalars_list)

    ; int_scalars_list is a linked list which we will change
    ; into a struct whose elements and values are the tag and
    ; and value names of the linked list
    int_scalars_list = H5D_READ(dataset)
    for i=1, (size(int_scalars_list))[3] do begin
      if ((size(int_scalars))[0] eq 0) then begin  ; if int_scalars doesn't exist yet, create it
        int_scalars = create_struct(strtrim(int_scalars_list[i-1].name), int_scalars_list[i-1].value)
      endif else begin  ; otherwise, append to it
        int_scalars = create_struct(int_scalars, strtrim(int_scalars_list[i-1].name), int_scalars_list[i-1].value)
      endelse
    endfor

    H5T_CLOSE, datatype
    H5D_CLOSE, dataset  ; "integer scalars"

    ; read in the real scalars
    dataset = H5D_OPEN(file_identifier, "real scalars")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=real_scalars_list)

    ; real_scalars_list is a linked list which we will change
    ; into a struct whose elements and values are the tag and
    ; and value names of the linked list
    real_scalars_list = H5D_READ(dataset)
    for  i=1, (size(real_scalars_list))[3] do  begin 
      if ((size(real_scalars))[0] EQ  0) then  begin ; if real_scalars doesn't exist yet, create it
        real_scalars = create_struct(strtrim(real_scalars_list[i-1].name), real_scalars_list[i-1].value)
      endif else begin  ; otherwise, append to it`
        real_scalars = create_struct(real_scalars, strtrim(real_scalars_list[i-1].name), real_scalars_list[i-1].value)
      endelse 
    endfor 
    H5T_CLOSE, datatype
    H5D_CLOSE, dataset  ; "real scalars"
  endif                 ; end of readin for FLASH3
  ; DEV don't know where these are read in FLASH2

  ; get the dimensionality
  ; NOTE: below is the old way of getting ndim based on coordinate data
  ;       new way is based off of gid.
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
    else: message, 'gid size per block does not match 1, 2 or 3d values. Read will Fail.'
    endswitch


  H5S_CLOSE, dataspace
  H5D_CLOSE,dataset
endif else begin  ; end of hdf5, start of netCDF

  ; read in the simulation parameters
  if (FLASH3) then begin
    ncdf_attget, file_identifier, 1, "globalnumblocks",globalNumBlocks
    ncdf_attget, file_identifier, 1, "time", time
    ncdf_attget, file_identifier, 1, "dt", dt
    ncdf_attget, file_identifier, 1, "nstep", nsteps 
    ncdf_attget, file_identifier, 1, "nxb" , nxb
    ncdf_attget, file_identifier, 1, "nyb" , nyb
    ncdf_attget, file_identifier, 1, "nzb" , nzb

    int_scalars = create_struct('globalNumBlocks' , globalNumBlocks, $
                                'nsteps', nsteps, 'redshift', 1.0, $		
                                'nxb', nxb, 'nyb', nyb, 'nzb', nzb)
    real_scalars = create_struct('time',time,'dt',dt,'timestep',1.0)
  endif else begin  ; FLASH2 file
                    ;DEV: kda - this is the old format way
                    ;ncdf_attget, file_identifier, /GLOBAL, "timestep", timestep
                    ;ncdf_attget, file_identifier, /GLOBAL, "redshift", redshift 
    ncdf_attget, file_identifier, /GLOBAL, "total_blocks",total_blocks
    ncdf_attget, file_identifier, /GLOBAL, "time", time
    ncdf_attget, file_identifier, /GLOBAL, "timestep", timestep
    ncdf_attget, file_identifier, /GLOBAL, "nsteps", nsteps
    ncdf_attget, file_identifier, /GLOBAL, "redshift", redshift
    ncdf_attget, file_identifier, /GLOBAL, "nxb" , nxb
    ncdf_attget, file_identifier, /GLOBAL, "nyb" , nyb
    ncdf_attget, file_identifier, /GLOBAL, "nzb" , nzb

    sim_params = create_struct('total_blocks' , total_blocks, $
                               'nsteps', nsteps, 'redshift', redshift, $		
                               'time', time, 'timestep', timestep, $
                               'nxb', nxb, 'nyb', nyb, 'nzb', nzb)
  endelse  ; end of flash2/3 if block for netcdf 

  ; get the dimensionality
  dimid = ncdf_dimid(file_identifier, 'dim_NDIM')
  ncdf_diminq, file_identifier, dimid, nndim, ndim
endelse  ; end of netcdf

; Common -- get dimensions
nfaces = 2*ndim
nchild = 2^ndim

;---------------------------------------------------------
; figure out if we are dealing with corners
; DEV a hack in FLASH2, FLASH3 should write them
;---------------------------------------------------------
corners = 0
if (strpos(filename, 'crn_') GT 0) then begin
  corners = 1
endif

;------------------------------------------------------
;  get GEOMETRY
;-------------------------------------------------------

if (n_elements(geometry) eq 0) then begin
  geometry = determine_geometry(filename)
endif 

;------------------------------------------------------------------------------
;   setup the structures to pass the data
;----------------------------------------------------------------------------

; simulation parameters and int_scalars are different in flash2/3
if (FLASH3) then begin
  totBlocks = int_scalars.globalNumBlocks
  nxb = int_scalars.nxb
  nyb = int_scalars.nyb
  nzb = int_scalars.nzb
  time = real_scalars.time  ; double precision
  dt = real_scalars.dt      ; double precision
endif else begin            ; end flash3, begin FLASH2
  totBlocks = sim_params.total_blocks
  nxb = sim_params.nxb
  nyb = sim_params.nyb
  nzb = sim_params.nzb
  time = sim_params.time    ; double precision
  dt = sim_params.timestep  ; double precision
endelse 

; call a standard routine to set up tree structures
define_tree,file_precision,nfaces,nchild,ndim, TREE=tree
tree = replicate(tree,totBlocks)

; set up params structure to hold everything
if (not double) then begin  ; single precision, plot files
  params = {totBlocks:totBlocks, $
            corners:corners, $
            ndim:ndim, $
            nvar:nvar, $
            nxb:nxb, $
            nyb:nyb, $
            nzb:nzb, $
            ntopx:1, $
            ntopy:1, $
            ntopz:1, $
            time:float(time), $
            dt: float(dt), $
            redshift:1.0, $
            geometry:"unknown"}
endif else begin   ; double precision, checkpoint files
  params = {totBlocks:totBlocks, $
            corners:corners, $
            ndim:ndim, $
            nvar:nvar, $
            nxb:nxb, $
            nyb:nyb, $
            nzb:nzb, $
            ntopx:1, $
            ntopy:1, $
            ntopz:1, $
            time:time, $
            dt:dt, $
            redshift:1.d0, $
            geometry:"unknown"}
endelse   ; end of double/single blocks
params.geometry = strcompress(geometry, /REMOVE_ALL)

; HACK FLASH2
; for some reason, the redshift field is not stored in the plotfiles
if (FLASH2) then begin 
  if ( (where(tag_names(sim_params) EQ "REDSHIFT"))[0] NE -1) then begin
    if (not double) then begin
      params.redshift = float(sim_params.redshift)
    endif else begin
      params.redshift = sim_params.redshift
    endelse
  endif                       
endif  ; end of flash2 redshift hack

; DEV -- where is redshift stored in FLASH3?
; DEV need to do param/string readin described by Katy

;----------------------------------------------------------------------------
; We now have the dimension, total number of blocks, and the number of
; zones per block in each direction.  We can use this to read the tree
; information, coordinate information, and unknowns
;-----------------------------------------------------------------------------

;----------------------------------------------------------------------------
; read in the tree information
;----------------------------------------------------------------------------
ngid = long(2^ndim + 1 + 2*ndim)
if (file_type EQ 1) then begin 
  dataset = H5D_OPEN(file_identifier, "refine level")
  refine_level = H5D_READ(dataset)
  H5D_CLOSE, dataset

  dataset = H5D_OPEN(file_identifier, "node type")
  node_type = H5D_READ(dataset)
  H5D_CLOSE, dataset

  dataset = H5D_OPEN(file_identifier, "gid")
  gid = H5D_READ(dataset)
  H5D_CLOSE, dataset

  file_contents = h5_parse(filename)

  if ( (where(tag_names(file_contents) EQ "PROCESSOR_NUMBER"))[0] NE -1) then begin
    dataset = H5D_OPEN(file_identifier, "processor number")
    proc_number = H5D_READ(dataset)
    H5D_CLOSE, dataset
  endif else begin
    proc_number = 0
  endelse
endif else begin  ; end hdf5, begin netcdf
  ncdf_varget, file_identifier, "lrefine", refine_level
  ncdf_varget, file_identifier, "nodetype", node_type
  ncdf_varget, file_identifier, "processor_number", proc_number
  ncdf_varget, file_identifier, "gid", gid
endelse

tree[*].lrefine = refine_level
tree[*].nodeType = node_type
tree[*].processorNumber = proc_number

for i = 0, ngid-1 do begin
  tree[*].gid[i] = reform(gid[i,*])
endfor

undefine, refine_level
undefine, node_type
undefine, gid
undefine, proc_number

;----------------------------------------------------------------------------
; read in the coordinate infomation
;----------------------------------------------------------------------------
if (file_type EQ 1) then begin 
  ;in here we have to set up memspaces and dataspaces to get the right sizes
  ;of these datasets. 

  start_2d = [0,0]
  stride_2d = [1, 1]
  ; initialize to a big number so space is allocated for long integers.
  count_2d = [1,33000]
  count_2d[0] = ndim
  count_2d[1] = totBlocks

  start_3d = [0,0,0]
  stride_3d = [1,1,1]
  count_3d = [1,1,33000]
  count_3d[0] = 2
  count_3d[1] = ndim
  count_3d[2] = totBlocks

  ;***********coordinates**************
  dataset = H5D_OPEN(file_identifier, "coordinates")
  dataspace = H5D_GET_SPACE(dataset)
  H5S_SELECT_HYPERSLAB, dataspace, start_2d, count_2d, STRIDE=stride_2d, /RESET
  memspace = H5S_CREATE_SIMPLE(count_2d)
  coord = H5D_READ(dataset, FILE_SPACE=dataspace, MEMORY_SPACE=memspace)
  H5S_CLOSE, memspace
  H5D_CLOSE, dataset
  H5S_CLOSE, dataspace

  ; **********block size**************
  dataset = H5D_OPEN(file_identifier, "block size")
  dataspace = H5D_GET_SPACE(dataset)
  H5S_SELECT_HYPERSLAB, dataspace, start_2d, count_2d, STRIDE=stride_2d, /RESET
  memspace = H5S_CREATE_SIMPLE(count_2d)
  size = H5D_READ(dataset, FILE_SPACE=dataspace, MEMORY_SPACE=memspace)
  H5S_CLOSE, memspace
  H5S_CLOSE, dataspace
  H5D_CLOSE, dataset

  ; **********bounding box************
  dataset = H5D_OPEN(file_identifier, "bounding box")
  dataspace = H5D_GET_SPACE(dataset)
  H5S_SELECT_HYPERSLAB, dataspace, start_3d, count_3d, STRIDE=stride_3d, /RESET
  memspace = H5S_CREATE_SIMPLE(count_3d)  
  bnd_box = H5D_READ(dataset, FILE_SPACE=dataspace, MEMORY_SPACE=memspace)
  H5S_CLOSE, memspace
  H5S_CLOSE, dataspace
  H5D_CLOSE, dataset



endif else begin  ;end hdf5, start netcdf
  ncdf_varget, file_identifier, "coordinates" , coord
  ncdf_varget, file_identifier, "blocksize", size
  ncdf_varget, file_identifier, "bndbox", bnd_box
endelse

; save the coordinate information into the tree structure
for i = 0, ndim-1 do begin
  tree[*].coord[i] = reform(coord[i,*])
  tree[*].size[i] = reform(size[i,*])
  tree[*].bndBox[0,i] = reform(bnd_box[0,i,*])
  tree[*].bndBox[1,i] = reform(bnd_box[1,i,*])
endfor

undefine, coord
undefine, size
undefine, bnd_box

;---------------------------------------------------------------
; read in the unknowns - incorporates all hdf5/netcdf and
;                        single/double precision differences
;---------------------------------------------------------------
if (var_name EQ 'none') then  begin  ; read all variables
  ; Initialize 'unk' and use reform() with additional arguments to *force*
  ; IDL to create a 5-D array regardless of whether some of the dimensions
  ; (like 'totBlocks' or 'nzb' in a 2d-simulation), have size=1. IDL's 
  ; array-creation functions normally collapse these size-1 dimensions.
  if (not double) then begin
    unk = reform(fltarr(nvar,params.nxb,params.nyb,params.nzb,params.totBlocks), $
                 nvar, params.nxb, params.nyb, params.nzb, params.totBlocks)
  endif else begin
    unk = reform(dblarr(nvar,params.nxb,params.nyb,params.nzb,params.totBlocks), $
                 nvar, params.nxb, params.nyb, params.nzb, params.totBlocks)
  endelse 

  for i = 0, nvar-1 do begin
    if (file_type EQ 1) then begin  ; hdf5 
      dataset = H5D_OPEN(file_identifier, unk_names[i])
      var_data = H5D_READ(dataset)
      H5D_CLOSE, dataset
    endif else begin  ; netcdf
      ncdf_varget, file_identifier, unk_names[i], var_data
    endelse
    if Keyword_set(verbose) then help, unk_names[i], var_data
    if (not double) then begin
      unk[i,*,*,*,*] = temporary(float(var_data))
    endif else begin  ; double
      unk[i,*,*,*,*] = temporary(var_data)
    endelse
  endfor
endif else begin  ; read single variable
  if (file_type EQ 1) then begin  ; hdf5
    dataset = H5D_OPEN(file_identifier, unk_names[var])
    var_data = H5D_READ(dataset)
    H5D_CLOSE, dataset
  endif  else begin  ;netcdf
    ncdf_varget,file_identifier,unk_names[var], var_data
  endelse
  if Keyword_Set(verbose) then help, unk_names[var], var_data
  if (not double) then begin 
    unk = reform(float(var_data), 1, params.nxb, params.nyb, params.nzb, params.totBlocks)
  endif else begin
    unk = reform(var_data, 1, params.nxb, params.nyb, params.nzb, params.totBlocks)
  endelse
endelse 

; switch the block index to the front
unk = transpose(temporary(unk),[0,4,1,2,3])

; get the number of particles
numParticles = 0l
; sometimes we don't store the particle information -- check first
is_particles = 0  ; false
if (file_type EQ 1) then begin 
  ;; DEV update read_particles to handle both flash2 AND flash3 cases
  if (FLASH3) then begin
    read_particles, filename, NUM_PARTICLES=numParticles, PARTICLES=particles
    if Keyword_set(verbose) then begin 
      print,'in read_amr.pro, next help is on particles'
      help,particles,/STRUCTURE
    endif  
  endif else begin 
    file_contents = h5_parse(filename)
    if ((where(tag_names(file_contents) EQ "TRACER_PARTICLES"))[0] NE -1) then begin
      dataset = H5D_OPEN(file_identifier, "tracer particles")
      datatype = H5D_GET_TYPE(dataset)
      idl_type = H5T_IDLTYPE(datatype, STRUCTURE=particles)
      particles = H5D_READ(dataset)
      H5T_CLOSE, datatype
      H5D_CLOSE, dataset
      numParticles = (size(particles))[1]
    endif 
  endelse 
  H5F_CLOSE, file_identifier

  if (numParticles GT 0) then is_particles = 1 ; true

endif else begin  ; end HDF5, begin netcdf
  dimid = ncdf_dimid(file_identifier, 'dim_particles')
  if (dimid LT 0) then begin
    is_particles = 0  ; false
  endif else begin
    is_particles = 1  ; true
    ncdf_diminq, file_identifier, dimid, name, size
    numParticles = size
    if (numParticles GT 0) then begin
      dimid = ncdf_dimid(file_identifier, 'dim_npart_props')
      ncdf_diminq, file_identifier, dimid, name, numPartProps
      create_struct_string = 'particle = {'
      partPropNames = strarr(numPartProps)
      for i = 0, numPartProps-1 do begin
        attname = 'particle_props_name_' + strtrim(string(i), 2)
        ncdf_attget, file_identifier, /global, attname, name
        name = string(name)
        partPropNames[i] = name
        if (i eq numPartProps - 1) then begin
          create_struct_string = create_struct_string + name + ':double(0)'
        endif else begin
          create_struct_string = create_struct_string + name + ':double(0), '
        endelse  
      endfor            
      create_struct_string = create_struct_string + '}'        
      result = execute(create_struct_string)

      particles = replicate(particle, numParticles)
      partPropData = dblarr(numPartProps, numParticles)
      ncdf_varget,file_identifier,'particles', partPropData

      for i = 0L, numParticles-1L do begin
        for j = 0L, numPartProps-1L do begin
          assign_string = 'particles[i].' + partPropNames[j] + ' = partPropData[j,i]'
          result = execute(assign_string)
        endfor
      endfor 
    endif   
  endelse   
  NCDF_CLOSE, file_identifier
endelse  ; end of hdf5/netcdf split

if (not is_particles) then  begin ; no particles defined
  numParticles = 0l
  particles = 1
endif 

; ----------------------------------------------------------------
;  end of file reading
;------------------------------------------------------------------

;------------------------------------------------------------------------------
; compute the number of top level blocks in each direction
;------------------------------------------------------------------------------
top_blocks = where(tree[*].lrefine EQ 1)

ntopx = 1
ntopy = 1
ntopz = 1

; for 1d problems
case ndim of
  1: begin  ; 1d
    ntopx = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                        min(tree[*].bndBox[0,0]))))[1]
  end

  2: begin  ; 2d
    ; find the number of level 1 blocks whose lower y coord is the
    ; bottom of the domain
    ntopx1 = (size(where(tree[top_blocks].bndBox[0,1] EQ $
                         min(tree[*].bndBox[0,1]))))[1]

    ; now find the number of top level blocks at the upper y coord
    ntopx2 = (size(where(tree[top_blocks].bndBox[1,1] EQ $
                         max(tree[*].bndBox[1,1]))))[1]

    ; take the max of these, so we can deal with L shaped domains
    ntopx = ntopx1 > ntopx2

    ; find the number of level 1 blocks whose min x coord is the
    ntopy = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                        min(tree[*].bndBox[0,0]))))[1]
  end

  3: begin  ; 3d

    ; find the number of level 1 blocks whose lower y coord is the minimum
    ; y value and whose lower z coord is the minimum z value
    ntopx = (size(where(tree[top_blocks].bndBox[0,1] EQ $
                        min(tree[*].bndBox[0,1]) AND $
                        tree[top_blocks].bndBox[0,2] EQ $
                        min(tree[*].bndBox[0,2]))))[1]

    ; find the number of level 1 blocks whose lower x coord is the minimum
    ; x value and whose lower z coord is the minimum z value
    ntopy = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                        min(tree[*].bndBox[0,0]) AND $
                        tree[top_blocks].bndBox[0,2] EQ $
                        min(tree[*].bndBox[0,2]))))[1]

    ; find the number of level 1 blocks whose lower x coord is the minimum
    ; x value and whose lower y coord is the minimum y value
    ntopz = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                        min(tree[*].bndBox[0,0]) AND $
                        tree[top_blocks].bndBox[0,1] EQ $
                        min(tree[*].bndBox[0,1]))))[1]
  end
endcase

params.ntopx = ntopx
params.ntopy = ntopy
params.ntopz = ntopz

end  ; of read_amr.pro    
