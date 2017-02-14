;------------------------------------------------------------------------------
; file_information.pro
;
; probe a file and dump out some information about it -- we want to do this
; independently of any other routines, so do all the reading necessary 
; here. 
;
; If the silent keyword is set, no printing will be done, but the
; information will be available from the other keywords in the calling
; sequence
;------------------------------------------------------------------------------

pro file_information, filename, SILENT=silent, $
  CORNERS=corners, $
  DATE=date, $
  HDF_TYPE=fileFormat, $
  RUN_COMMENT=runComment, $
  FLASH_VERSION=flashVersion, $
  FILE_FORMAT_VERSION=fileFormatVersion, $
  DIMENSIONALITY=ndim, $
  NXB=nxb, NYB=nyb, NZB=nzb, $
  FILE_TYPE=fileType, $
  PRECISION=variablePrecision, $
  NUMBER_OF_VARIABLES=nvar, $
  NUMBER_OF_PARTICLES=numParticles, $
  BUILD_DATE=buildDate, $
  BUILD_DIR=buildDir, $
  BUILD_MACHINE=buildMachine, $
  SETUP_CALL=setupCall

; Define non-binary defaults so return and prints are always full
; these must
; be unique, otherwise IDL will give them all the same memory address,
; and each will contain the setupCall at the return -- silly IDL.
corners = "corners undefined"   ; "yes" or "no" -- badly determined
date = "date undefined"         ; -- not stored
; run comment isn't stored in flash2.5 or 3
runComment = "runComment undefined" ; -- not stored
flashVersion = "flashVersion undefined" ; -- not stored
fileFormat = "fileFormat undefined" ; HDF5 or NCMPI (NETCDF?)
fileFormatVersion = "fileFormatVersion undefined" ; version 7=2.5, version 8=3.0
fileType = "fileType undefined" ; "checkpoint" or "plotfile" -- badly determined
variablePrecision = "variablePrecision undefined" ; "single precision" or "double precision" -- appears incorrect (both double in flash3)
ndim = -1
nxb = -1 
nyb = -1 
nzb = -1
nvar = -1
numParticles = -1
buildDate = "buildDate undefined"
buildDir = "buildDir undefined"
buildMachine = "buildMachine undefined"
setupCall = "setupCall undefined"


xflash_dir = get_xflash_path()

; start by determining the file type -- whether HDF5 or NCMPI
type = determine_file_type(filename)

case type of
    1: fileFormat = "HDF5"
    2: fileFormat = "netCDF"	
    else: fileFormat = "something I don't know about"
ENDCASE

IF (type EQ -1) THEN BEGIN
    print, 'Sorry, file ',filename,' is of type ',fileFormat
    print, 'Bailing out'
    return      
ENDIF 
IF (type EQ 2) THEN BEGIN
    print, 'Sorry, I cannot read files of type netCDF yet'
    print, 'Bailing out (but other routines may still work)'
    return
ENDIF 

; determine precision -- this doesn't appear to be correct in flash 3
;  h5ls implies all are 64 bit variables
if (strpos(filename, 'plt_') GT 0) then begin
    fileType = 'plotfile (hack from filename)'
    variablePrecision = 'single precision (hack from filename)'
endif else begin
    fileType = 'checkpoint (hack from filename)'
    variablePrecision = 'double precision (hack from filename)'
ENDELSE


pathEnd = strpos(filename, '/', /reverse_search)

; in IDL, the default integer size is 2-bytes.  We need to use 4-byte
; integers when passing from IDL to C.
file_identifier = 0l

; open the file
if(type EQ 1) then begin
    file_identifier = H5F_OPEN(filename)

; determine the file format version FIRST as attributes vary with
; versions of FLASH
; version 7 in flash2.5, version 8 in (later) flash3

; set up a boolean for Flash Version
    FLASH2 = 0
    FLASH3 = 0

; we've got to look for the file format version first because
; in flash3, we've moved it into the sim info dataset

    fileFormatVersion = determine_file_version(filename)
    IF (fileFormatVersion GE 8) THEN FLASH3 = 1
    IF (fileFormatVersion LE 7) THEN FLASH2 = 1

    IF (FLASH2) THEN BEGIN 
        group = H5G_OPEN(file_identifier, "/")

        attribute = H5A_OPEN_NAME(group, "run comment")
        runComment = strtrim( H5A_READ(attribute), 2)
        H5A_CLOSE, attribute

        attribute = H5A_OPEN_NAME(group, "FLASH version")
        flashVersion = strtrim( H5A_READ(attribute), 2)
        H5A_CLOSE, attribute

        attribute = H5A_OPEN_NAME(group, "file creation time") ;
        date = strtrim( H5A_READ(attribute), 2)
        H5A_CLOSE, attribute

        attribute = H5A_OPEN_NAME(group, "FLASH build date")
        buildDate = strtrim( H5A_READ(attribute),2)
        H5A_CLOSE, attribute

        attribute = H5A_OPEN_NAME(group, "FLASH build directory")
        buildDir = strtrim( H5A_READ(attribute), 2)
        H5A_CLOSE, attribute

        attribute = H5A_OPEN_NAME(group, "FLASH build machine")
        buildMachine = strtrim( H5A_READ(attribute), 2)
        H5A_CLOSE, attribute

        attribute = H5A_OPEN_NAME(group, "FLASH setup call")
        setupCall = strtrim( H5A_READ(attribute), 2)
        H5A_CLOSE, attribute

        H5G_CLOSE, group
    ENDIF ; end of FLASh2
    IF (FLASH3) THEN BEGIN
        dataset = H5D_OPEN(file_identifier,"sim info")
        datatype = H5D_GET_TYPE(dataset)
        idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_info) ;; this line is not strictly necessary
        sim_info = H5D_READ(dataset)
        fileFormatVersion = strtrim(sim_info.file_format_version, 2)
        flashVersion = strtrim(sim_info.flash_version, 2)
        date = strtrim(sim_info.file_creation_time,2)
        buildDate = strtrim(sim_info.build_date,2)
        buildDir = strtrim(sim_info.build_dir,2)
        buildMachine = strtrim(sim_info.build_machine,2)
        setupCall = strtrim(sim_info.setup_call, 2)
        H5T_CLOSE, datatype
        H5D_CLOSE, dataset
    ENDIF ; end of FLASH3

;  corners.  There are lots of other things stored in logical
;  parameters, too
    dataset = H5D_OPEN(file_identifier, "logical runtime parameters")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=logical_parameters)
    logical_parameters = H5D_READ(dataset)
    H5D_CLOSE, dataset
    H5T_CLOSE, datatype

    for i=1, (size(logical_parameters))[3] do begin
        if (stregex(logical_parameters[i-1].name, '^corners', /BOOLEAN)) then corners_bool = logical_parameters[i-1].value
    endfor
    if (corners_bool) then begin
        corners = "yes"
    endif else begin
        corners = "no"
    endelse
; end of corners

; get the dimension
    ndim = 00
    dataset = H5D_OPEN(file_identifier, "coordinates")
    dataspace = H5D_GET_SPACE(dataset)
    dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)
    ndim = dims[0]

    H5D_CLOSE,dataset
    H5S_CLOSE, dataspace

    
; get the number of unknowns
    nvar = 0l
    dataset = H5D_OPEN(file_identifier, "unknown names")
    if (dataset LT 0) then fail = -1

    varnames =  H5D_READ(dataset)
    varnames = reform(varnames)
    nvar = (size(varnames))[1]
    H5D_CLOSE, dataset

; get the geometry
    geometry = determine_geometry(filename)


; read the header
    tot_blocks = 0l
    time       = 0.d0
    timestep   = 0.d0
    redshift   = 0.d0
    nsteps     = 0l
    nzones_per_block = lonarr(3)


    dataset = H5D_OPEN(file_identifier, "simulation parameters")
    datatype = H5D_GET_TYPE(dataset)
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_params)
    sim_params = H5D_READ(dataset)

    H5D_CLOSE, dataset
    H5T_CLOSE, datatype

    tot_blocks = sim_params.total_blocks
    time	   = sim_params.time
    timestep   = sim_params.timestep
    dt = timestep
    redshift   = sim_params.redshift
    nxb = sim_params.nxb
    nyb = sim_params.nyb
    nzb = sim_params.nzb

; get the number of particles
    numParticles = 0l

    IF (FLASH3) THEN BEGIN 
        dataset = H5D_OPEN(file_identifier, "integer scalars")
        datatype = H5D_GET_TYPE(dataset)
        idl_type = H5T_IDLTYPE(datatype, STRUCTURE=int_scalars)

        int_scalars = H5D_READ(dataset)
        H5D_CLOSE, dataset
        H5T_CLOSE, datatype

        for i=1, (size(int_scalars))[3] do begin
            if (stregex(int_scalars[i-1].name, '^nxb', /BOOLEAN)) then nxb = int_scalars[i-1].value
            if (stregex(int_scalars[i-1].name, '^nyb', /BOOLEAN)) then nyb = int_scalars[i-1].value
            if (stregex(int_scalars[i-1].name, '^nzb', /BOOLEAN)) then nzb = int_scalars[i-1].value
        endfor

        dataset = H5D_OPEN(file_identifier, "real scalars")
        datatype = H5D_GET_TYPE(dataset)
        idl_type = H5T_IDLTYPE(datatype, STRUCTURE=real_scalars)

        real_scalars = H5D_READ(dataset)
        H5D_CLOSE, dataset
        H5T_CLOSE, datatype

        for i=1, (size(real_scalars))[3] do begin
            if (stregex(real_scalars[i-1].name, '^time', /BOOLEAN)) then time = real_scalars[i-1].value
            if (stregex(real_scalars[i-1].name, '^dt', /BOOLEAN)) then dt = real_scalars[i-1].value
        endfor

;        time = real_scalars.time
;        dt = real_scalars.dt
    ENDIF ; end of FLASH3

    H5F_CLOSE, file_identifier
    
ENDIF ; end of HDF5 -- no netCDF implemented

; See Coyote webpage for good explanation of binary keywords
; http://www.dfanning.com/tips/keyword_check.html
if (NOT Keyword_set(silent)) then begin

    print, ' '
    print, '-------------------------------------------------------------------------------'
    print, ' file name  = ',  strmid(filename,pathEnd+1,strlen(filename)-pathEnd)
    print, ' '
    print, ' FLASH version:       ', flashVersion
    print, ' file format:         ', fileFormat
    print, ' file format version: ', fileFormatVersion
    print, ' '
    print, ' build date:          ', buildDate
    print, ' build directory:     ', buildDir
    print, ' build machine:       ', buildMachine
    print, ' '
    print, ' setup call:          ', setupCall
    print, ' '
    print, ' execution date:      ', date
    print, ' run comment:         ', runComment
    print, ' '
    print, ' dimension:           ', ndim
    print, ' geometry detected:   ', geometry
    print, ' '
    print, ' type of file:        ', fileType
    print, ' '
    print, ' number of variables: ', nvar
    print, ' variable precision:  ', variablePrecision
    print, ' '
    print, ' number of particles: ', numParticles, ' not implemented!'
    print, ' '
    print, ' nxb, nyb, nzb:       ', nxb, nyb, nzb
    print, ' corners stored:      ', corners
    print, ' '
    print, ' time:                ', time
    print, ' dt:                  ', dt
    print, '-------------------------------------------------------------------------------'
endif

end



