g;=================================================
; Copyright Mark Krumholz (2001) 
; Modified for FLASH use, Robert Fisher (2005)
;=================================================

;==========================================================================
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
;============================================================================

; Routine to read a single amr component from a file 
; and store it in an IDL structure. This structure can then be passed 
; to other amrlib routines to do various useful things. The argument
; filename is the name of the file to be read, and componentName
; is the name of the component to be read.

; Modified to handle FLASH files, based upon the FLASH-supplied FIDLR program
; read_amr, and using the FLASH-supplied FIDLR function determine_file_type.
; HDF5 files use built-in IDL HDF5 routines, which requires v 5.6 of IDL
; or later.

;-----------------------------------
; IMPORTANT ASSUMPTIONS/RESTRICTIONS
;-----------------------------------

; For FLASH, we assume
;   (*) All cell-centered quantities. The "corner" option in FLASH files
;         is set to a hard zero.
;   (*) Further, for the purposes of internal representation, the levels are
;         defined from 0 to maxlevel - 1, as opposed to 1 to maxlevel.
;   (*) The number of timesteps is not stored; evidently it is not in 
;         older FLASH files (<= 2.4), although it is in the 2.5 manual.
;   (*) Lastly, we assume all datasets are Cartesian. No support exists
;         for spherical or cylindrical geometry.

; Only HDF 5 currently supported; CDF lines not yet adopted from read_amr.

; For reference, these limitations are denoted throughout the code by
;  comments marked "LIMITATIONS".

;-------------------
; Internal storage
;-------------------

; The highest level object is an amr descriptor. It consistents of a
; series of levels, and also stores global information about the AMR grid.

; The next level object is a level. It consists of a series of fabs (
; or "blocks" in PARAMESH-parlance.) It also stores various information
; about that level.
                                                                                
; The lowest level structure is a fab (or "block" in PARAMESH-parlance).
; It holds a single block with a specified index range, and physical range,
; and data.
                                                                                

;==========================================================================
;
; Begin code. First, inspect and read in file.
;
;==========================================================================

if not keyword_set(verbose) then begin
  verbose = 1
endif

if (verbose) then begin
  print, 'Verbose output.'
endif

;-------------------------------------
; check to see which variable to read
;-------------------------------------

if n_elements(componentName) EQ 0 then componentName = 'none'

if (componentName eq 'none') then begin
  print, 'Error : no component name to read in specified.'
  return, -1
endif

;---------------------------------------------------------
; Get the file type.
;  itype = 1 HDF5; itype = 2, NCDF. Otherwise -1.   
;---------------------------------------------------------
                                      
itype = determine_file_type(filename)
                                                                                
if(itype EQ -1) then begin
  print,'Error: Problem determining filetype of file ',filename
  exit 
endif

if (verbose) then begin
  if (itype EQ 1) then print, 'Recognized HDF5 file format.'
  if (itype EQ 2) then print, 'Recognized NCDF file format.'
endif

;------------------------------------------------------------
; read in FLASH data using code from FIDLR program read_amr.
;------------------------------------------------------------

if(itype EQ 1) then begin
                                                                                
;------------------------------------------------------------------
; open up the file for read now and read in the header information
;------------------------------------------------------------------                                                                             
  file_identifier = H5F_OPEN(filename)

;--------------------------------------------------------------
; Read in FLASH annotations.
;--------------------------------------------------------------
   
flashVersion = '                    '
date = '                                        '
runComment = 'run comment is not stored               '
buildDate = '1 build date not defined              '
buildDir =  '2 build directory not defined              '
buildMachine = '3 build machine not defined                 '
setupCall = '4 setup call not defined              '
                                                                                
group = H5G_OPEN(file_identifier, "/")
attribute = H5A_OPEN_NAME(group, "run comment")
runComment = H5A_READ(attribute)
H5A_CLOSE, attribute
                                                                                
attribute = H5A_OPEN_NAME(group, "FLASH version")
flashversion = H5A_READ(attribute)
H5A_CLOSE, attribute

attribute = H5A_OPEN_NAME(group, "file creation time");
date = H5A_READ(attribute)
                                                                                
H5A_CLOSE, attribute
H5G_CLOSE, group
                                                                                
group = H5G_OPEN(file_identifier, "/")
attribute = H5A_OPEN_NAME(group, "FLASH build date")
buildDate = H5A_READ(attribute)
H5A_CLOSE, attribute
                                                                                
attribute = H5A_OPEN_NAME(group, "FLASH build directory")
buildDir = H5A_READ(attribute)
H5A_CLOSE, attribute
                                                                                
attribute = H5A_OPEN_NAME(group, "FLASH build machine")
buildMachine = H5A_READ(attribute)
H5A_CLOSE, attribute
                                                                                
attribute = H5A_OPEN_NAME(group, "FLASH setup call")
setupCall = H5A_READ(attribute)
H5A_CLOSE, attribute
                                                                                
H5G_CLOSE, group

if (verbose) then begin
  print, 'Flash version = ', flashVersion
  print, 'Run comment = ', runComment
  print, 'File creation date = ', date
  print, 'FLASH build date = ', buildDate
  print, 'FLASH build directory = ', buildDir
  print, 'FLASH build machine = ', buildMachine
  print, 'FLASH setup call = ', setupCall
endif

;--------------------------------------------------------------
; grab the list of variables, and make sure that our requested
; variable (if any) exists in the dataset
;--------------------------------------------------------------
                                                                                
  dataset = H5D_OPEN(file_identifier, "unknown names")
                                                                                
  unk_names =  H5D_READ(dataset)
  unk_names = reform(temporary(unk_names))
                                                                                
  H5D_CLOSE, dataset

;---------------------------------------
; get the dimensionality of the dataset
;---------------------------------------
                                                                                
  dataset = H5D_OPEN(file_identifier, "coordinates")
  dataspace = H5D_GET_SPACE(dataset)
  dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)
                                                                                
  ndim = dims[0]
  
  nfaces = 2*ndim
  nchild = 2^ndim
                                                                              
  H5S_CLOSE, dataspace
  H5D_CLOSE, dataset

  if (ndim ne 1) and (ndim ne 2) and $
  (ndim ne 3) then begin
    print, 'Error: only 1, 2, or 3 dimensional plots are supported.'
    return, -1
  endif

  if (verbose) then begin
    print, '# Dimensions = ', ndim
  endif                                                                                
;------------------------------------------------------------------
; we can also now set the number of variables and determine index
;  for component to be read in
;------------------------------------------------------------------

  nvar = (size(unk_names))[1]

  if n_elements(componentName) EQ 0 then componentName = 'none'
                                                                                
  if (componentName NE 'none') then begin
    index = (where(unk_names EQ componentName))[0]
                                                                                
    if (index EQ -1) then begin
      print, 'ERROR: requested variable not found in dataset'
      return, -1
    endif

    if (verbose) then begin
      print, 'stored components : ', unk_names                    
      print, 'reading in only the variable ', componentName
      print, 'component # ', index
    endif

  endif
                                                                                
;----------------------------------                                            
; read in the simulation paramters
;----------------------------------

  dataset = H5D_OPEN(file_identifier, "simulation parameters")
  datatype = H5D_GET_TYPE(dataset)
                                                                                
  idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_params)
                                                                                
  sim_params = H5D_READ(dataset)
                                                                                
  H5T_CLOSE, datatype
  H5D_CLOSE, dataset

;------------------
; read in geometry
;------------------

  if (n_elements(geometry) eq 0) then BEGIN
    geometry = determine_geometry(filename)
;    group = H5G_OPEN(file_identifier, "/")
;    attribute = H5A_OPEN_NAME(group, "geometry name")
;                                                                                
;    geometry = H5A_READ(attribute)
;                                                                                
;    H5A_CLOSE, attribute
;    H5G_CLOSE, group
  endif

;----------------------------------------------------------------------------
; Following FIDLR read_amr, build a "tree" structure store AMR data internally.
;
; Note that this term is ill-chosen, as the structure is not actually
; a tree, but simply a record consisting of various arrays.
;-----------------------------------------------------------------------------

  if (n_elements(double) EQ 0) then double = 0
                                                                                
                                                                                
  if (NOT double) then begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            processorNumber:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:fltarr(ndim), $
            size:fltarr(ndim), $
            bndBox:fltarr(2,ndim)}
  endif else begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            processorNumber:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:dblarr(ndim), $
            size:dblarr(ndim), $
            bndBox:dblarr(2,ndim)}
  endelse
                                                                                
  tree = replicate(tree,sim_params.total_blocks)

; LIMITATION
; set corners = 0 by default

  corners = 0

  if (NOT double) then begin
    params = {totBlocks:sim_params.total_blocks, $
              corners:corners, $
              ndim:ndim, $
              nvar:nvar, $
              nxb:sim_params.nxb, $
              nyb:sim_params.nyb, $
              nzb:sim_params.nzb, $
              ntopx:1, $
              ntopy:1, $
              ntopz:1, $
              time:float(sim_params.time), $
              dt:float(sim_params.timestep), $
              redshift:1.0, $
              geometry:"unknown"}
  endif else begin
      params = {totBlocks:sim_params.total_blocks, $
              corners:corners, $
              ndim:ndim, $
              nvar:nvar, $
              nxb:sim_params.nxb, $
              nyb:sim_params.nyb, $
              nzb:sim_params.nzb, $
              ntopx:1, $
              ntopy:1, $
              ntopz:1, $
              time:sim_params.time, $
              dt:sim_params.timestep, $
              redshift:1.d0, $
              geometry:"unknown"}
  endelse
                                                                                
  params.geometry = strcompress(geometry, /REMOVE_ALL)
                                                                                
; for some reason, the redshift field is not stored in the plotfiles
  if ( (where(tag_names(sim_params) EQ "REDSHIFT"))[0] NE -1) then begin
                                                                                
      if (NOT double) then begin
        params.redshift = float(sim_params.redshift)
      endif else begin
        params.redshift = sim_params.redshift
      endelse
                                                                                
  endif
                                                                                
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; We now have the dimension, total number of blocks, and the number of
; zones per block in each direction.  We can use this to read the tree
; information, coordinate information, and unknowns
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                                                                                
;----------------------------------------------------------------------------
; read in the tree information
;----------------------------------------------------------------------------

  ngid = long(2^ndim + 1 + 2*ndim)
                                                                                
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
                                                                                
  dataset = H5D_OPEN(file_identifier, "coordinates")
  coord = H5D_READ(dataset)
  H5D_CLOSE, dataset
                                                                                
  dataset = H5D_OPEN(file_identifier, "block size")
  size = H5D_READ(dataset)
  H5D_CLOSE, dataset
                                                                                
  dataset = H5D_OPEN(file_identifier, "bounding box")
  bnd_box = H5D_READ(dataset)
  H5D_CLOSE, dataset

;----------------------------------------------------------
; save the coordinate information into the tree structure
;----------------------------------------------------------

  for i = 0, ndim-1 do begin
    tree[*].coord[i] = reform(coord[i,*])
    tree[*].size[i] = reform(size[i,*])
    tree[*].bndBox[0,i] = reform(bnd_box[0,i,*])
    tree[*].bndBox[1,i] = reform(bnd_box[1,i,*])
  endfor

  undefine, coord
  undefine, size
  undefine, bnd_box

;-------------------------------------------------------------------------
; LIMITATION
; Assume Cartesian geometry, and determine the problem domain box extent
;  by inspecting the bounding box for all blocks.
;-------------------------------------------------------------------------

  boxmin = fltarr (ndim)
  boxmax = fltarr (ndim)

  for i = 0, ndim - 1 do begin
    boxmin [i] = min (tree [*].bndBox [0, i] )
    boxmax [i] = max (tree [*].bndBox [1, i] )
  endfor

;-------------------------
; read in the unknowns
;-------------------------

  if (verbose) then begin
    print, 'Total number of blocks = ', params.totBlocks
  endif                                                                                
  if (NOT double) then begin
        unk = fltarr(params.nxb,params.nyb,params.nzb,params.totBlocks)
                                                                                
        dataset = H5D_OPEN(file_identifier, unk_names[index])
        var = H5D_READ(dataset)
        H5D_CLOSE, dataset

        if (verbose) then begin                                            
          help, var
        endif
                                                                        
        unk = reform(temporary(float(var)), $
                     1,params.nxb,params.nyb,params.nzb,params.totBlocks)
                                                                                
  endif else begin
        unk = dblarr(params.nxb,params.nyb,params.nzb,params.totBlocks)
                                                                                
        dataset = H5D_OPEN(file_identifier, unk_names[var])
        var = H5D_READ(dataset)
        H5D_CLOSE, dataset
                                                                                
        help, var
                                                                                
        unk = reform(temporary(unk), $
                     1,params.nxb,params.nyb,params.nzb,params.totBlocks)
  endelse

; sometimes, a simulation exists that has only a single block.  IDL
; has this annoying habit or automagically dropping a dimension from
; an array if it's index is 1.  A lot of other routines want this
; index to be here, so restore it.

  if (params.totBlocks EQ 1) then begin
        unk = reform(unk, 1,params.nxb,params.nyb,params.nzb,1)
  endif
                                                                                
; switch the block index to the front
  unk = transpose(temporary(unk),[0,4,1,2,3])
                                                                                

; LIMITATION - We read in the particles, but do not store them yet.        

; get the number of particles
  numParticles = 0l

; sometimes we don't store the particle information -- check first
  file_contents = h5_parse(filename)

  if ( (where(tag_names(file_contents) EQ $
            "TRACER_PARTICLES"))[0] NE -1) then begin
                                                                                
    dataset = H5D_OPEN(file_identifier, "tracer particles")
    datatype = H5D_GET_TYPE(dataset)
                                                                                
    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=particles)
                                                                                
    particles = H5D_READ(dataset)
                                                                                
    H5T_CLOSE, datatype
    H5D_CLOSE, dataset
                                                                                
    numParticles = (size(particles))[1]
  endif else begin
    numParticles = 0l
    particles = 1
  endelse
                                                                                
                                                                                
  H5F_CLOSE, file_identifier

endif  ; end if itype = 1

if (verbose) then begin
  print, 'File read complete.'
endif

;--------------------------------------------------------------------------
; 
; File read complete. Now, translate FLASH data into internal IDLAMRLIB
;  data representation.
;
;--------------------------------------------------------------------------

; Convert the levels from [1, LMAX] to [0, (LMAX - 1)] notation. This is
;   simply a difference between PARAMESH and Berger/Colella conventions.

maxlevel = max (tree [*].lrefine) - 1  

;------------------------------------------------------------------------------
; Compute the number of top level blocks in each direction. Assuming
;  a logically Cartesian domain, this gives us the base problem resolution.
;  Also, for purposes of the gridDescriptor array defined below, we need
;  to know the maximum number of grids on any level - store this in maxgrids.
;------------------------------------------------------------------------------

maxgrids = 0

for n = 1, maxlevel + 1 do begin
   blocks_on_this_level = (size (where (tree [*].lrefine EQ n) ) ) [1]
   maxgrids = max ([blocks_on_this_level, maxgrids])
endfor

if (verbose) then begin
  print, 'maxlevel (base level 0) = ', maxlevel
  print, 'Maximum number of blocks per level = ', maxgrids
endif

top_blocks = where(tree[*].lrefine EQ 1)
                                                                                
ntopx = 1
ntopy = 1
ntopz = 1
                                                                                
; for 1d problems
case ndim of
    1: begin
        ntopx = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                            min(tree[*].bndBox[0,0]))))[1]
        ntop  = [ntopx]
    end
                                                                                
;for 2d problems
    2: begin
                                                                                
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
        ntop = [ntopx, ntopy]
    end
                                                                                
                                                                                
;for 3d problems
    3: begin
                                                                                
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
        ntop = [ntopx, ntopy, ntopz]
                                                                                
    end
endcase

;------------------------------------------------------------------- 
; Create sructures in which to hold the amr data
;-------------------------------------------------------------------

; establish a structure to hold sink particles
sinkParticleDescriptor = { x:dblarr(ndim), p:dblarr(ndim), $
	j:dblarr(ndim), m:double(0.0) }

; the lowest level structure is a fab. It holds a single grid with a
; specified index range, and physical range, and data.

fabDescriptor = { idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	idxtype:lonarr(ndim), xlo:dblarr(ndim), xhi:dblarr(ndim), $
	dataptr:ptr_new(/allocate_heap) }
ptr_free, fabDescriptor.dataptr

; note: dataptr:ptr_new(/allocate_heap) is roughly equivalent to void *dataptr
; However, IDL insists on allocating valid memory to the pointer, which we must
; free up since we don't want it and will do our own memory assignment later.

; the next level object is a level. It consists of a series of fabs. It
; also stores various information about that level.

; create a generic grid descriptor
gridDescriptor = { xlo:dblarr(ndim), xhi:dblarr(ndim) }

levelDescriptor = { level:0, ngrids:0, leveltime:0.0d, levelsteps:0L,  $
	idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	periodicity:lonarr(ndim), gridspacing:dblarr(ndim), $
	grids:replicate(gridDescriptor,maxgrids), $
	nfab:0, fabptr:ptr_new(/allocate_heap) }
ptr_free, levelDescriptor.fabptr	; free up dummy memory

; The highest level object is an amr descriptor. It consists of a
; series of levels, and also stores global information about the amr grid.

; We have to have four cases here because IDL doesn't have a goddamn void *
; pointer! 

; RTF -- this has been reduced to 2 for FLASH, since there are no star
; particles.

nsink = 0 ; LIMITATION skip particles for now in FLASH
nstar = 0 ; no star particles in FLASH

idxlo=lonarr(ndim,maxlevel+1)
idxhi=lonarr(ndim,maxlevel+1)

;----------------------------------------------------
; Set refinement ratio for datasets with refinement
;----------------------------------------------------

if maxlevel gt 0 then refratio=lonarr(maxlevel+1) else refratio=0

; Because PARAMESH imposes equal block sizes over all levels, we can
; compute the refinement ratio from the spatial extent of any two boxes
; on two different levels in any direction.

; Compute refinement ratio across all levels.

for n = 1, maxlevel  do begin

  refratio [n] =  $
 ( $
  (tree ( (where(tree[*].lrefine EQ n) )      [0]).bndBox [1, 0] - $
   tree ( (where (tree[*].lrefine EQ n) )     [0]).bndBox [0, 0] ) / $
  (tree ( (where(tree[*].lrefine EQ n + 1) )  [0]).bndBox [1, 0] - $
   tree ( (where (tree[*].lrefine EQ n + 1) ) [0]).bndBox [0, 0] ) $
 )

endfor

product = 1.

if (verbose) then begin

  for n = 1, maxlevel + 1 do begin
    blocks_on_this_level = (size (where (tree [*].lrefine EQ n) ) ) [1]
    if n gt 1 then product = product / (refratio [n - 1])^ndim
    if n gt 1 then filling_factor = filling_factor * product * blocks_on_this_level $
     else filling_factor = 1.
    if (verbose) then begin
      print, 'blocks on level ', n, ' = ', blocks_on_this_level, filling_factor * 100., ' % of domain'
    endif

  endfor
endif
                                                                                                           
nb = [sim_params.nxb, sim_params.nyb, sim_params.nzb]

for l=0, ndim-1 do begin
  idxlo [l, 0] = 0 
  idxhi [l, 0] = ntop [l] * nb [l] - 1
endfor

for n = 1, maxlevel do begin

  for l=0, ndim-1 do begin
    idxlo [l, n] = 0
    idxhi [l, n] = (refratio [n] * (idxhi [l, n - 1] + 1) ) - 1
  endfor

endfor

nsteps = 0

dx=fltarr(ndim,maxlevel+1) 

dx (*, 0) = (boxmax (*) - boxmin (*) ) / (idxhi (*, 0) + 1.)

for n = 1, maxlevel do begin
  dx (*, n) = dx (*, n - 1) / refratio [n]
endfor

if (verbose) then begin
  print, 'grid spacing over all directions'
  print,  dx
endif
 
if (nsink eq 0) then begin
	if (nstar eq 0) then begin
		; no stars or sinks
		amrDescriptor = { name:filename, $
		componentName:componentName, $
		version:flashversion, $
		ndim:ndim, $
		time:sim_params.time, $
		maxlevel:maxlevel, $
		boxmin:boxmin, boxmax:boxmax, $
		refratio:refratio, $
		idxlo:idxlo, idxhi:idxhi, $
		idxtype:intarr (ndim), $
		periodicity:0, $
		levelsteps:nsteps, $
		gridspacing:dx, $
		coordtype:0, $
		levels:replicate(levelDescriptor,maxlevel+1), $
		nsink:0, sinkparticles:ptr_new(), $
		nstar:0, starparticles:ptr_new() }
	endif
endif else begin
	if (nstar eq 0) then begin
		; sinks but no stars
		amrDescriptor = { name:filename, $
		componentName:componentName, $
		version:flashversion, $
		ndim:ndim, $
		time:sim_params.time, $
		maxlevel:maxlevel, $
		boxmin:boxmin, boxmax:boxmax, $
		refratio:refratio, $
		idxlo:idxlo, idxhi:idxhi, $
		idxtype:intarr(ndim), $
		periodicity:0, $
		levelsteps:nsteps, $
		gridspacing:dx, $
		coordtype:0, $
		levels:replicate(levelDescriptor,maxlevel+1), $
		nsink:nsink, $
		sinkparticles:replicate(sinkParticleDescriptor,nsink), $
		nstar:0, starparticles:ptr_new() }
        endif
endelse

;----------------------------------------------------------------------
; The next thing we need is the byte ordering. From what I can tell,
; the object in the second set of parentheses is always EITHER
; (nbytes, (1 2 3 ... nbytes)) OR
; (nbytes, (nbytes nbytes-1 nbytes-2 ... 1)).
; The first of these means little endian (same as the Sun and Intel),
; the second means big endian (same as AIX).
;-----------------------------------------------------------------------                                                                                
temp=''

temp=strmid(temp, strpos(temp, "))")+2)
temp=strmid(temp, strpos(temp, "(")+1)
temp1=strmid(temp, 0, strpos(temp, ","))
nbytes2=long(temp1)
temp=strmid(temp, strpos(temp, "(")+1)
temp1=strmid(temp, 0, strpos(temp, " "))
firstnum=long(temp1)
if firstnum ne nbytes2 then bigendian=1 else bigendian=0

;------------------------------------------------------------------------    
; This uses the non-standard routine is_ieee_big to determine the
; endian-ness of the local machine. A copy of the routine is included
; with amrlib, but I didn't write it.
;------------------------------------------------------------------------
                                                                                
if bigendian ne is_ieee_big() then do_swap_endian = 1 else do_swap_endian = 0

;---------------------------------------------------------
;  convert data from FLASH "tree" to amrDescriptor format
;---------------------------------------------------------

for n=0, maxlevel do begin

	; Fill out level descriptors

	amrDescriptor.levels[n].level=n
	amrDescriptor.levels[n].ngrids = $
                (size (where(tree[*].lrefine EQ n + 1) )) [1]
	amrDescriptor.levels[n].leveltime = $
		sim_params.time
	amrDescriptor.levels[n].levelsteps = $
	        nsteps	
	amrDescriptor.levels[n].idxlo = $
		idxlo (*, n)
	amrDescriptor.levels[n].idxhi = $
		idxhi (*, n)
	amrDescriptor.levels[n].periodicity = $
		0
	amrDescriptor.levels[n].gridspacing = $
	        dx (*, n)	
	amrDescriptor.levels[n].nfab = $
                (size (where (tree [*].lrefine EQ n + 1) )) [1]

	; we now know the number of fabs, so create an array of
	; fab descriptors for the fabptr to point to.

	amrDescriptor.levels[n].fabptr = $
		ptr_new(replicate(fabDescriptor,amrDescriptor.levels[n].nfab))

	; loop through the fabs and fill them

	for k=0, amrDescriptor.levels[n].nfab-1 do begin

                fabno = (where (tree [*].lrefine EQ n + 1)) [k] 

		; record properties of each fab. idxlo and idxhi
		; are the index limits, idxtype is the index type
		; (0 = cell centered, 1 = edge centered), and
		; xlo and xhi are the physical limits

                ; LIMITATION
                ; For FLASH, we assume cell-centered quantities for now

		(*amrDescriptor.levels[n].fabptr)[k].idxtype = $
			0 
		(*amrDescriptor.levels[n].fabptr)[k].xlo = $
                   tree [fabno].bndBox [0, *]
		(*amrDescriptor.levels[n].fabptr)[k].xhi = $
                   tree [fabno].bndBox [1, *]

                amrDescriptor.levels[n].grids[k].xlo = $
                   tree [fabno].bndBox [0, *]
                amrDescriptor.levels[n].grids[k].xhi = $
                   tree [fabno].bndBox [1, *]

                (*amrDescriptor.levels[n].fabptr)[k].idxlo = $
                  floor( (*amrDescriptor.levels[n].fabptr)[k].xlo / dx (*, n) )
                (*amrDescriptor.levels[n].fabptr)[k].idxhi = $
                  ceil ( (*amrDescriptor.levels[n].fabptr)[k].xhi / dx (*, n) ) - 1

                localfab = fltarr(params.nxb,params.nyb,params.nzb)

                localfab (*, *, *) = unk (0, fabno, *, *, *)  

		; read in fab data, simultaneously allocating memory to hold it
		(*amrDescriptor.levels[n].fabptr)[k].dataptr = $
                   ptr_new (localfab)

	endfor

endfor

; Set the index type at the top level. In the plot file format, this
; is stored fab by fab, since different components can have different
; types (cell vs. edge centered) within the same plot file. However,
; since amrDescriptor applies to only a single component, it makes more
; sense for the index type to be stored at the top level -- it will be
; the same for every fab in any case.

amrDescriptor.idxtype = (*amrDescriptor.levels[0].fabptr)[0].idxtype

; return data
return, amrDescriptor
end


; Copyright Mark Krumholz (2001)
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

; routine to read a single amr component from a plot file / directory
; and store it in an IDL structure. This structure can then be passed 
; to other amrlib routines to do various useful things. The argument
; dirname is the name of the plot directory to be read, and componetName
; is the name of the component to be read.

; read plotfile header
plotdescriptor = read_amrheader(dirname)
if (plotdescriptor.ndim ne 1) and (plotdescriptor.ndim ne 2) and $
	(plotdescriptor.ndim ne 3) then begin
	print, 'Error: only 1, 2, or 3 dimensional plots are supported.'
	return, -1
endif
ndim = plotdescriptor.ndim

; get component number and check that it exists
component = where(plotDescriptor.quantities eq componentName)
component = component[0]
if component eq -1 then begin
	print, 'Component '+componentName+' does not exist.'
	print, 'Valid components are:'
	print, plotDescriptor.quantities
	return, -1
endif

; figure out which family of multifab files contains this component,
; and which component number within that family it is
fileNameIndex = 0
while component ge $
	plotDescriptor.levels[0].nfilecomp[fileNameIndex] do begin
	component = component - $
		plotDescriptor.levels[0].nfilecomp[fileNameIndex]
	fileNameIndex = fileNameIndex + 1
endwhile

; read number of sink particles if sink file exists
nsink=0
if strpos(dirname, '/', /reverse_search) eq strlen(dirname)-1 then $
	sinkname=dirname+'SinkParticles' $
else sinkname=dirname+'/SinkParticles'
if (file_test(sinkname)) then begin
	openr, fp, sinkname, /get_lun
	readf, fp, nsink
	free_lun, fp
endif

; read number of star particles if star file exists
nstar=0
if strpos(dirname, '/', /reverse_search) eq strlen(dirname)-1 then $
	starname=dirname+'StarParticles' $
else starname=dirname+'/StarParticles'
if (file_test(starname)) then begin
	openr, fp, starname, /get_lun
	readf, fp, nstar
	free_lun, fp
endif

; create a new structure in which to hold the amr data

; establish a structure to hold sink particles
sinkParticleDescriptor = { x:dblarr(ndim), p:dblarr(ndim), $
	j:dblarr(ndim), m:double(0.0) }

; establish a structure to hold star particles
starParticleDescriptor = { m:double(0.0), x:dblarr(ndim), p:dblarr(ndim), $
	j:dblarr(ndim), mlast:double(0.0), r:double(0.0), mdeut:double(0.0), $
	n:double(0.0), mdot:double(0.0), burnState:0 }

; the lowest level structure is a fab. It holds a single grid with a
; specified index range, and physical range, and data.

fabDescriptor = { idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	idxtype:lonarr(ndim), xlo:dblarr(ndim), xhi:dblarr(ndim), $
	dataptr:ptr_new(/allocate_heap) }
ptr_free, fabDescriptor.dataptr

; note: dataptr:ptr_new(/allocate_heap) is roughly equivalent to void *dataptr
; However, IDL insists on allocating valid memory to the pointer, which we must
; free up since we don't want it and will do our own memory assignment later.

; the next level object is a level. It consists of a series of fabs. It
; also stores various information about that level.

levelDescriptor = { level:0, ngrids:0, leveltime:0.0d, levelsteps:0L,  $
	idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	periodicity:lonarr(ndim), gridspacing:dblarr(ndim), $
	grids:plotdescriptor.levels[0].grids, $
	nfab:0, fabptr:ptr_new(/allocate_heap) }
ptr_free, levelDescriptor.fabptr	; free up dummy memory

; The highest level object is an amr descriptor. It consistents of a
; series of levels, and also stores global information about the amr grid.

; We have to have four cases here because IDL doesn't have a goddamn void *
; pointer!

if (nsink eq 0) then begin
	if (nstar eq 0) then begin
		; no stars or sinks
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:0, sinkparticles:ptr_new(), $
		nstar:0, starparticles:ptr_new() }
	endif else begin
		; stars but no sinks
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:0, sinkparticles:ptr_new(), $
		nstar:nstar, $
		starparticles:replicate(starParticleDescriptor,nstar) }
	endelse
endif else begin
	if (nstar eq 0) then begin
		; sinks but no stars
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:nsink, $
		sinkparticles:replicate(sinkParticleDescriptor,nsink), $
		nstar:0, starparticles:ptr_new() }
	endif else begin
		; both sinks and stars -- this probably shouldn't ever
		; happen, but we include the possibility anyway
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:nsink, $
		sinkparticles:replicate(sinkParticleDescriptor,nstar), $
		nstar:nstar, $
		starparticles:replicate(starParticleDescriptor,nstar) }
	endelse
endelse

; now read all the data

for n=0, amrDescriptor.maxlevel do begin

	; read header file
	mfhdr = read_multifabHeader(n, fileNameIndex, plotdescriptor)

	; store useful info
	amrDescriptor.levels[n].level=n
	amrDescriptor.levels[n].ngrids = $
		plotdescriptor.levels[n].ngrids
	amrDescriptor.levels[n].leveltime = $
		plotdescriptor.levels[n].leveltime
	amrDescriptor.levels[n].levelsteps = $
		plotdescriptor.levels[n].levelsteps
	amrDescriptor.levels[n].idxlo = $
		plotdescriptor.levels[n].idxlo
	amrDescriptor.levels[n].idxhi = $
		plotdescriptor.levels[n].idxhi
	amrDescriptor.levels[n].periodicity = $
		plotdescriptor.levels[n].periodicity
	amrDescriptor.levels[n].gridspacing = $
		plotdescriptor.levels[n].gridspacing
	amrDescriptor.levels[n].grids = $
		plotdescriptor.levels[n].grids
	amrDescriptor.levels[n].nfab = mfhdr.nfab

	; we now know the number of fabs, so create an array of
	; fab descriptors for the fabptr to point to.
	amrDescriptor.levels[n].fabptr = $
		ptr_new(replicate(fabDescriptor,mfhdr.nfab))

	; loop through the fabs
	for l=0, amrDescriptor.levels[n].nfab-1 do begin

		; record properties of each fab. idxlo and idxhi
		; are the index limits, idxtype is the index type
		; (0 = cell centered, 1 = edge centered), and
		; xlo and xhi are the physical limits
		(*amrDescriptor.levels[n].fabptr)[l].idxlo = $
			mfhdr.fabs[l].idxlo
		(*amrDescriptor.levels[n].fabptr)[l].idxhi = $
			mfhdr.fabs[l].idxhi
		(*amrDescriptor.levels[n].fabptr)[l].idxtype = $
			mfhdr.fabs[l].idxtype
		(*amrDescriptor.levels[n].fabptr)[l].xlo = $
			plotdescriptor.boxmin + $
			plotdescriptor.gridspacing[*,n] * mfhdr.fabs[l].idxlo
		(*amrDescriptor.levels[n].fabptr)[l].xhi = $
			plotdescriptor.boxmin + $
			plotdescriptor.gridspacing[*,n] * $
			(mfhdr.fabs[l].idxhi + $
				(mfhdr.fabs[l].idxtype eq 0))

		; read in fab data, simultaneously allocating memory to hold it
		(*amrDescriptor.levels[n].fabptr)[l].dataptr = ptr_new( $
			read_fabcomponent(l, component, plotdescriptor, mfhdr))

	endfor

endfor

; Set the index type at the top level. In the plot file format, this
; is stored fab by fab, since different components can have different
; types (cell vs. edge centered) within the same plot file. However,
; since amrDescriptor applies to only a single component, it makes more
; sense for the index type to be stored at the top level -- it will be
; the same for every fab in any case.
amrDescriptor.idxtype = (*amrDescriptor.levels[0].fabptr)[0].idxtype

; read the sink particle data
if (nsink ne 0) then begin
	openr, fp, sinkname, /get_lun
	readf, fp, nsink
	sinkdata=dblarr(3*ndim+1)
	sinkptr=0
	for sinkptr=0, nsink-1 do begin
		readf, fp, sinkdata
		amrDescriptor.sinkparticles[sinkptr].m = sinkdata[0]
		amrDescriptor.sinkparticles[sinkptr].x = sinkdata[1:ndim]
		amrDescriptor.sinkparticles[sinkptr].p = $
			sinkdata[ndim+1:2*ndim]
		amrDescriptor.sinkparticles[sinkptr].j = $
			sinkdata[2*ndim+1:3*ndim]
	endfor
	free_lun, fp
endif

; read the star particle data
if (nstar ne 0) then begin
	openr, fp, starname, /get_lun
	readf, fp, nstar
	stardata=dblarr(3*ndim+7)
	starptr=0
	for sinkptr=0, nstar-1 do begin
		readf, fp, stardata
		amrDescriptor.starparticles[sinkptr].m = stardata[0]
		amrDescriptor.starparticles[sinkptr].x = stardata[1:ndim]
		amrDescriptor.starparticles[sinkptr].p = $
			stardata[ndim+1:2*ndim]
		amrDescriptor.starparticles[sinkptr].j = $
			stardata[2*ndim+1:3*ndim]
		amrDescriptor.starparticles[sinkptr].mlast = stardata[3*ndim+1]
		amrDescriptor.starparticles[sinkptr].r = stardata[3*ndim+2]
		amrDescriptor.starparticles[sinkptr].mdeut = stardata[3*ndim+3]
		amrDescriptor.starparticles[sinkptr].n = stardata[3*ndim+4]
		amrDescriptor.starparticles[sinkptr].mdot = stardata[3*ndim+5]
		amrDescriptor.starparticles[sinkptr].burnstate = $
			stardata[3*ndim+6]
	endfor
	free_lun, fp
endif

; return data
return, amrDescriptor
end

	verbose=verbose, maxlevel=maxlevel

; Copyright Mark Krumholz (2002)
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel

; figure out which direction we're summing along
if amr.ndim lt 3 then begin
	print, 'Error: amr object must have 3 dimensions.'
	return, -1
endif
if plane eq 0 then begin
	dim1=1
	dim2=2
	dim3=0
endif else if plane eq 1 then begin
	dim1=0
	dim2=2
	dim3=1
endif else if plane eq 2 then begin
	dim1=0
	dim2=1
	dim3=2
endif else begin
	print, 'Error: plane must be 0, 1, or 2'
	return, -1
endelse

; set range of plot range
if not keyword_set(xrange) then xlim=[amr.boxmin[dim1], amr.boxmax[dim1]] $
else xlim=xrange
if not keyword_set(yrange) then ylim=[amr.boxmin[dim2], amr.boxmax[dim2]] $
else ylim=yrange

; if asked to do so, snap range to grid on finest level
if keyword_set(xrange) and (n_elements(snaptogrid) ne 0) then begin
	xlim[0] = floor( (xlim[0] - amr.boxmin[dim1]) / $
			 amr.gridspacing[dim1,maxlevel] ) * $
		  amr.gridspacing[dim1, maxlevel] + $
		  amr.boxmin[dim1]
	xlim[1] = ceil( (xlim[1] - amr.boxmin[dim1]) / $
			 amr.gridspacing[dim1,maxlevel] ) * $
		  amr.gridspacing[dim1, maxlevel] + $
		  amr.boxmin[dim1]
	if keyword_set(verbose) then print, 'New x range = ', xlim
endif
if keyword_set(yrange) and (n_elements(snaptogrid) ne 0) then begin
	ylim[0] = floor( (ylim[0] - amr.boxmin[dim2]) / $
			 amr.gridspacing[dim2,maxlevel] ) * $
		  amr.gridspacing[dim2, maxlevel] + $
		  amr.boxmin[dim2]
	ylim[1] = ceil( (ylim[1] - amr.boxmin[dim2]) / $
			 amr.gridspacing[dim2,maxlevel] ) * $
		  amr.gridspacing[dim2, maxlevel] + $
		  amr.boxmin[dim2]
	if keyword_set(verbose) then print, 'New y range: ', ylim
endif

; set column image size equal to number of cells on maxlevel
; covered by the requested range
imglo = [floor((xlim[0]-amr.boxmin[dim1]) / $
		amr.gridspacing[dim1,maxlevel]), $
	 floor((ylim[0]-amr.boxmin[dim2]) / $
		amr.gridspacing[dim2,maxlevel])]
imghi = [ceil((xlim[1]-amr.boxmin[dim1]) / $
		amr.gridspacing[dim1,maxlevel]), $
	 ceil((ylim[1]-amr.boxmin[dim2]) / $
		amr.gridspacing[dim2,maxlevel])] - 1
img=fltarr(imghi[0]-imglo[0]+1, imghi[1]-imglo[1]+1)

; go through amr structure, filling appropriate values into
; column image box
boxmin = [amr.boxmin[dim1], amr.boxmin[dim2]]
boxmax = [amr.boxmax[dim1], amr.boxmax[dim2]]
for n=maxlevel, 0, -1 do begin

	; set the refinement ratio between this level and the finest level
	if n ne 0 then refratio = amr.refratio[n-1]

	; create a shorthand for gridspacing
	dx = [amr.gridspacing[dim1,n], amr.gridspacing[dim2,n]]

	; figure out the indices on this level corresponding the the
	; physical limits given
	xidxlim = lonarr(2)
	yidxlim = lonarr(2)
	xidxlim[0] = floor((xlim[0] - boxmin[0]) / dx[0])
	xidxlim[1] = ceil((xlim[1] - boxmin[0]) / dx[0]) - 1
	yidxlim[0] = floor((ylim[0] - boxmin[1]) / dx[1])
	yidxlim[1] = ceil((ylim[1] - boxmin[1]) / dx[1]) - 1

	; loop through fabs
	for m=0, amr.levels[n].nfab-1 do begin

		if keyword_set(verbose) then print, 'Level ', $
			strtrim(string(n),2), ', fab ', strtrim(string(m),2)

		; set up some shorthands
		fabidxmin = (*amr.levels[n].fabptr)[m].idxlo
		fabidxmax = (*amr.levels[n].fabptr)[m].idxhi

		; find the overlap between our target index range and the
		; index range stored in this fab
		overlapmin = lonarr(2)
		overlapmax = lonarr(2)
		overlapmin[0] = (fabidxmin[dim1] gt xidxlim[0]) * $
			fabidxmin[dim1] + $
			(fabidxmin[dim1] le xidxlim[0]) * xidxlim[0]
		overlapmin[1] = (fabidxmin[dim2] gt yidxlim[0]) * $
			fabidxmin[dim2] + $
			(fabidxmin[dim2] le yidxlim[0]) * yidxlim[0]
		overlapmax[0] = (fabidxmax[dim1] le xidxlim[1]) * $
			fabidxmax[dim1] + $
			(fabidxmax[dim1] gt xidxlim[1]) * xidxlim[1]
		overlapmax[1] = (fabidxmax[dim2] le yidxlim[1]) * $
			fabidxmax[dim2] + $
			(fabidxmax[dim2] gt yidxlim[1]) * yidxlim[1]

		; if there is no overlap between this fab and the
		; image box, then move to the next fab
		if (overlapmin[0] gt overlapmax[0]) or $
		   (overlapmin[1] gt overlapmax[1]) then continue

		; We want to get the list of all fine fabs that overlay
		; this current coarse fab. We will store the result in
		; overlay_list. Don't do this if we're on maxlevel,
		; though.
		if n lt maxlevel then begin
		   overlay_list = lonarr(amr.levels[n+1].nfab+1) - 1
		   overlay_list_ptr = 0
		   for i=0, amr.levels[n+1].nfab-1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*amr.levels[n+1].fabptr)[i].idxlo $
			   / refratio
			overlaymax = ((*amr.levels[n+1].fabptr)[i].idxhi+1) $
			   / refratio - 1

			; check if this fine fab overlaps our current fab
			if total(fabidxmin gt overlaymax) ne 0 then continue
			if total(fabidxmax lt overlaymin) ne 0 then continue

			; if we're here, this is an overlaying fab, so
			; record its number
			overlay_list[overlay_list_ptr] = i
			overlay_list_ptr = overlay_list_ptr + 1
		   endfor
		endif

		; grab the data for this region
		data = *(*amr.levels[n].fabptr)[m].dataptr

		; If there isn't an overlay, we can skip this next part.
		; If there is, we construct a mask to block out cells that
		; are overlayed by finer data.
		if n lt maxlevel then begin

		   ; initialize the mask
		   mask = data * 0

		   ; loop through the overlaying fabs
		   overlay_list_ptr=0
		   while overlay_list[overlay_list_ptr] ne -1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*amr.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxlo $
			   / refratio
			overlaymax = ((*amr.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxhi+1) $
			   / refratio - 1

			; create an object to record the intersection limits
			intersectmin = lonarr(amr.ndim)
			intersectmax = lonarr(amr.ndim)

			; loop through dimensions
			for i=0, amr.ndim-1 do begin

			   ; figure out the limits of the intersection region
			   ; in this dimesion
			   intersectmin[i] = max([fabidxmin[i], overlaymin[i]])
			   intersectmax[i] = min([fabidxmax[i], overlaymax[i]])

			endfor

			; convert to mask / data indices
			maskmin = intersectmin - fabidxmin
			maskmax = intersectmax - fabidxmin

			; add 1 to every mask cell for each dimension where
			; that cell is within the intersection limits
			if amr.ndim eq 1 then begin
			   mask[maskmin[0]:maskmax[0]] = $
				mask[maskmin[0]:maskmax[0]] + 1
			endif
			if amr.ndim eq 2 then begin
			   mask[maskmin[0]:maskmax[0],*] = $
				mask[maskmin[0]:maskmax[0],*] + 1
			   mask[*,maskmin[1]:maskmax[1]] = $
				mask[*,maskmin[1]:maskmax[1]] + 1
			endif
			if amr.ndim eq 3 then begin
			   mask[maskmin[0]:maskmax[0],*,*] = $
				mask[maskmin[0]:maskmax[0],*,*] + 1
			   mask[*,maskmin[1]:maskmax[1],*] = $
				mask[*,maskmin[1]:maskmax[1],*] + 1
			   mask[*,*,maskmin[2]:maskmax[2]] = $
				mask[*,*,maskmin[2]:maskmax[2]] + 1
			endif

			; put a 1 in mask cells that are inside the overlap
			; region in every dimension, a 0 otherwise
			mask = (mask eq amr.ndim)

			; apply the mask to the region
			data = (1 - mask) * data

			; increment the pointer
			overlay_list_ptr = overlay_list_ptr + 1

		   endwhile
		endif

		; extract the portion of the data that intersects the
		; requested image box
		xidxlist = lindgen(overlapmax[0]-overlapmin[0]+1) $
			+ overlapmin[0]
		yidxlist = lindgen(overlapmax[1]-overlapmin[1]+1) $
			+ overlapmin[1]
		zidxlist = lindgen(fabidxmax[dim3]-fabidxmin[dim3]+1) $
			+ fabidxmin[dim3]
		datasz = size(data)
		data = reform(data, 1, datasz[1], datasz[2], datasz[3], $
			      /overwrite)
		if plane eq 0 then begin
			data =  reform(data[0, zidxlist-fabidxmin[dim3], $
				xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(zidxlist), $
				n_elements(xidxlist), $
				n_elements(yidxlist))
		endif else if plane eq 1 then begin
			data =  reform(data[0, xidxlist-fabidxmin[dim1], $
				zidxlist-fabidxmin[dim3], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(xidxlist), $
				n_elements(zidxlist), $
				n_elements(yidxlist))
		endif else begin
			data =  reform(data[0, xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2], $
				zidxlist-fabidxmin[dim3]], $
				n_elements(xidxlist), $
				n_elements(yidxlist), $
				n_elements(zidxlist))
		endelse

		; now sum in the appropriate direction
		column = reform(amr.gridspacing[dim3,n] * $
				total(data, dim3 + 1, /double), $
			        n_elements(xidxlist), n_elements(yidxlist))

		; refine the overlap index range to maxlevel
		overlapminref = overlapmin
		overlapmaxref = overlapmax + 1
		for l=maxlevel-1,n,-1 do begin
		   overlapminref = overlapminref * amr.refratio[l]
		   overlapmaxref = overlapmaxref * amr.refratio[l]
		endfor
		overlapmaxref = overlapmaxref - 1

		imgidxlo = long(overlapminref-imglo)
		imgidxhi = long(overlapmaxref-imglo)

		; crop the image box index range to fit what we
		; have available
		imgidxlo = (imgidxlo gt 0) * imgidxlo
		imgidxhi = (imgidxhi le (imghi-imglo)) * imgidxhi + $
			   (imgidxhi gt (imghi-imglo)) * $
				(imghi-imglo)

		; add data to image array
		img[imgidxlo[0]:imgidxhi[0], $
		    imgidxlo[1]:imgidxhi[1]] = $
		   img[imgidxlo[0]:imgidxhi[0], $
		       imgidxlo[1]:imgidxhi[1]] + $
			congrid(column, imgidxhi[0]-imgidxlo[0]+1, $
				imgidxhi[1]-imgidxlo[1]+1)

	endfor
endfor

; return the column image box
return, img

end


; Procedure to compute the cellwise angular momentum of amr data about
; an arbitrary point. The keyword comp specifies the component of j to
; calculate. Setting comp=0 gives jx, comp=1 gives jy, comp=2 gives
; jz. The component is always returned in its proper slot -- the other
; slots are unaltered if comp is set. If comp is not set, all three
; components are calculated. The position to compute the angular
; momentum about is specified by pos. If it is unset, the origin
; is used.

; Read keywords
if not keyword_set(pos) then pos=fltarr(px.ndim)
if n_elements(comp) eq 0 then comp=-1 else begin
    if ((comp lt 0) or (comp ge amr.ndim)) then begin
	print, 'Error: comp must be in the range 0 to ', amr.ndim-1
	return
    endif
endelse

; Create x, y, and z vectors
if (comp ne 0) then xvec=make_rad_amr(px, comp=0, pos=pos)
if (comp ne 1) then yvec=make_rad_amr(px, comp=1, pos=pos)
if (comp ne 2) then zvec=make_rad_amr(px, comp=2, pos=pos)

; Compute the components
if ((comp eq -1) or (comp eq 0)) then begin
	rypz = amr_multiply(yvec, pz)
	rzpy = amr_multiply(zvec, py)
	jx = amr_subtract(rypz, rzpy)
	amr_free, rypz
	amr_free, rzpy
endif
if ((comp eq -1) or (comp eq 1)) then begin
	rzpx = amr_multiply(zvec, px)
	rxpz = amr_multiply(xvec, pz)
	jy = amr_subtract(rzpx, rxpz)
	amr_free, rzpx
	amr_free, rxpz
endif
if ((comp eq -1) or (comp eq 2)) then begin
	rxpy = amr_multiply(xvec, py)
	rypx = amr_multiply(yvec, px)
	jz = amr_subtract(rxpy, rypx)
	amr_free, rxpy
	amr_free, rypx
endif

; Free memory
if (comp ne 0) then amr_free, xvec
if (comp ne 1) then amr_free, yvec
if (comp ne 2) then amr_free, zvec

return
end
	verbose=verbose, maxlevel=maxlevel, nvelbin=nvelbin,	$
	velbinlim=velbinlim, velres=velres, rhocut=rhocut,	$
	velbinctr=velbinctr, smearsize=smearsize,		$
	smeardist=smeardist, csound=csound, zrange=zrange

; Copyright Mark Krumholz (2002)
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

; This routine observes an amr object along the direction specified by
; dir and returns ab observation in the form of a ppv (position-position
; -velocity) cube. The user can specify the number of velocity bins (via
; the nvelbin keyword), can explicitly specify a series of velocity bins
; (by setting velbinlim equal to an array containing the bin limits), and
; can specify the minimum density rho for which the observation will
; see anything (corresponding to an excitation threshhold). All
; observations are assumed to be in the optically thin limit. If velbinlim
; is not set, it returns the velocity bin limits used. The keyword
; velbinctr returns the centers of the velocity bins.

; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=rho.maxlevel
if not keyword_set(csound) then csound=0.0

; figure out which direction we're observing along
if rho.ndim lt 3 then begin
	print, 'Error: rho object must have 3 dimensions.'
	return, -1
endif
if dir eq 0 then begin
	dim1=1
	dim2=2
	dim3=0
endif else if dir eq 1 then begin
	dim1=0
	dim2=2
	dim3=1
endif else if dir eq 2 then begin
	dim1=0
	dim2=1
	dim3=2
endif else begin
	print, 'Error: direction must be 0, 1, or 2'
	return, -1
endelse

; set up the velocity bins
if n_elements(velbinlim) lt 2 then begin
	velmax = max_amr(v, velmin)
	velmax = velmax + 0.001*(velmax-velmin) + 3*csound
	velmin = velmin - 0.001*(velmax-velmin) - 3*csound
	if not keyword_set(velres) then begin
		if not keyword_set(nvelbin) then nvelbin = 20
		velbinlim = velmin + dindgen(nvelbin+1)*(velmax-velmin)/nvelbin
	endif else begin
		nvelbin = ceil((velmax-velmin)/velres)
		velbinlim = velmin + dindgen(nvelbin+1)*velres
	endelse
endif else nvelbin=n_elements(velbinlim)-1
velbinctr=0.5*(velbinlim[0:nvelbin-1]+velbinlim[1:nvelbin])

; set range we're going to extract
if not keyword_set(xrange) then xlim=[rho.boxmin[dim1], rho.boxmax[dim1]] $
else xlim=xrange
if not keyword_set(yrange) then ylim=[rho.boxmin[dim2], rho.boxmax[dim2]] $
else ylim=yrange
if not keyword_set(zrange) then zlim=[rho.boxmin[dim3], rho.boxmax[dim3]] $
else zlim=zrange

; if asked to do so, snap range to grid on finest level
if keyword_set(xrange) and (n_elements(snaptogrid) ne 0) then begin
	xlim[0] = floor( (xlim[0] - rho.boxmin[dim1]) / $
			 rho.gridspacing[dim1,maxlevel] ) * $
		  rho.gridspacing[dim1, maxlevel] + $
		  rho.boxmin[dim1]
	xlim[1] = ceil( (xlim[1] - rho.boxmin[dim1]) / $
			 rho.gridspacing[dim1,maxlevel] ) * $
		  rho.gridspacing[dim1, maxlevel] + $
		  rho.boxmin[dim1]
	if keyword_set(verbose) then print, 'New x range = ', xlim
endif
if keyword_set(yrange) and (n_elements(snaptogrid) ne 0) then begin
	ylim[0] = floor( (ylim[0] - rho.boxmin[dim2]) / $
			 rho.gridspacing[dim2,maxlevel] ) * $
		  rho.gridspacing[dim2, maxlevel] + $
		  rho.boxmin[dim2]
	ylim[1] = ceil( (ylim[1] - rho.boxmin[dim2]) / $
			 rho.gridspacing[dim2,maxlevel] ) * $
		  rho.gridspacing[dim2, maxlevel] + $
		  rho.boxmin[dim2]
	if keyword_set(verbose) then print, 'New y range: ', ylim
endif
if keyword_set(zrange) and (n_elements(snaptogrid) ne 0) then begin
	zlim[0] = floor( (zlim[0] - rho.boxmin[dim3]) / $
			 rho.gridspacing[dim3,maxlevel] ) * $
		  rho.gridspacing[dim3, maxlevel] + $
		  rho.boxmin[dim3]
	zlim[1] = ceil( (zlim[1] - rho.boxmin[dim3]) / $
			 rho.gridspacing[dim3,maxlevel] ) * $
		  rho.gridspacing[dim3, maxlevel] + $
		  rho.boxmin[dim3]
	if keyword_set(verbose) then print, 'New z range: ', zlim
endif

; set column image size equal to number of cells on maxlevel
; covered by the requested range
collo = [floor((xlim[0]-rho.boxmin[dim1]) / $
		rho.gridspacing[dim1,maxlevel]), $
	 floor((ylim[0]-rho.boxmin[dim2]) / $
		rho.gridspacing[dim2,maxlevel])]
colhi = [ceil((xlim[1]-rho.boxmin[dim1]) / $
		rho.gridspacing[dim1,maxlevel]), $
	 ceil((ylim[1]-rho.boxmin[dim2]) / $
		rho.gridspacing[dim2,maxlevel])] - 1
ppv=fltarr(colhi[0]-collo[0]+1, colhi[1]-collo[1]+1, nvelbin)

; go through amr structure, filling appropriate values into
; column image box
boxmin = [rho.boxmin[dim1], rho.boxmin[dim2], rho.boxmin[dim3]]
boxmax = [rho.boxmax[dim1], rho.boxmax[dim2], rho.boxmax[dim3]]
for n=maxlevel, 0, -1 do begin

	; set the refinement ratio between this level and the finest level
	if n ne 0 then refratio = rho.refratio[n-1]

	; create a shorthand for gridspacing
	dx = [rho.gridspacing[dim1,n], rho.gridspacing[dim2,n], $
	      rho.gridspacing[dim3,n]]

	; figure out the indices on this level corresponding the the
	; physical limits given
	xidxlim = lonarr(2)
	yidxlim = lonarr(2)
	zidxlim = lonarr(2)
	xidxlim[0] = floor((xlim[0] - boxmin[0]) / dx[0])
	xidxlim[1] = ceil((xlim[1] - boxmin[0]) / dx[0]) - 1
	yidxlim[0] = floor((ylim[0] - boxmin[1]) / dx[1])
	yidxlim[1] = ceil((ylim[1] - boxmin[1]) / dx[1]) - 1
	zidxlim[0] = floor((zlim[0] - boxmin[2]) / dx[2])
	zidxlim[1] = ceil((zlim[1] - boxmin[2]) / dx[2]) - 1

	; loop through fabs
	for m=0, rho.levels[n].nfab-1 do begin

		if keyword_set(verbose) then print, 'Level ', $
			strtrim(string(n),2), ', fab ', strtrim(string(m),2)

		; set up some shorthands
		fabidxmin = (*rho.levels[n].fabptr)[m].idxlo
		fabidxmax = (*rho.levels[n].fabptr)[m].idxhi

		; find the overlap between our target index range and the
		; index range stored in this fab
		overlapmin = lonarr(3)
		overlapmax = lonarr(3)
		overlapmin[0] = (fabidxmin[dim1] gt xidxlim[0]) * $
			fabidxmin[dim1] + $
			(fabidxmin[dim1] le xidxlim[0]) * xidxlim[0]
		overlapmin[1] = (fabidxmin[dim2] gt yidxlim[0]) * $
			fabidxmin[dim2] + $
			(fabidxmin[dim2] le yidxlim[0]) * yidxlim[0]
		overlapmin[2] = (fabidxmin[dim3] gt zidxlim[0]) * $
			fabidxmin[dim3] + $
			(fabidxmin[dim3] le zidxlim[0]) * zidxlim[0]
		overlapmax[0] = (fabidxmax[dim1] le xidxlim[1]) * $
			fabidxmax[dim1] + $
			(fabidxmax[dim1] gt xidxlim[1]) * xidxlim[1]
		overlapmax[1] = (fabidxmax[dim2] le yidxlim[1]) * $
			fabidxmax[dim2] + $
			(fabidxmax[dim2] gt yidxlim[1]) * yidxlim[1]
		overlapmax[2] = (fabidxmax[dim3] le zidxlim[1]) * $
			fabidxmax[dim3] + $
			(fabidxmax[dim3] gt zidxlim[1]) * zidxlim[1]

		; if there is no overlap between this fab and the
		; image box, then move to the next fab
		if (overlapmin[0] gt overlapmax[0]) or $
		   (overlapmin[1] gt overlapmax[1]) or $
		   (overlapmin[2] gt overlapmax[2]) then continue

		; We want to get the list of all fine fabs that overlay
		; this current coarse fab. We will store the result in
		; overlay_list. Don't do this if we're on maxlevel,
		; though.
		if n lt maxlevel then begin
		   overlay_list = lonarr(rho.levels[n+1].nfab+1) - 1
		   overlay_list_ptr = 0
		   for i=0, rho.levels[n+1].nfab-1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*rho.levels[n+1].fabptr)[i].idxlo $
			   / refratio
			overlaymax = ((*rho.levels[n+1].fabptr)[i].idxhi+1) $
			   / refratio - 1

			; check if this fine fab overlaps our current fab
			if total(fabidxmin gt overlaymax) ne 0 then continue
			if total(fabidxmax lt overlaymin) ne 0 then continue

			; if we're here, this is an overlaying fab, so
			; record its number
			overlay_list[overlay_list_ptr] = i
			overlay_list_ptr = overlay_list_ptr + 1
		   endfor
		endif

		; grab the data for this region
		rhofab = *(*rho.levels[n].fabptr)[m].dataptr
		vfab = *(*v.levels[n].fabptr)[m].dataptr

		; If there isn't an overlay, we can skip this next part.
		; If there is, we construct a mask to block out cells that
		; are overlayed by finer data.
		if n lt maxlevel then begin

		   ; initialize the mask
		   mask = rhofab * 0

		   ; loop through the overlaying fabs
		   overlay_list_ptr=0
		   while overlay_list[overlay_list_ptr] ne -1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*rho.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxlo $
			   / refratio
			overlaymax = ((*rho.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxhi+1) $
			   / refratio - 1

			; create an object to record the intersection limits
			intersectmin = lonarr(rho.ndim)
			intersectmax = lonarr(rho.ndim)

			; loop through dimensions
			for i=0, rho.ndim-1 do begin

			   ; figure out the limits of the intersection region
			   ; in this dimesion
			   intersectmin[i] = max([fabidxmin[i], overlaymin[i]])
			   intersectmax[i] = min([fabidxmax[i], overlaymax[i]])

			endfor

			; convert to mask / data indices
			maskmin = intersectmin - fabidxmin
			maskmax = intersectmax - fabidxmin

			; add 1 to every mask cell for each dimension where
			; that cell is within the intersection limits
			if rho.ndim eq 1 then begin
			   mask[maskmin[0]:maskmax[0]] = $
				mask[maskmin[0]:maskmax[0]] + 1
			endif
			if rho.ndim eq 2 then begin
			   mask[maskmin[0]:maskmax[0],*] = $
				mask[maskmin[0]:maskmax[0],*] + 1
			   mask[*,maskmin[1]:maskmax[1]] = $
				mask[*,maskmin[1]:maskmax[1]] + 1
			endif
			if rho.ndim eq 3 then begin
			   mask[maskmin[0]:maskmax[0],*,*] = $
				mask[maskmin[0]:maskmax[0],*,*] + 1
			   mask[*,maskmin[1]:maskmax[1],*] = $
				mask[*,maskmin[1]:maskmax[1],*] + 1
			   mask[*,*,maskmin[2]:maskmax[2]] = $
				mask[*,*,maskmin[2]:maskmax[2]] + 1
			endif

			; put a 1 in mask cells that are inside the overlap
			; region in every dimension, a 0 otherwise
			mask = (mask eq rho.ndim)

			; apply the mask to the region
			rhofab = (1 - mask) * rhofab

			; increment the pointer
			overlay_list_ptr = overlay_list_ptr + 1

		   endwhile
		endif

		; extract the portion of the data that intersects the
		; requested image box
		xidxlist = lindgen(overlapmax[0]-overlapmin[0]+1) $
			+ overlapmin[0]
		yidxlist = lindgen(overlapmax[1]-overlapmin[1]+1) $
			+ overlapmin[1]
		zidxlist = lindgen(overlapmax[2]-overlapmin[2]+1) $
			+ overlapmin[2]
		datasz = size(rhofab)
		rhofab = reform(rhofab, 1, datasz[1], datasz[2], datasz[3], $
			      /overwrite)
		vfab = reform(vfab, 1, datasz[1], datasz[2], datasz[3], $
			      /overwrite)
		if dir eq 0 then begin
			rhofab = reform(rhofab[0, zidxlist-fabidxmin[dim3], $
				xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(zidxlist), $
				n_elements(xidxlist), $
				n_elements(yidxlist))
			vfab = reform(vfab[0, zidxlist-fabidxmin[dim3], $
				xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(zidxlist), $
				n_elements(xidxlist), $
				n_elements(yidxlist))
		endif else if dir eq 1 then begin
			rhofab = reform(rhofab[0, xidxlist-fabidxmin[dim1], $
				zidxlist-fabidxmin[dim3], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(xidxlist), $
				n_elements(zidxlist), $
				n_elements(yidxlist))
			vfab = reform(vfab[0, xidxlist-fabidxmin[dim1], $
				zidxlist-fabidxmin[dim3], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(xidxlist), $
				n_elements(zidxlist), $
				n_elements(yidxlist))
		endif else begin
			rhofab = reform(rhofab[0, xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2], $
				zidxlist-fabidxmin[dim3]], $
				n_elements(xidxlist), $
				n_elements(yidxlist), $
				n_elements(zidxlist))
			vfab = reform(vfab[0, xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2], $
				zidxlist-fabidxmin[dim3]], $
				n_elements(xidxlist), $
				n_elements(yidxlist), $
				n_elements(zidxlist))
		endelse

		; set densities below the cutoff to zero
		if keyword_set(rhocut) then $
			rhofab = rhofab * (rhofab gt rhocut)

		; refine the overlap index range to maxlevel
		overlapminref = overlapmin
		overlapmaxref = overlapmax + 1
		for l=maxlevel-1,n,-1 do begin
		   overlapminref = overlapminref * rho.refratio[l]
		   overlapmaxref = overlapmaxref * rho.refratio[l]
		endfor
		overlapmaxref = overlapmaxref - 1

		imgidxlo = long(overlapminref[0:1]-collo)
		imgidxhi = long(overlapmaxref[0:1]-collo)

		; crop the image box index range to fit what we
		; have available
		imgidxlo = (imgidxlo gt 0) * imgidxlo
		imgidxhi = (imgidxhi le (colhi-collo)) * imgidxhi + $
			   (imgidxhi gt (colhi-collo)) * $
				(colhi-collo)

		; loop through velocity bins
		for i=0, nvelbin-1 do begin

			; create a mask for this bin
			vsize=size(vfab)
			if not keyword_set(csound) then begin
				; just bin assuming no maxwellian distrib.
				vmask = (vfab ge velbinlim[i]) and $
					(vfab lt velbinlim[i+1])
			endif else begin
				; bin by maxwellian distribution
				vmat1 = (velbinlim[i]-vfab) / $
					(sqrt(2.)*csound)
				vmat2 = (velbinlim[i+1]-vfab) / $
					(sqrt(2.)*csound)
				vmask = 0.5 * (errorf(vmat2)-errorf(vmat1))
			endelse
			vmask=reform(vmask, vsize[1], vsize[2], vsize[3])

			; now sum in the appropriate direction
			rhofabtmp=reform(rhofab*vmask, vsize[1], vsize[2], $
					 vsize[3])
			column = reform(rho.gridspacing[dim3,n] * $
					total(rhofabtmp, dim3 + 1, $
					/double), $
				        n_elements(xidxlist), $
					n_elements(yidxlist))

			; add data to ppv array
			ppv[imgidxlo[0]:imgidxhi[0], $
		    	    imgidxlo[1]:imgidxhi[1], i] = $
		 	  ppv[imgidxlo[0]:imgidxhi[0], $
		    	      imgidxlo[1]:imgidxhi[1], i] + $
				congrid(column, imgidxhi[0]-imgidxlo[0]+1, $
					imgidxhi[1]-imgidxlo[1]+1)

		endfor

	endfor
endfor

; smear the data with a beam if requested
if keyword_set(smearsize) then begin
	for n=0,nvelbin-1 do begin
		ppv[*,*,n] = beamsmear(ppv[*,*,n], $
				       rho.gridspacing[0,maxlevel], $
				       smeardist, smearsize)
	endfor
endif

; return the ppv data
return, ppv

end

