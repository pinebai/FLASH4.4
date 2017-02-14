function get_amrcomponent_flash, filename, componentName, verbose = verbose

;=================================================
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
  if (itype EQ 2) then pring, 'Recognized NCDF file format.'
endif

;------------------------------------------------------------
; read in FLASH data using code from FIDLR program read_amr.
;------------------------------------------------------------

if(itype EQ 1) then begin

  if (verbose) then begin
    print, 'Opening file for input...'
  endif
                                                                                
;------------------------------------------------------------------
; open up the file for read now and read in the header information
;------------------------------------------------------------------                                                                             
;--------------------------------------------------------------
; Read in FLASH annotations and print
;--------------------------------------------------------------

  file_information,filename,FLASH_VERSION=flashversion

; do full read_amr

  read_amr,filename,VAR_NAME=componentName,$
           TREE=tree,DATA=unk,PARAMETERS=params,STORED_VARS=unk_names,$
           NUM_PARTICLES=numParticle,PARTICLES=particles,INT_PROP_NAMES=IntPropNames,$
           REAL_PROP_NAMES=realPropNames,GEOMETRY=geometry
                                                                                

;-------------------------------------------------------------------------
; LIMITATION
; Assume Cartesian geometry, and determine the problem domain box extent
;  by inspecting the bounding box for all blocks.
;-------------------------------------------------------------------------

  boxmin = fltarr (params.ndim)
  boxmax = fltarr (params.ndim)
  ndim = params.ndim

  for i = 0, ndim - 1 do begin
    boxmin [i] = min (tree [*].bndBox [0, i] )
    boxmax [i] = max (tree [*].bndBox [1, i] )
  endfor

; LIMITATION - We read in the particles, but do not store them yet.        
                                                                              
                                                                                
endif  ; end if itype = 1

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

periodicity=intarr(ndim,maxlevel+1)
periodicity (*, *) = 0
;----------------------------------------------------
; Set refinement ratio for datasets with refinement
;----------------------------------------------------

if maxlevel gt 0 then refratio=lonarr(maxlevel+1) else refratio=0

; Because PARAMESH imposes equal block sizes over all levels, we can
; compute the refinement ratio from the spatial extent of any two boxes
; on two different levels in any direction.

refratio [0] = 1

; Compute refinement ratio across all levels.

for n = 1, maxlevel do begin

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

nb = [params.nxb, params.nyb, params.nzb]

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
		time:params.time, $
		maxlevel:maxlevel, $
		boxmin:boxmin, boxmax:boxmax, $
		refratio:refratio, $
		idxlo:idxlo, idxhi:idxhi, $
		idxtype:intarr (ndim), $
		periodicity:periodicity, $
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
		periodicity:periodicity, $
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
                                                                                
if (bigendian ne is_ieee_big()) then BEGIN
    do_swap_endian = 1 
ENDIF else BEGIN
    do_swap_endian = 0
ENDELSE 

;---------------------------------------------------------
;  convert data from FLASH "tree" to amrDescriptor format
;---------------------------------------------------------

for n=0, maxlevel do begin

	; Fill out level descriptors

	amrDescriptor.levels[n].level=n
	amrDescriptor.levels[n].ngrids = $
                (size (where(tree[*].lrefine EQ n + 1) )) [1]
	amrDescriptor.levels[n].leveltime = $
		params.time
	amrDescriptor.levels[n].levelsteps = $
	        nsteps	
	amrDescriptor.levels[n].idxlo = $
		idxlo (*, n)
	amrDescriptor.levels[n].idxhi = $
		idxhi (*, n)
;	amrDescriptor.levels[n].periodicity = $
;		0
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
