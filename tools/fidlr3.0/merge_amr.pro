;------------------------------------------------------------------------------
; merge_amr.pro
;
; take a hydro variable in AMR block format and merge it into a
; single, uniformly gridded array at the desired AMR resolution.
;
; this version is for 2-d or 3-d data.  The dimensionality of the
; dataset is determined from the PARAMETERS structure that is passed
; in along with the dataset.
;
; corner data is supported automatically, if params.corners is set to
; 1.  The data will be interpolated to the cell centers, before
; scaling to a uniform grid.
;
; arguments:
;
;     tvar --  the variable to be scaled, tvar[maxblocks,nxb,nyb
;
;     TREE -- the tree information structure, as returned from the
;             read_amr routine.
;
;     PARAMETERS -- the parameter structure, as returned from the
;                   read_amr routine.
;
; optional arguments
;
;     XRANGE -- the minimum and maximum x coord of the desired uniform
;               domain 
;
;     YRANGE -- the minimum and maximum y coord of the desired uniform
;               domain 
;
;     ZRANGE -- the minimum and maximum z coord of the desired uniform
;               domain (3-D only).
;
;    Note, for 3-d data, if the range in any 1 coordinate direction is
;    set equal to each other, then a 2-d plane will be produced,
;    passing through that coord point.  
;
;    Ex: XRANGE=[0.,1.], YRANGE=[0.,1.], ZRANGE=[0.5,0.5]
; 
;      will merge the data onto a 2-d plane in the x-y plane passing
;      through z=0.5.
; 
;  
;     XMERGE  \  return the x, y, and z (3-d only) positions of the cells 
;     YMERGE   > in the new master scaled array
;     ZMERGE  /
;
;     sample  -- # of levels down from finest resolution to sample
;                into the new array.  Blocks with finer resolution
;                will be averaged.  Note: sample = 0 is the same as
;                omitting the sample keyword.
;
;     DOUBLE -- use double precision -- this requires that the AMR
;               data be double precision, usually the result of
;               passing /DOUBLE to the read_amr routine.
;
;
; This routine can be used to select a subset of the domain for
; uniform gridding.  Only this region is allocated as a uniform box.
; If the minimum and maximum range in a direction are equal, then the
; data is sampled into a plane passing through that coordinate value.
;
; Description of variables:
;
; scaling is the factor that a block is expanded (or shruck) by to put
; it at the same resolution as the scaled domain.  Coarsely refined
; blocks need to be expanded greatly when mapping into the uniform domain
; 
; {x,y,z}range_min_index is the offset, in units of uniform zones from
; the computational zone boundary, and the boundary of the domain we
; are scaling.
;
; {x,y,z}start are the location into the scaled block to copy.  In
; most cases, this will be 0, as the beginning of the block normally
; falls into the scaled domain.  If the domain intersects a block,
; then this will be 0 < {x,y,z}start < scaled*n{x,y,z}b - 1
;
; {x,y,z}index is the offset from the block lower boundary to the
; lower boundary of the computation domain
;
; 
;     +--------------------------------------------------------+
;     |                                                        |
;     |                                                        |
;     |       ................................                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .      xxxxxxxx                .                 |
;     |       .      x      x                .                 |
;     |       .      x      x                .                 |
;     |       .      x      x                .                 |
;     |       .      x      x                .                 |
;     |       .      xxxxxxxx                .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       ................................                 |
;     |                                                        |
;     |                                                        |
;     +--------------------------------------------------------+
;     
;     |<- A ->|
;             |<-B ->|
;
;
;    |  marks the computation domain boundary
;    .  marks the region we are scaling to a uniform mesh
;    x  marks the current block
;
;    A  is xrange_min_index
;    B  is xindex
;
;
; The case is more complicated when the block goes out of the domain
; we are scaling
;
;
;     +--------------------------------------------------------+
;     |                                                        |
;     |                                                        |
;     |       ................................                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |    xxxxxxxxxxxxxxxxxx                .                 |
;     |    x  .             x                .                 |
;     |    x  .             x                .                 |
;     |    x  .             x                .                 |
;     |    x  .             x                .                 |
;     |    xxxxxxxxxxxxxxxxxx                .                 |
;     |       .                              .                 |
;     |       .                              .                 |
;     |       ................................                 |
;     |                                                        |
;     |                                                        |
;     +--------------------------------------------------------+
;
;     |<- A ->|
;        ->| B|<-
;
; now, xindex (B) is negative
; and xstart is -B -- we start in the center of the block when copying
;                     data into the uniform grid
; 
; 
;------------------------------------------------------------------------------

function merge_amr, temp_arr, $
                    TREE = tree, PARAMETERS = params, $
                    DOUBLE = double, $
                    XRANGE = txrange, YRANGE = tyrange, ZRANGE = tzrange, $
                    XMERGE = xout, YMERGE = yout, ZMERGE = zout, $
                    SAMPLE = sample, DEBUG = debug, REFORM=reformit

; allow flag for reforming to be undefined
if (not Keyword_Set(reformit)) then reformit = 0

start_time = systime(/seconds)

; find the minimum and maximum refinement levels of the good data
index_good = where(tree[*].nodetype EQ 1)

lmax = max(tree[index_good].lrefine)
lmin = min(tree[index_good].lrefine)

if n_elements(sample) EQ 0 then begin
  lwant = lmax
endif else begin
  lwant = lmax - sample
  if lwant LE 0 then begin
    print, 'Warning in merge_amr, cannot sample below coarsest resolution'
    lwant = 1
  endif
endelse
    
if n_elements(double) EQ 0 then double = 0

; if no domain ranges are specified, then set them to the extreme of
; the domain    
if n_elements(txrange) NE 2 then begin
  xrange = fltarr(2)
  xrange[0] = min(tree[*].bndBox[0,0])     ; set to minimum of x coord
  xrange[1] = max(tree[*].bndBox[1,0])     ; set to maximum of x coord
endif else begin
  xrange = float(txrange)
endelse

if n_elements(tyrange) NE 2 then begin
  yrange = fltarr(2)
  yrange[0] = min(tree[*].bndBox[0,1])     ; set to minimum of y coord
  yrange[1] = max(tree[*].bndBox[1,1])     ; set to maximum of y coord
endif else begin
  yrange = float(tyrange)
endelse

if (params.ndim EQ 3) then begin
  if n_elements(tzrange) NE 2 then begin
    zrange = fltarr(2)
    zrange[0] = min(tree[*].bndBox[0,2]) ; set to minimum of z coord
    zrange[1] = max(tree[*].bndBox[1,2]) ; set to maximum of z coord
  endif else begin
    zrange = float(tzrange)
  endelse
endif


if params.corners EQ 1 then begin
  if Keyword_Set(debug) then print, 'corner data detected, interpolating to cell centers'

  nxb_centered = params.nxb - 1
  nyb_centered = params.nyb - 1
  if (params.ndim EQ 3) then nzb_centered = params.nzb - 1
endif

if Keyword_Set(debug) then begin
  print, 'zooming between, x: ', xrange[0], xrange[1]
  print, '                 y: ', yrange[0], yrange[1]
  if (params.ndim EQ 3) then print, '                 z: ', yrange[0], yrange[1]
endif 


;------------------------------------------------------------------------------
; compute the number of zones in each direction for the user's domain
; limits, assuming uniform grid
;------------------------------------------------------------------------------
case params.corners of

  0: begin   ; no corners
    dx_fine = (max(tree[*].bndBox[1,0]) - $
               min(tree[*].bndBox[0,0]))/ $
      (params.ntopx*params.nxb*2^(lwant-1))

    if xrange[1] EQ xrange[0] then xrange[1] = xrange[0] + dx_fine
    dx_merge = xrange[1] - xrange[0]


    dy_fine = (max(tree[*].bndBox[1,1]) - $
               min(tree[*].bndBox[0,1]))/ $
      (params.ntopy*params.nyb*2^(lwant-1))

    if yrange[1] EQ yrange[0] then yrange[1] = yrange[0] + dy_fine
    dy_merge = yrange[1] - yrange[0]

        
    if (params.ndim EQ 3) then begin
      dz_fine = (max(tree[*].bndBox[1,2]) - $
                 min(tree[*].bndBox[0,2]))/ $
        (params.ntopz*params.nzb*2^(lwant-1))

      if zrange[1] EQ zrange[0] then zrange[1] = zrange[0] + dz_fine
      dz_merge = zrange[1] - zrange[0]
    endif
  end

  1: begin  ; corners
    dx_fine = (max(tree[*].bndBox[1,0]) - $
               min(tree[*].bndBox[0,0]))/ $
      (params.ntopx*nxb_centered*2^(lwant-1))

    if xrange[1] EQ xrange[0] then xrange[1] = xrange[0] + dx_fine
    dx_merge = xrange[1] - xrange[0]


    dy_fine = (max(tree[*].bndBox[1,1]) - $
               min(tree[*].bndBox[0,1]))/ $
      (params.ntopy*nyb_centered*2^(lwant-1))

    if yrange[1] EQ yrange[0] then yrange[1] = yrange[0] + dy_fine
    dy_merge = yrange[1] - yrange[0]

        
    if (params.ndim EQ 3) then begin
      dz_fine = (max(tree[*].bndBox[1,2]) - $
                 min(tree[*].bndBox[0,2]))/ $
        (params.ntopz*nzb_centered*2^(lwant-1))

      if zrange[1] EQ zrange[0] then zrange[1] = zrange[0] + dz_fine
      dz_merge = zrange[1] - zrange[0]
    endif
  end

endcase

nx = long(dx_merge/dx_fine) > 1
ny = long(dy_merge/dy_fine) > 1
if params.ndim EQ 3 then nz = long(dz_merge/dz_fine) > 1


;------------------------------------------------------------------------------
; declare the storage for the uniformly gridded result
;------------------------------------------------------------------------------
case params.ndim of
  2: begin
    if (double) then begin
      temp_merge = dblarr(nx,ny)
    endif else begin
      temp_merge = fltarr(nx,ny)
    endelse

    ; initialize temp_merge to negative numbers, if there are any holes in
    ; the domain (i.e. non-rectangular), then these will remain negative
    ; and scale will color these black
    temp_merge[*,*] = -1.e30

    if keyword_set(debug) then print, 'uniform grid size = ', nx, ny

  end

  3: begin
    if (double) then begin
      temp_merge = dblarr(nx,ny,nz)
    endif else begin
      temp_merge = fltarr(nx,ny,nz)
    endelse

    ; initialize temp_merge to negative numbers
    temp_merge[*,*] = -1.e30

    if keyword_set(debug) then print, 'uniform grid size = ', nx, ny, nz

  end
endcase


;------------------------------------------------------------------------------
; compute the x, y, and z coordinates of the data for the entire domain
;------------------------------------------------------------------------------
xmax = max(tree[*].bndBox[1,0])
xmin = min(tree[*].bndBox[0,0])

ymax = max(tree[*].bndBox[1,1])
ymin = min(tree[*].bndBox[0,1])

if (params.ndim EQ 3) then begin
  zmax = max(tree[*].bndBox[1,2])
  zmin = min(tree[*].bndBox[0,2])
endif

case params.corners of

  0: begin  ; no corners
    nx_domain = params.ntopx*params.nxb*2^(lwant-1)
    ny_domain = params.ntopy*params.nyb*2^(lwant-1)
    if params.ndim EQ 3 then nz_domain = $
      params.ntopz*params.nzb*2^(lwant-1)
  end

  1: begin  ; corners
    nx_domain = params.ntopx*nxb_centered*2^(lwant-1)
    ny_domain = params.ntopy*nyb_centered*2^(lwant-1)
    if params.ndim EQ 3 then nz_domain = $
      params.ntopz*nzb_centered*2^(lwant-1)
  end

endcase


dx = (xmax - xmin)/nx_domain
x = (findgen(nx_domain) + .5)*dx + xmin

dy = (ymax - ymin)/ny_domain
y = (findgen(ny_domain) + .5)*dy + ymin

if (params.ndim EQ 3) then begin
  dz = (zmax - zmin)/nz_domain
  z = (findgen(nz_domain) + .5)*dz + zmin
endif

;------------------------------------------------------------------------------
; find out where the range extreme would fall in a uniform grid of the 
; total domain -- compute these using nx and ny above to avoid errors
; from rounding
;------------------------------------------------------------------------------
xrange_min_index = (where(xrange[0] LE x))[0]
xrange_max_index = xrange_min_index + nx - 1

yrange_min_index = (where(yrange[0] LE y))[0]
yrange_max_index = yrange_min_index + ny - 1

if (params.ndim EQ 3) then begin
  zrange_min_index = (where(zrange[0] LE z))[0]
  zrange_max_index = zrange_min_index + nz - 1
endif


;------------------------------------------------------------------------------
; get the list of blocks that fall in the zoomed domain and are leaf
; blocks
;------------------------------------------------------------------------------
case params.ndim of
  2: plotBlocks = where(((tree[*].nodeType EQ 1 AND $
                          tree[*].lrefine LT lwant) OR $
                         (tree[*].lrefine EQ lwant)) AND $
                        tree[*].bndBox[1,0] GE xrange[0] AND $
                        tree[*].bndBox[0,0] LE xrange[1] AND $
                        tree[*].bndBox[1,1] GE yrange[0] AND $
                        tree[*].bndBox[0,1] LE yrange[1])

  3: plotBlocks = where(tree[*].nodeType EQ 1 AND $
                        tree[*].bndBox[1,0] GE xrange[0] AND $
                        tree[*].bndBox[0,0] LE xrange[1] AND $
                        tree[*].bndBox[1,1] GE yrange[0] AND $
                        tree[*].bndBox[0,1] LE yrange[1] AND $
                        tree[*].bndBox[1,2] GE zrange[0] AND $
                        tree[*].bndBox[0,2] LE zrange[1])
endcase

numPlotBlocks = (size(plotBlocks))[1]


;------------------------------------------------------------------------------
; loop over the blocks and store them at the resolution of the 
; finest block into the master array -- difficulties arise when
; only a part of a block falls into the user's domain subset
;------------------------------------------------------------------------------
i = (lonarr(1))[0]

for i = 0l, numPlotBlocks-1 do begin

  curBlk = plotBlocks[i] 

  ; compute the scale factor to bring it to the resolution of the finest block
  scaling = 2.^(lwant - tree[curBlk].lrefine)

  case params.ndim of
    2: begin
      sub_array = congrid(reform(temp_arr[curBlk,*,*]), $
                          scaling*params.nxb, $
                          scaling*params.nyb)
    end
    3: begin

      case params.corners of
        0: begin        ; no corners
          sub_array = rebin(reform(temp_arr[curBlk,*,*,*]), $
                            scaling*params.nxb, $
                            scaling*params.nyb, $
                            scaling*params.nzb, /SAMPLE)
        end

        1: begin        ; corners

          ; create a temporary array to hold the interpolated block's data
          case double of
            0: var_cc = fltarr(nxb_centered, $
                               nyb_centered, $
                               nzb_centered)
                        
            1: var_cc = dblarr(nxb_centerd, $
                               nyb_centerd, $
                               nzb_centered)
          endcase

          ; get the temp array by averaging to cell centers
          for kz = 0, nzb_centered-1 do begin
            for jz = 0, nyb_centered-1 do begin
              for iz = 0, nxb_centered-1 do begin
                                
                var_cc[iz,jz,kz] = .125* $
                  (temp_arr[curBlk,iz,  jz,  kz  ] + $
                   temp_arr[curBlk,iz+1,jz,  kz  ] + $
                   temp_arr[curBlk,iz,  jz+1,kz  ] + $
                   temp_arr[curBlk,iz+1,jz+1,kz  ] + $
                   temp_arr[curBlk,iz,  jz,  kz+1] + $
                   temp_arr[curBlk,iz+1,jz,  kz+1] + $
                   temp_arr[curBlk,iz,  jz+1,kz+1] + $
                   temp_arr[curBlk,iz+1,jz+1,kz+1])

              endfor
            endfor
          endfor
                    
          sub_array = rebin(reform(var_cc), $
                            scaling*nxb_centered, $
                            scaling*nyb_centered, $
                            scaling*nzb_centered, /SAMPLE)

        end
      endcase
    end
  endcase

  ; find out where it should live in the master array
  xindex = (where(tree[curBlk].bndBox[0,0] LE x))[0] - xrange_min_index
  yindex = (where(tree[curBlk].bndBox[0,1] LE y))[0] - yrange_min_index
  if (params.ndim EQ 3) then $
    zindex = (where(tree[curBlk].bndBox[0,2] LE z))[0] - zrange_min_index
    
  ; set the limits for the local block array
  case params.corners of

    0: begin                ; no corners
                
      xstart = 0
      xspan = scaling*params.nxb-1
            
      ystart = 0
      yspan = scaling*params.nyb-1
            
      if (params.ndim EQ 3) then begin
        zstart = 0
        zspan = scaling*params.nzb-1
      endif
            
    end
        
    1: begin                ; corners
            
      xstart = 0
      xspan = scaling*nxb_centered-1
            
      ystart = 0
      yspan = scaling*nyb_centered-1
            
      if (params.ndim EQ 3) then begin
        zstart = 0
        zspan = scaling*nzb_centered-1
      endif
                
    end
            
  endcase


  ; first check if we go over the maximum 
  if (xindex + xspan GT (xrange_max_index-xrange_min_index)) then $
    xspan = (xrange_max_index-xrange_min_index) - xindex

  ; now check if we go below the minimum
  if (xindex LT 0 AND xindex + xspan GE 0) then begin
    xstart = -xindex
    xspan = xindex+xspan
    xindex = 0
  endif

  if (yindex + yspan GT (yrange_max_index-yrange_min_index)) then $
    yspan = (yrange_max_index-yrange_min_index) - yindex
    
  if (yindex LT 0 AND yindex + yspan GE 0) then begin
    ystart = -yindex
    yspan = yindex+yspan
    yindex = 0
  endif
    
  if (params.ndim EQ 3) then begin
    if (zindex + zspan GT (zrange_max_index-zrange_min_index)) then $
      zspan = (zrange_max_index-zrange_min_index) - zindex
        
    if (zindex LT 0 AND zindex + zspan GE 0) then begin
      zstart = -zindex
      zspan = zindex+zspan
      zindex = 0
    endif
  endif

  case params.ndim of
        
    2: begin
      if (xindex + xspan GE 0 AND $
          yindex + yspan GE 0 AND $
          xindex LE (xrange_max_index-xrange_min_index) AND $
          yindex LE (yrange_max_index-yrange_min_index)) then begin

        ; store it
        temp_merge[xindex:xindex+xspan, $
                   yindex:yindex+yspan] = $
          sub_array[xstart:xstart+xspan, $
                    ystart:ystart+yspan]

                
      endif
    end

    3: begin
      if (xindex + xspan GE 0 AND $
          yindex + yspan GE 0 AND $
          zindex + zspan GE 0 AND $
          xindex LE (xrange_max_index-xrange_min_index) AND $
          yindex LE (yrange_max_index-yrange_min_index) AND $
          zindex LE (zrange_max_index-zrange_min_index)) then begin

        ; store it
        temp_merge[xindex:xindex+xspan, $
                   yindex:yindex+yspan, $
                   zindex:zindex+zspan] = $
          sub_array[xstart:xstart+xspan, $
                    ystart:ystart+yspan, $
                    zstart:zstart+zspan]
      endif
    end
  endcase
endfor

if Keyword_Set(debug) then print, 'time to merge = ', systime(/seconds) - start_time

xout = x[xrange_min_index:xrange_max_index]
yout = y[yrange_min_index:yrange_max_index]  
if params.ndim EQ 3 then zout = z[zrange_min_index:zrange_max_index]


; many routines immediately issue 
;variable = merge_amr(....)
;variable = reformat(temporary(variable)
; the keyword REFORM=reformit does this automatically
if Keyword_Set(reformit) then begin
  temp_merge = reform(temporary(temp_merge))
endif 

return, temp_merge    

end
