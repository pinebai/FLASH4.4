function loaddata, filename, var, $
                   SAMPLE=sample, $
                   DOUBLE=double, $
                   XCOORDS=x, YCOORDS=y, ZCOORDS=z, $
                   XLCOORD=xl, XRCOORD=xr, $
                   XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange, $
                   TIME=time, UNIFORM_1D=uniform_1d

;
; simple proceedure to read a single variable (or derived variable), var,
; from the file, filename.  The data is read in, and put onto a
; uniformly gridded mesh at the resolution of the finest AMR mesh (or 
; subsampled by sample levels if desired).  The coordinates of the
; data are returned through the optional keywords, XCOORDS, YCOORDS,
; and ZCOORDS.
;
; USAGE
; 
; Here is the basic way to read a file in and grab the energy.
;
; energy = loaddata('myCheckPointFile','ener')
; 
; optionally
; 
; energy = loaddata('myCheckPointFile','ener', XCOORDS=x, YCOORDS=y,
; TIME=t)
;



if (n_elements(filename) EQ 0) then begin
    print, 'ERROR: no filename specified to get_var_list'
    return, -1
endif

if n_elements(var) EQ 0 then begin
    print, 'ERROR: no variable specified'
    return, -1
endif

if n_elements(sample) EQ 0 then sample = 0

if n_elements(uniform_1d) EQ 0 then uniform_1d = 0

if n_elements(double) EQ 0 then double = 0

;------------------------------------------------------------------------------
; read in the data
;------------------------------------------------------------------------------
itype = determine_file_type(filename)

if (double) then begin
    read_amr, filename, VAR_NAME=var, $
      TREE=tree, DATA=unk, PARAMETERS=params
endif else begin
    read_amr, filename, VAR_NAME=var, $
      TREE=tree, DATA=unk, PARAMETERS=params
endelse

time = params.time

; set the ranges
if n_elements(xrange) EQ 0 then begin
    xrange = fltarr(2)
    xrange[0] = min(tree[*].bndBox[0,0])     ; set to minimum of x coord
    xrange[1] = max(tree[*].bndBox[1,0])     ; set to maximum of x coord
endif

if (params.ndim GE 2) then begin
    if n_elements(yrange) EQ 0 then begin
        yrange = fltarr(2)
        yrange[0] = min(tree[*].bndBox[0,1]) ; set to minimum of x coord
        yrange[1] = max(tree[*].bndBox[1,1]) ; set to maximum of x coord
    endif
endif

if (params.ndim EQ 3) then begin
    if n_elements(zrange) EQ 0 then begin
        zrange = fltarr(2)
        zrange[0] = min(tree[*].bndBox[0,2]) ; set to minimum of x coord
        zrange[1] = max(tree[*].bndBox[1,2]) ; set to maximum of x coord
    endif
endif
   

;------------------------------------------------------------------------------
; merge it onto a uniform grid
;------------------------------------------------------------------------------
case params.ndim of
    1: begin

        if (uniform_1d) then begin
            print, 'ERROR: merging 1-d data onto a uniform grid not supported.'
        endif else begin

            print, '1-d data will not be uniformly spaced.'

; allocate storage for the vector
            top_blocks = (size(where(tree[*].nodeType EQ 1)))[1]

            if (double) then begin
                sData = dblarr(params.nxb*top_blocks)
                x = dblarr(params.nxb*top_blocks)
                xl = dblarr(params.nxb*top_blocks)
                xr = dblarr(params.nxb*top_blocks)
            endif else begin
                sData = fltarr(params.nxb*top_blocks)
                x = fltarr(params.nxb*top_blocks)
                xl = fltarr(params.nxb*top_blocks)
                xr = fltarr(params.nxb*top_blocks)
            endelse

            ipos = 0

            for block = 0, params.totBlocks-1 do begin

                if (tree[block].nodeType EQ 1) then begin
                    xmin = tree[block].bndBox[0,0]
                
                    x_block = (tree[block].size[0]/params.nxb)* $
                      (findgen(params.nxb) + 0.5) + xmin

                    x[ipos:ipos+params.nxb-1] = x_block

                    xl_block = (tree[block].size[0]/params.nxb)* $
                      (findgen(params.nxb)) + xmin

                    xl[ipos:ipos+params.nxb-1] = xl_block

                    xr_block = (tree[block].size[0]/params.nxb)* $
                      (findgen(params.nxb) + 1.0) + xmin

                    xr[ipos:ipos+params.nxb-1] = xr_block

                    sData[ipos:ipos+params.nxb-1] = reform(unk[*,block,*])
                    
                    ipos = ipos + params.nxb
                    
                endif
                
            endfor

; the data may not be in order if it was run on multiple processors,
; so sort it.
            sort_indices = sort(x)

            x = temporary(x[sort_indices])
            xl = temporary(xl[sort_indices])
            xr = temporary(xr[sort_indices])
            sData = temporary(sData[sort_indices])


        endelse

    end

    2: begin

        if (params.geometry EQ "CARTESIAN" OR $
            params.geometry EQ "CYLINDRICAL") then begin
        
            if (double) then begin
                help, unk
                sData = merge_amr(reform(unk, params.totBlocks, $
                                              params.nxb, $
                                              params.nyb, $
                                              params.nzb), $
                                  XMERGE=x, YMERGE=y, $
                                  XRANGE=xrange, YRANGE=yrange, $
                                  TREE=tree, PARAMETERS=params, $
                                  SAMPLE=sample, /DOUBLE)
            endif else begin
                print, tree

                sData = merge_amr(reform(unk, params.totBlocks, $
                                              params.nxb, $
                                              params.nyb, $
                                              params.nzb), $
                                  XMERGE=x, YMERGE=y, $
                                  XRANGE=xrange, YRANGE=yrange, $
                                  TREE=tree, PARAMETERS=params, $
                                  SAMPLE=sample)
            endelse
        endif else if (params.geometry EQ "SPHERICAL") then begin
                sData = merge_polar(reform(unk, params.totBlocks, $
                                           params.nxb, $
                                           params.nyb, $
                                           params.nzb), $
                                        XMERGE=x, YMERGE=y, $
                                        RRANGE=xrange, TRANGE=yrange, $
                                        TREE=tree, PARAMETERS=params, $
                                        SAMPLE=sample)
        endif

    end
    3: begin
        if (double) then begin
            sData = merge_amr(reform(unk, params.totBlocks, $
                                          params.nxb, $
                                          params.nyb, $
                                          params.nzb), $
                              XMERGE=x, YMERGE=y, ZMERGE=z, $
                              XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange, $
                              TREE=tree, PARAMETERS=params, $
                              SAMPLE=sample, /DOUBLE)
        endif else begin
            sData = merge_amr(reform(unk, params.totBlocks, $
                                          params.nxb, $
                                          params.nyb, $
                                          params.nzb), $
                              XMERGE=x, YMERGE=y, ZMERGE=z, $
                              XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange, $
                              TREE=tree, PARAMETERS=params, $
                              SAMPLE=sample)
        endelse
    end
endcase


return, sData
end


