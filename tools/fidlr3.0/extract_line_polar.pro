; given a point, a direction, and the block structured data for a
; single variable, return the line of data through that point in the 
; specified direction.

function extract_line_polar, temp_arr, $
                             TREE = tree, PARAMETERS = params, $
                             DIRECTION = dir, $
                             POINT = point, $
                             COORDS = coords1d


help, temp_arr

; compute the maximum number of point we are going to return (assuming
; a completely uniform grid)
lrefine_max = max(tree[*].lrefine)

case dir of 
    0: max_points = params.ntopx*params.nxb*2^(lrefine_max - 1)
    1: max_points = params.ntopy*params.nyb*2^(lrefine_max - 1)
endcase

coords = fltarr(max_points)
data1d = fltarr(max_points)


; find the leaf blocks where we pass through the given point in the
; specified direction

case dir of

; x-direction
; 
;   if we are going along the x-direction, then we don't care about
;   the x coord, but we want to target zones the contain our y point
;
    0: begin
        index_good = where(tree[*].nodetype EQ 1 AND $
                           tree[*].bndBox[0,1] LE point[1] AND $
                           tree[*].bndBox[1,1] GT point[1])


; now loop over all these good blocks and, if a zone contains our
; line, store it's center and the data.  Since we are Cartesian, and 
; our slice is along the x direction, we only need to check the first 
; row of y coords, and we get params.nxb points all at once
        nstored = 0

        index_base = where(tree[*].lrefine EQ 1)
        min_size = min(tree[index_base].size[0])
        max_size = max(tree[index_base].size[0])
        
        if min_size NE max_size then begin
            log_space = 1
            print, 'logarithmically spaced grid found'
        endif else begin
            log_space = 0
        endelse

        for n = 0, (size(index_good))[1] - 1 do begin

            blk = index_good[n]
 
; compute the y coords of the zone left and right edge
            yl = (tree[blk].bndBox[1,1] - tree[blk].bndBox[0,1])* $
              ((findgen(params.nyb))/float(params.nyb)) + $
              tree[blk].bndBox[0,1]

            yr = (tree[blk].bndBox[1,1] - tree[blk].bndBox[0,1])* $
              ((findgen(params.nyb) + 1.0)/float(params.nyb)) + $
              tree[blk].bndBox[0,1]
           
            j = where(yl[*] LE point[1] AND yr[*] GT point[1])

; we know the x coords -- they are independent of the y position 
; the block.  Sometimes, these may be log spaced
            
            if (log_space) then begin
                dlogr = (alog10(double(tree[blk].bndBox[1,0])) - $
                         alog10(double(tree[blk].bndBox[0,0])))/ $
                  double(params.nxb)

                rlcoords = 10.^(alog10(double(tree[blk].bndBox[0,0])) + $
                                dindgen(params.nxb)*dlogr)
                rrcoords = 10.^(alog10(double(tree[blk].bndBox[0,0])) + $
                                (dindgen(params.nxb)+1.0)*dlogr)
                x = 0.5*(rlcoords + rrcoords)

            endif else begin
                x = (tree[blk].bndBox[1,0] - tree[blk].bndBox[0,0])* $
                  ((findgen(params.nxb) + 0.5)/float(params.nxb)) + $
                  tree[blk].bndBox[0,0]
            endelse

            coords[nstored:nstored+params.nxb-1] = x[*]
            data1d[nstored:nstored+params.nxb-1] = temp_arr[blk,*,j]

            nstored = nstored + params.nxb

        endfor
    end
    
; y-direction
; 
;   if we are going along the y-direction, then we don't care about
;   the y coord, but we want to target zones the contain our x point
;    
    1: begin


        index_good = where(tree[*].nodetype EQ 1 AND $
                           tree[*].bndBox[0,0] LE point[0] AND $
                           tree[*].bndBox[1,0] GT point[0])


        index_base = where(tree[*].lrefine EQ 1)
        min_size = min(tree[index_base].size[0])
        max_size = max(tree[index_base].size[0])
        
        if min_size NE max_size then begin
            log_space = 1
            print, 'logarithmically spaced grid found'
        endif else begin
            log_space = 0
        endelse

; now loop over all these good blocks and, if a zone contains our
; line, store it's center and the data.  Since we are Cartesian, and
; our slice is along the y direction, we only need to check the first
; row of x coords, and we get params.nyb points all at once
        nstored = 0

        for n = 0, (size(index_good))[1] - 1 do begin

            blk = index_good[n]

; compute the x coords of the zone left and right edge
            if (log_space) then begin
                dlogr = (alog10(double(tree[blk].bndBox[1,0])) - $
                         alog10(double(tree[blk].bndBox[0,0])))/ $
                  double(params.nxb)

                xl = 10.^(alog10(double(tree[blk].bndBox[0,0])) + $
                          dindgen(params.nxb)*dlogr)
                xr = 10.^(alog10(double(tree[blk].bndBox[0,0])) + $
                          (dindgen(params.nxb)+1.0)*dlogr)

            endif else begin
                xl = (tree[blk].bndBox[1,0] - tree[blk].bndBox[0,0])* $
                  ((findgen(params.nxb))/float(params.nxb)) + $
                  tree[blk].bndBox[0,0]
                
                xr = (tree[blk].bndBox[1,0] - tree[blk].bndBox[0,0])* $
                  ((findgen(params.nxb) + 1.0)/float(params.nxb)) + $
                  tree[blk].bndBox[0,0]
            endelse

            i = where(xl[*] LE point[0] AND xr[*] GT point[0])

; we know the y coords -- they are independent of the y position
; the block
            y = (tree[blk].bndBox[1,1] - tree[blk].bndBox[0,1])* $
              ((findgen(params.nyb) + 0.5)/float(params.nyb)) + $
              tree[blk].bndBox[0,1]

            coords[nstored:nstored+params.nyb-1] = y[*]
            data1d[nstored:nstored+params.nyb-1] = temp_arr[blk,i,*]

            nstored = nstored + params.nyb

        endfor


    end

endcase

coords1d = coords[0:nstored-1]
data1d =  data1d[0:nstored-1]

index_sorted = sort(coords1d)

coords1d = coords1d[index_sorted]
data1d = data1d[index_sorted]

return, data1d

end






