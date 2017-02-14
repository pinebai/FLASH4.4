;
; returns a 1d slice through 2d data                      2003/02 robi banerjee
;
; TO DO: consider problem.orientation

function extract3d_line, data2d, $
                       POINT = point, $
                       XMERGE = x, $
                       YMERGE = y, $
                       ZMERGE = z, $
                       SLICE_DIR = islice_dir, $
                       DIRECTION = dir, $
                       PROBLEM_INFO = problem, $
                       COORDS = coords1d


   help, data2d

; default return values
   coords1d = 0.0
   data1d   = 0.0

   nx = n_elements(x)
   ny = n_elements(y)
   nz = n_elements(z)
;
; assume equidistant spacing
;
   dxh = abs((x[nx-1]-x[0])/float(nx-1))/2.
   dyh = abs((y[ny-1]-y[0])/float(ny-1))/2.
   dzh = abs((z[nz-1]-z[0])/float(nz-1))/2.

   case islice_dir of
       0: begin ; x-y plane
           case dir of ; along x direction
               0: begin
                   idx = where( y-dyh le point[1] AND y+dyh gt point[1], rc)
                   if rc gt 0 then begin 
                       print, 'draw x-line at (y,z) = ',y[idx[0]],z[0]
                       data1d = data2d(*,idx[0])
                       coords1d = x
                   endif
               end
               1: begin ; along y direction
                   idx = where( x-dxh le point[0] AND x+dxh gt point[0], rc)
                   if rc gt 0 then begin
                       print, 'draw y-line at (x,z) = ',x[idx[0]],z[0]
                       data1d = data2d(idx[0],*)
                       coords1d = y
                   endif
               end
           endcase
       end

       1: begin ; x-z plane
           case dir of ; along x direction
               0: begin
                   idx = where( z-dzh le point[1] AND z+dzh gt point[1], rc)
                   if rc gt 0 then begin
                       print, 'draw x-line at (y,z) = ',y[0],z[idx[0]]
                       data1d = data2d(*,idx[0])
                       coords1d = x
                   endif
               end
               1: begin ; along z direction
                   idx = where( x-dxh le point[0] AND x+dxh gt point[0], rc)
                   if rc gt 0 then begin
                       print, 'draw z-line at (x,y) = ',x[idx[0]],y[0]
                       data1d = data2d(idx[0],*)
                       coords1d = z
                   endif
               end
           endcase
       end

       2: begin ; y-z plane
           case dir of ; along y direction
               0: begin
                   idx = where( z-dzh le point[1] AND z+dzh gt point[1], rc)
                   if rc gt 0 then begin
                       print, 'draw y-line at (x,z) = ',x[0],z[idx[0]]
                       data1d = data2d(*,idx[0])
                       coords1d = y
                   endif
               end
               1: begin ; along z direction
                   idx = where( y-dyh le point[0] AND y+dyh gt point[0], rc)
                   if rc gt 0 then begin
                       print, 'draw z-line at (x,y) = ',x[0],y[idx[0]]
                       data1d = data2d(idx[0],*)
                       coords1d = z
                   endif
               end
           endcase
       end
   endcase

   help, data1d, coords1d

   return, data1d
end






