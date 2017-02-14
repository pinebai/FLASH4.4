pro draw_blocks_polar, color, ORIENTATION=orientation, TREE=tree, PARAMETERS=params
;
; draw the block boundaries on a plot
;


if n_elements(color) EQ 0 then color = 0

; loop over the blocks and plot the good data

i = (lonarr(1))[0]

for i = 1l, params.totBlocks do begin
    
    if tree[i-1].nodetype EQ 1 then begin


        case orientation of
            0: begin

              rt_to_xy, tree[i-1].bndBox[0,0], tree[i-1].bndBox[0,1], X=x1, Y=y1
              rt_to_xy, tree[i-1].bndBox[1,0], tree[i-1].bndBox[0,1], X=x2, Y=y2

; do several intermediate points
              deltat = tree[i-1].bndBox[1,1] - tree[i-1].bndBox[0,1]
              nstep = 5
              dt = deltat/nstep

              t_sub = findgen(nstep+0.5)*dt + tree[i-1].bndBox[0,1]
              
              x_sub = fltarr(nstep)
              y_sub = fltarr(nstep)

              for j = 0, nstep - 1 do begin
                  rt_to_xy, tree[i-1].bndBox[1,0], t_sub[j], X=xtmp, Y=ytmp
                  x_sub[j] = xtmp
                  y_sub[j] = ytmp
              endfor
               
;              rt_to_xy, tree[i-1].bndBox[1,0], 0.5*(tree[i-1].bndBox[1,1] + tree[i-1].bndBox[0,1]), X=x2a, Y=y2a
              
              rt_to_xy, tree[i-1].bndBox[1,0], tree[i-1].bndBox[1,1], X=x3, Y=y3
              rt_to_xy, tree[i-1].bndBox[0,0], tree[i-1].bndBox[1,1], X=x4, Y=y4

; only draw the outer radial arc
                plots, [x1, x2, x_sub, x3, x4], $
                  [y1, y2, y_sub, y3, y4], $
                  color = color, thick = 1, noclip=0

;                plots, [x1, x2, x2a, x3, x4], $
;                  [y1, y2, y2a, y3, y4], $
;                  color = color, thick = 1, noclip=0

            end
        endcase

;        xyouts, tree[i-1].coord[1], tree[i-1].coord[0], $
;           string(i,format = '(i3)'), align = .5, color = ltblue, $
;           charthick = 1.5
                
    endif
endfor

end


