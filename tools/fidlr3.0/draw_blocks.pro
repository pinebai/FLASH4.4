function do_draw, i, tree, level
    return, tree[i-1].lrefine EQ level OR (tree[i-1].lrefine LT level AND tree[i-1].nodetype EQ 1)
end



pro draw_blocks, color, ORIENTATION=orientation, TREE=tree, PARAMETERS=params, LEVEL=level
;
; draw the block boundaries on a plot
;


if n_elements(color) EQ 0 then color = 0
if n_elements(level) EQ 0 then level = 100000

xmax = max(tree[*].bndBox[1,0])
ymax = max(tree[*].bndBox[1,1])

; loop over the blocks and plot the good data

i = (lonarr(1))[0]

for i = 1l, params.totBlocks do begin

    if do_draw(i, tree, level) then begin
        ;(tree[i-1].bndBox[0,0] LE xmax AND tree[i-1].bndBox[0,1] LE ymax) then begin

        case orientation of
            0: begin
                plots, [tree[i-1].bndBox[0,0], tree[i-1].bndBox[0,0], $
                        tree[i-1].bndBox[1,0], tree[i-1].bndBox[1,0], $
                        tree[i-1].bndBox[0,0]], $
                  [tree[i-1].bndBox[0,1], tree[i-1].bndBox[1,1], $
                   tree[i-1].bndBox[1,1], tree[i-1].bndBox[0,1], $
                   tree[i-1].bndBox[0,1]], $
                  color = color, thick = 1, noclip=0
            end
            1: begin
                plots, [tree[i-1].bndBox[0,1], tree[i-1].bndBox[0,1], $
                        tree[i-1].bndBox[1,1], tree[i-1].bndBox[1,1], $
                        tree[i-1].bndBox[0,1]], $
                  [tree[i-1].bndBox[0,0], tree[i-1].bndBox[1,0], $
                   tree[i-1].bndBox[1,0], tree[i-1].bndBox[0,0], $
                   tree[i-1].bndBox[0,0]], $
                  color = color, thick = 1, noclip=0
            end
            2: begin
                plots, [tree[i-1].bndBox[0,0], tree[i-1].bndBox[0,0], $
                        tree[i-1].bndBox[1,0], tree[i-1].bndBox[1,0], $
                        tree[i-1].bndBox[0,0]], $
                  [tree[i-1].bndBox[0,1], tree[i-1].bndBox[1,1], $
                   tree[i-1].bndBox[1,1], tree[i-1].bndBox[0,1], $
                   tree[i-1].bndBox[0,1]], $
                  color = color, thick = 1, noclip=0
            end
        endcase

;        xyouts, tree[i-1].coord[1], tree[i-1].coord[0], $
;           string(i,format = '(i3)'), align = .5, color = ltblue, $
;           charthick = 1.5
                
    endif
endfor

; This is code for drawing the morton curve through the blocks.  There
; is no widget to turn this off or on, but that's what it needs.  Just
; always do it so we don't lose the code.  The gui should turn on and
; off including parents or just leaf blocks in the curve-- this can
; be changed by altering the do_draw routine, above.  This code right
; now.  Also the gui could control whether numbers are printed along
; the curve
                                                   
; Some potential problems: 1) I don't know if I'm drawing things
; correctly in the orientation 1 case below.  I just swapped x and y
; coordinates, but I don't know if that's correct. 2) The numbers
; printed along this line simply overlap with proc numbers if they're
; both printed out.  


;drawn_count = 0
;first = 1
;morton_color = color('dkblue')
;for i = 1l, params.totBlocks do begin

;    if do_draw(i, tree, level) then begin
;        drawn_count = drawn_count + 1

;        case orientation of
;            0: begin
;                if (first eq 1) then begin
;                    first = 0
;                    plots, tree[i-1].coord[0], tree[i-1].coord[1], noclip=0, color=morton_color
;                endif else begin
;                    plots, tree[i-1].coord[0], tree[i-1].coord[1], /CONTINUE, noclip=0, color=morton_color
;                endelse 
;                xyouts, tree[i-1].coord[0], tree[i-1].coord[1], $
;                        string(drawn_count), $
;                        alignment=1., noclip=0
;            end
;            1: begin
;                if (first eq 1) then begin
;                    first = 0
;                    plots, tree[i-1].coord[1], tree[i-1].coord[0], noclip=0, color=morton_color
;                endif else begin
;                    plots, tree[i-1].coord[1], tree[i-1].coord[0], /CONTINUE, noclip=0, color=morton_color
;                endelse 
;                xyouts, tree[i-1].coord[1], tree[i-1].coord[0], $
;                        string(drawn_count), $
;                        alignment=1., noclip=0

;            end
;            2: begin
;                if (first eq 1) then begin
;                    first = 0
;                    plots, tree[i-1].coord[0], tree[i-1].coord[1], noclip=0, color=morton_color
;                endif else begin
;                    plots, tree[i-1].coord[0], tree[i-1].coord[1], /CONTINUE, noclip=0, color=morton_color
;                endelse 
;                xyouts, tree[i-1].coord[0], tree[i-1].coord[1], $
;                        string(drawn_count), $
;                        alignment=1., noclip=0
;            end
;        endcase

;    endif
;endfor

end

