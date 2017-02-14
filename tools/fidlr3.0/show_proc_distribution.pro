pro show_proc_distribution, LEVEL=level $
  , ORIENTATION=orientation $
  , TREE=tree $
  , PARAMETERS=params $
  , SHOW_NUMBERS=show_numbers

;
; color each block based on which processor it lived
;

if (n_elements(show_numbers) EQ 0) then begin
    show_numbers = 0
end

xmax = max(tree[*].bndBox[1,0])
ymax = max(tree[*].bndBox[1,1])

; loop over the blocks and plot the good data

i = (lonarr(1))[0]

colors = ['red', 'green', 'orange', 'blue', 'yellow', 'purple']
numberOfColors = (size(colors))[1] 
for i = 1l, params.totBlocks do begin
    
    if tree[i-1].lrefine EQ level OR (tree[i-1].lrefine LT level AND tree[i-1].nodetype EQ 1) then begin

        blockColor = tree[i-1].processorNumber+1

        case orientation of
            0: begin
                
                polyfill, [tree[i-1].bndBox[0,0], tree[i-1].bndBox[0,0], $
                           tree[i-1].bndBox[1,0], tree[i-1].bndBox[1,0], $
                           tree[i-1].bndBox[0,0]], $
                          [tree[i-1].bndBox[0,1], tree[i-1].bndBox[1,1], $
                           tree[i-1].bndBox[1,1], tree[i-1].bndBox[0,1], $
                           tree[i-1].bndBox[0,1]], $
                          color = blockColor, noclip=0
                
                
                if (show_numbers EQ 1) then begin
                    xyouts, 0.5 * (tree[i-1].bndBox[0,0] + tree[i-1].bndBox[1,0]), $
                            0.5 * (tree[i-1].bndBox[0,1] + tree[i-1].bndBox[1,1]), $
                            string(tree[i-1].processorNumber), $
                            alignment=1., noclip=0
                end

            end
            1: begin
                polyfill, [tree[i-1].bndBox[0,1], tree[i-1].bndBox[0,1], $
                           tree[i-1].bndBox[1,1], tree[i-1].bndBox[1,1], $
                           tree[i-1].bndBox[0,1]], $
                          [tree[i-1].bndBox[0,0], tree[i-1].bndBox[1,0], $
                           tree[i-1].bndBox[1,0], tree[i-1].bndBox[0,0], $
                           tree[i-1].bndBox[0,0]], $
                          color = blockColor, noclip=0
                if (show_numbers EQ 1) then begin
                    xyouts, 0.5 * (tree[i-1].bndBox[0,0] + tree[i-1].bndBox[1,0]), $
                            0.5 * (tree[i-1].bndBox[0,1] + tree[i-1].bndBox[1,1]), $
                            string(tree[i-1].processorNumber), noclip=0, alignment=1.
                end
            end
            2: begin
                polyfill, [tree[i-1].bndBox[0,0], tree[i-1].bndBox[0,0], $
                           tree[i-1].bndBox[1,0], tree[i-1].bndBox[1,0], $
                           tree[i-1].bndBox[0,0]], $
                          [tree[i-1].bndBox[0,1], tree[i-1].bndBox[1,1], $
                           tree[i-1].bndBox[1,1], tree[i-1].bndBox[0,1], $
                           tree[i-1].bndBox[0,1]], $
                          color = blockColor, noclip=0
                if (show_numbers EQ 1) then begin
                    xyouts, 0.5 * (tree[i-1].bndBox[0,0] + tree[i-1].bndBox[1,0]), $
                            0.5 * (tree[i-1].bndBox[0,1] + tree[i-1].bndBox[1,1]), $
                            string(tree[i-1].processorNumber), noclip=0, alignment=1.
                end
                
                
            end
        endcase

;        xyouts, tree[i-1].coord[1], tree[i-1].coord[0], $
;           string(i,format = '(i3)'), align = .5, color = ltblue, $
;           charthick = 1.5
                
    endif
endfor

end
