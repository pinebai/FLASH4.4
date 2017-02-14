pro draw3d_blocks , color $
                , ORIENTATION=orientation $
                , SLICE_DIR = islice_dir $
                , TREE=tree $
                , PARAMETERS=params $
                , XRANGE = xrange $
                , YRANGE = yrange $
                , ZRANGE = zrange $
                , LEVEL = level

;
; draw the block boundaries on a plot
;


if n_elements(color) EQ 0 then color = 0
if n_elements(level) EQ 0 then level = 10000
xmax = max(tree[*].bndBox[1,0])
ymax = max(tree[*].bndBox[1,1])

; loop over the blocks and plot the good data

i = (lonarr(1))[0]

for i = 0l, params.totBlocks-1 do begin
    
    if tree[i].lrefine EQ level OR (tree[i].lrefine LT level AND tree[i].nodetype EQ 1) then begin

        case islice_dir of
            
            0: begin            ; x-y plane
                
                ix = 0
                iy = 1
                iz = 2
                zz = zrange[0]
                
            end

            1: begin            ; x-z plane
                
                ix = 0
                iy = 2
                iz = 1
                zz = yrange[0]
                
            end
            
            2: begin            ; y-z plane
                
                ix = 1
                iy = 2
                iz = 0
                zz = xrange[0]
                
            end
            
        endcase                 ; islice_dir

        if tree[i].bndBox[0,iz] le zz AND $
          zz le tree[i].bndBox[1,iz] then begin

            case orientation of
                0: begin
                    plots, [tree[i].bndBox[0,ix], tree[i].bndBox[0,ix], $
                            tree[i].bndBox[1,ix], tree[i].bndBox[1,ix], $
                            tree[i].bndBox[0,ix]], $
                      [tree[i].bndBox[0,iy], tree[i].bndBox[1,iy], $
                       tree[i].bndBox[1,iy], tree[i].bndBox[0,iy], $
                       tree[i].bndBox[0,iy]], $
                      color = color, thick = 1, noclip=0
                end
                1: begin
                    plots, [tree[i].bndBox[0,iy], tree[i].bndBox[0,iy], $
                            tree[i].bndBox[1,iy], tree[i].bndBox[1,iy], $
                            tree[i].bndBox[0,iy]], $
                      [tree[i].bndBox[0,ix], tree[i].bndBox[1,ix], $
                       tree[i].bndBox[1,ix], tree[i].bndBox[0,ix], $
                       tree[i].bndBox[0,ix]], $
                      color = color, thick = 1, noclip=0
                end
                2: begin
                    plots, [tree[i].bndBox[0,ix], tree[i].bndBox[0,ix], $
                            tree[i].bndBox[1,ix], tree[i].bndBox[1,ix], $
                            tree[i].bndBox[0,ix]], $
                      [tree[i].bndBox[0,iy], tree[i].bndBox[1,iy], $
                       tree[i].bndBox[1,iy], tree[i].bndBox[0,iy], $
                       tree[i].bndBox[0,iy]], $
                      color = color, thick = 1, noclip=0
                end

            endcase             ; orientation

        endif                   ; plane found
        
;        xyouts, tree[i].coord[1], tree[i].coord[0], $
;           string(i,format = '(i3)'), align = .5, color = ltblue, $
;           charthick = 1.5
                
    endif
endfor

end
