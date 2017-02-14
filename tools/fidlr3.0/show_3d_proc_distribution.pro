pro show_3d_proc_distribution, LEVEL=level $
  , ORIENTATION=orientation $
  , SLICE_DIR = islice_dir $
  , TREE=tree $
  , PARAMETERS=params $
  , XRANGE = xrange $
  , YRANGE = yrange $
  , ZRANGE = zrange $
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
            
            blockColor = tree[i].processorNumber+1 mod 256
            
            case orientation of
                0: begin

                    polyfill, [tree[i].bndBox[0,ix], tree[i].bndBox[0,ix], $
                               tree[i].bndBox[1,ix], tree[i].bndBox[1,ix], $
                               tree[i].bndBox[0,ix]], $
                              [tree[i].bndBox[0,iy], tree[i].bndBox[1,iy], $
                               tree[i].bndBox[1,iy], tree[i].bndBox[0,iy], $
                               tree[i].bndBox[0,iy]], $
                              color = blockColor, noclip=0
                    
                    if (show_numbers EQ 1) then begin
                        xyouts, 0.5 * (tree[i].bndBox[0,ix] + tree[i].bndBox[1,ix]), $
                                0.5 * (tree[i].bndBox[0,iy] + tree[i].bndBox[1,iy]), $
                                string(tree[i].processorNumber), $
                                alignment=1.,noclip=0
                    end
                            
                 end
                1: begin

                    polyfill, [tree[i].bndBox[0,iy], tree[i].bndBox[0,iy], $
                               tree[i].bndBox[1,iy], tree[i].bndBox[1,iy], $
                               tree[i].bndBox[0,iy]], $
                              [tree[i].bndBox[0,ix], tree[i].bndBox[1,ix], $
                               tree[i].bndBox[1,ix], tree[i].bndBox[0,ix], $
                               tree[i].bndBox[0,ix]], $
                              color = blockColor, noclip=0
                    if (show_numbers EQ 1) then begin
                        xyouts, 0.5 * (tree[i].bndBox[0,ix] + tree[i].bndBox[1,ix]), $
                                0.5 * (tree[i].bndBox[0,iy] + tree[i].bndBox[1,iy]), $
                                string(tree[i].processorNumber), noclip=0, alignment=1.
                    end
                    
                end
                2: begin

                    polyfill, [tree[i].bndBox[0,ix], tree[i].bndBox[0,ix], $
                               tree[i].bndBox[1,ix], tree[i].bndBox[1,ix], $
                               tree[i].bndBox[0,ix]], $
                              [tree[i].bndBox[0,iy], tree[i].bndBox[1,iy], $
                               tree[i].bndBox[1,iy], tree[i].bndBox[0,iy], $
                               tree[i].bndBox[0,iy]], $
                              color = blockColor, noclip=0
                    
                    if (show_numbers EQ 1) then begin
                        xyouts, 0.5 * (tree[i].bndBox[0,ix] + tree[i].bndBox[1,ix]), $
                                0.5 * (tree[i].bndBox[0,iy] + tree[i].bndBox[1,iy]), $
                                string(tree[i].processorNumber), noclip=0, alignment=1.
                    end
                    
                end
            endcase
        endif                
    endif
endfor

end
