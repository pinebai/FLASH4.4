pro data_browse_event, ev

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

common save, current_block, current_var
common save_ptr, tree_ptr, data_ptr, params_ptr, xmerge_ptr, ymerge_ptr, zmerge_ptr

; take action pased on which widget was changed
case uval of

    'var': begin
        print, 'in here'
        ivar = widget_info(info.vars, /droplist_select)
        update_data, $
          WIDGET_INFORMATION=info, $
          NEW_BLOCK=current_block, $
          VARIABLE=ivar
                 
    end

    'jump': begin

; get the block number from the field
        widget_control, info.blockJump, get_value=temp
        new_block = temp[0]

        if (new_block GT params_ptr.totBlocks) then begin
            result = $
              dialog_message(['ERROR: requested block > total # of blocks'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
            widget_control, info.blockJump, set_value=current_block
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'find': begin
        widget_control, info.xblock, get_value=temp
        xcoord = temp[0]

        widget_control, info.yblock, get_value=temp
        ycoord = temp[0]

        index = where(tree_ptr[*].nodeType EQ 1 AND $
                      tree_ptr[*].bndBox[0,0] LT xcoord AND $
                      tree_ptr[*].bndBox[1,0] GE xcoord AND $
                      tree_ptr[*].bndBox[0,1] LT ycoord AND $
                      tree_ptr[*].bndBox[1,1] GE ycoord)

        new_block = index[0] 

        if (new_block EQ -1) then begin
            result = $
              dialog_message(['ERROR: requested coords are outside domain'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        
            widget_control, info.xblock, $
              set_value=tree_ptr[current_block-1].coord[0]

            widget_control, info.yblock, $
              set_value=tree_ptr[current_block-1].coord[1]

        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block+1, $
              VARIABLE=current_var
            
        endelse
        
    end

    'parent': begin
        new_block = tree_ptr[current_block-1].gid[4]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: parent block does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'xlneigh': begin
        new_block = tree_ptr[current_block-1].gid[0]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: -X neighbor does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'xrneigh': begin
        new_block = tree_ptr[current_block-1].gid[1]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: +X neighbor does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'ylneigh': begin
        new_block = tree_ptr[current_block-1].gid[2]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: -Y neighbor does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'yrneigh': begin
        new_block = tree_ptr[current_block-1].gid[3]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: +Y neighbor does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'child1': begin
        new_block = tree_ptr[current_block-1].gid[5]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: -X/-Y child does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'child2': begin
        new_block = tree_ptr[current_block-1].gid[6]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: +X/-Y child does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'child3': begin
        new_block = tree_ptr[current_block-1].gid[7]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: -X/+Y child does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end

    'child4': begin
        new_block = tree_ptr[current_block-1].gid[8]

        if (new_block LT 0) then begin
            result = $
              dialog_message(['ERROR: +X/+Y child does not exist'], $
                             title='ERROR', $
                             dialog_parent=info.mainBase)
        endif else begin
            update_data, $
              WIDGET_INFORMATION=info, $
              NEW_BLOCK=new_block, $
              VARIABLE=current_var
        endelse

    end


    'exit': begin
        widget_control, ev.top, /destroy
    end


endcase

end


pro update_data, WIDGET_INFORMATION=info, $
                 NEW_BLOCK=block, $
                 VARIABLE=ivar
                 
common save, current_block, current_var
common save_ptr, tree_ptr, data_ptr, params_ptr

widget_control, info.coord_x, set_value = tree_ptr[block-1].coord[0]
widget_control, info.coord_y, set_value = tree_ptr[block-1].coord[1]

widget_control, info.size_x, set_value = tree_ptr[block-1].size[0]
widget_control, info.size_y, set_value = tree_ptr[block-1].size[1]

widget_control, info.refine, set_value = tree_ptr[block-1].lrefine
widget_control, info.nodetype, set_value = tree_ptr[block-1].nodeType

widget_control, info.data, set_value = reverse(reform(data_ptr[ivar,block-1,*,*,*]),2)

widget_control, info.blockJump, set_value = block
widget_control, info.xblock, set_value = tree_ptr[block-1].coord[0]
widget_control, info.yblock, set_value = tree_ptr[block-1].coord[1]

if (tree_ptr[block-1].gid[4] LT 0) then begin
    widget_control, info.parent, sensitive=0
endif else begin
    widget_control, info.parent, sensitive=1
endelse

if (tree_ptr[block-1].gid[0] LT 0) then begin
    widget_control, info.xlneigh, sensitive=0
endif else begin
    widget_control, info.xlneigh, sensitive=1
endelse

if (tree_ptr[block-1].gid[1] LT 0) then begin
    widget_control, info.xrneigh, sensitive=0
endif else begin
    widget_control, info.xrneigh, sensitive=1
endelse

if (tree_ptr[block-1].gid[2] LT 0) then begin
    widget_control, info.ylneigh, sensitive=0
endif else begin
    widget_control, info.ylneigh, sensitive=1
endelse

if (tree_ptr[block-1].gid[3] LT 0) then begin
    widget_control, info.yrneigh, sensitive=0
endif else begin
    widget_control, info.yrneigh, sensitive=1
endelse

if (tree_ptr[block-1].gid[5] LT 0) then begin
    widget_control, info.child1, sensitive=0
endif else begin
    widget_control, info.child1, sensitive=1
endelse

if (tree_ptr[block-1].gid[6] LT 0) then begin
    widget_control, info.child2, sensitive=0
endif else begin
    widget_control, info.child2, sensitive=1
endelse

if (tree_ptr[block-1].gid[7] LT 0) then begin
    widget_control, info.child3, sensitive=0
endif else begin
    widget_control, info.child3, sensitive=1
endelse

if (tree_ptr[block-1].gid[8] LT 0) then begin
    widget_control, info.child4, sensitive=0
endif else begin
    widget_control, info.child4, sensitive=1
endelse

current_block = block
current_var = ivar

grid_visualize, BLOCK=current_block


end



pro grid_visualize, BLOCK=block

common save_ptr, tree_ptr, data_ptr, params_ptr

; we want to visualize the vicinity around the selected block.  We
; will do this by making sure that atleast 3 levels of refinement are
; shown, and that the current block is not on the boundary of our
; visualization, unless it is on the boundary of the computational
; domain

lmax = max(tree_ptr[*].lrefine)

lblock = tree_ptr[block-1].lrefine


labove = lmax - lblock

lmin = lmax - 4 < lblock-1

lmin = lmin > 1

; if the current block is not already at the minimum level of
; refinement, then find the parent of it that is
if (lblock NE lmin) then begin
    current_block = block

    for n = 1, lblock do begin
        parent = tree_ptr[current_block-1].gid[4]
        if tree_ptr[parent-1].lrefine EQ lmin then break
        current_block = parent
    endfor

    topmost_block = parent
endif else begin
    topmost_block = block
endelse


if (params_ptr.geometry EQ "CARTESIAN") then begin

; the range of our plot will be the bounding box of the topmost block
    plot_x_min = tree_ptr[topmost_block-1].bndBox[0,0]
    plot_x_max = tree_ptr[topmost_block-1].bndBox[1,0]
    
    plot_y_min = tree_ptr[topmost_block-1].bndBox[0,1]
    plot_y_max = tree_ptr[topmost_block-1].bndBox[1,1]

    plot, [plot_x_min,plot_x_max], [plot_y_min, plot_y_max], psym=2, $
      XRANGE=[plot_x_min, plot_x_max], YRANGE=[plot_y_min,plot_y_max], $
      XSTYLE=1, YSTYLE=1

; find all the blocks that fall within this range
    index_plot = where(tree_ptr[*].bndBox[0,0] GE plot_x_min AND $
                       tree_ptr[*].bndBox[1,0] LE plot_x_max AND $
                       tree_ptr[*].bndBox[0,1] GE plot_y_min AND $
                       tree_ptr[*].bndBox[1,1] LE plot_y_max)


    num_plot = (size(index_plot))[1]

    for i = 0, num_plot-1 do begin
        plots, [tree_ptr[index_plot[i]].bndBox[0,0], $
                tree_ptr[index_plot[i]].bndBox[0,0], $
                tree_ptr[index_plot[i]].bndBox[1,0], $
                tree_ptr[index_plot[i]].bndBox[1,0], $
                tree_ptr[index_plot[i]].bndBox[0,0]], $
               [tree_ptr[index_plot[i]].bndBox[0,1], $
                tree_ptr[index_plot[i]].bndBox[1,1], $
                tree_ptr[index_plot[i]].bndBox[1,1], $
                tree_ptr[index_plot[i]].bndBox[0,1], $
                tree_ptr[index_plot[i]].bndBox[0,1]], $
          thick = 1

    endfor


    plots, [tree_ptr[block-1].bndBox[0,0], $
            tree_ptr[block-1].bndBox[0,0], $
            tree_ptr[block-1].bndBox[1,0], $
            tree_ptr[block-1].bndBox[1,0], $
            tree_ptr[block-1].bndBox[0,0]], $
           [tree_ptr[block-1].bndBox[0,1], $
            tree_ptr[block-1].bndBox[1,1], $
            tree_ptr[block-1].bndBox[1,1], $
            tree_ptr[block-1].bndBox[0,1], $
            tree_ptr[block-1].bndBox[0,1]], $
      thick = 3
endif else if params_ptr.geometry EQ "SPHERICAL" then begin

    r_min = tree_ptr[topmost_block-1].bndBox[0,0]
    r_max = tree_ptr[topmost_block-1].bndBox[1,0]
    t_min = tree_ptr[topmost_block-1].bndBox[0,1]
    t_max = tree_ptr[topmost_block-1].bndBox[1,1]

; the range of our plot will be the bounding box of the topmost block
    rt_to_xy, tree_ptr[topmost_block-1].bndBox[0,0], $
              tree_ptr[topmost_block-1].bndBox[0,1], $
      X=x1, Y=y1
    rt_to_xy, tree_ptr[topmost_block-1].bndBox[1,0], $
              tree_ptr[topmost_block-1].bndBox[0,1], $
      X=x2, Y=y2
    rt_to_xy, tree_ptr[topmost_block-1].bndBox[0,0], $
              tree_ptr[topmost_block-1].bndBox[1,1], $
      X=x3, Y=y3
    rt_to_xy, tree_ptr[topmost_block-1].bndBox[1,0], $
              tree_ptr[topmost_block-1].bndBox[1,1], $
      X=x4, Y=y4

    plot_x_min = min([x1, x2, x3, x4])
    plot_x_max = max([x1, x2, x3, x4])

    plot_y_min = min([y1, y2, y3, y4])
    plot_y_max = max([y1, y2, y3, y4])

    dx_plot = plot_x_max - plot_x_min
    dy_plot = plot_y_max - plot_y_min

    if (dx_plot GT dy_plot) then begin
        plot_y_max = plot_y_min + dx_plot
    endif else if (dy_plot GT dx_plot) then begin
        plot_x_max = plot_x_min + dy_plot
    endif

    plot, [plot_x_min,plot_x_max], [plot_y_min, plot_y_max], psym=2, $
      XRANGE=[plot_x_min, plot_x_max], YRANGE=[plot_y_min,plot_y_max], $
      XSTYLE=1, YSTYLE=1
     
; find all the blocks that fall within this range
    index_plot = where(tree_ptr[*].bndBox[0,0] GE r_min AND $
                       tree_ptr[*].bndBox[1,0] LE r_max AND $
                       tree_ptr[*].bndBox[0,1] GE t_min AND $
                       tree_ptr[*].bndBox[1,1] LE t_max)

    num_plot = (size(index_plot))[1]

    for i = 0, num_plot-1 do begin

        rt_to_xy, tree_ptr[index_plot[i]].bndBox[0,0], $
                  tree_ptr[index_plot[i]].bndBox[0,1], $
          X=x1, Y=y1

        rt_to_xy, tree_ptr[index_plot[i]].bndBox[1,0], $
                  tree_ptr[index_plot[i]].bndBox[0,1], $
          X=x2, Y=y2

        deltat = tree_ptr[index_plot[i]].bndBox[1,1] - $
                 tree_ptr[index_plot[i]].bndBox[0,1]

        nstep = floor(10*2.^(lblock - tree_ptr[index_plot[i]].lrefine)) > 1

        dt = deltat/nstep

        t_sub = (findgen(nstep)+0.5)*dt + tree_ptr[index_plot[i]].bndBox[0,1]

        x_sub = fltarr(nstep)
        y_sub = fltarr(nstep)

        for j = 0, nstep-1 do begin
            rt_to_xy, tree_ptr[index_plot[i]].bndBox[1,0], t_sub[j], X=xtmp, Y=ytmp
            x_sub[j] = xtmp
            y_sub[j] = ytmp
        endfor

        rt_to_xy, tree_ptr[index_plot[i]].bndBox[1,0], $
                  tree_ptr[index_plot[i]].bndBox[1,1], $
          X=x3, Y=y3

        rt_to_xy, tree_ptr[index_plot[i]].bndBox[0,0], $
                  tree_ptr[index_plot[i]].bndBox[1,1], $
          X=x4, Y=y4

        if (index_plot[i] EQ block-1) then begin

            x_sub2 = fltarr(nstep)
            y_sub2 = fltarr(nstep)

            for j = nstep-1, 0, -1 do begin
                rt_to_xy, tree_ptr[index_plot[i]].bndBox[0,0], t_sub[j], $
                  X=xtmp, Y=ytmp
                x_sub2[nstep-j-1] = xtmp
                y_sub2[nstep-j-1] = ytmp
            endfor

            plots, [x1, x2, x_sub, x3, x4, x_sub2, x1], $
                   [y1, y2, y_sub, y3, y4, y_sub2, y1], thick=3
        endif else begin
            plots, [x1, x2, x_sub, x3, x4], $
                   [y1, y2, y_sub, y3, y4]
        endelse

    endfor
          
        
endif


; make the current block thickly outlined

end


pro data_browse, filename

; some variables will be saved
common save_ptr, tree_ptr, data_ptr, params_ptr
common save, current_block, current_var

; read in the file, and make sure it is a real FLASH file
itype = determine_file_type(filename)

if (itype EQ -1) then begin
    print, 'ERROR: invalid file'
    return
endif

ndim = determine_file_dimensionality(filename)

tree = {lrefine:0l, $
        nodeType:0l, $
        gid:lonarr(2*ndim+1+2^ndim), $
        coord:fltarr(ndim), $
        size:fltarr(ndim), $
        bndBox:fltarr(2,ndim)}
tree_ptr = ptr_new(tree)

params = {totBlocks:0, $
          corners:0, $
          ndim:0, $
          nvar:0, $
          nxb:0, $
          nyb:0, $
          nzb:0, $
          ntopx:1, $
          ntopy:1, $
          ntopz:1, $
          time:0.0, $
          dt:0.0}
params_ptr = ptr_new(params)

data_ptr = ptr_new(0.0)


; read in the data
read_amr, filename, $
  PARAMETERS=params_ptr, TREE=tree_ptr, DATA=data_ptr, STORED_VARS=varnames

current_block = 1
current_var = 0

; setup the widget
mainBase = widget_base(/column, title = 'FLASH data browser')

blockBase = widget_base(mainBase, /row)
dataBase = widget_base(mainBase, /column)

coordsBase = widget_base(blockBase, /column, /frame)
sizeBase = widget_base(blockBase, /column, /frame)

refineBase = widget_base(blockBase, /column, /frame)

varSelectBase = widget_base(blockBase, /row, /frame)

tableBase = widget_base(dataBase, /column)

etcBase = widget_base(dataBase, /row)
gidBase = widget_base(etcBase, /row, /frame)
zoomABase = widget_base(etcBase, /row, /frame)
zoomBBase = widget_base(etcBase, /column, /frame)
zoomCBase = widget_base(etcBase, /column)

parentBase = widget_base(gidBase, /column)
neighBase = widget_base(gidBase, /column)
childBase = widget_base(gidBase, /column)
childBaseA = widget_base(childBase, /row)
childBaseB = widget_base(childBase, /row)


coord_x = cw_field(coordsBase, title = "block coords, x:", $
                   value = tree_ptr[current_block-1].coord[0], /NOEDIT)

coord_y = cw_field(coordsBase, title = "              y:", $
                   value = tree_ptr[current_block-1].coord[1], /NOEDIT)

size_x = cw_field(sizeBase, title = "block size, x:", $
                   value = tree_ptr[current_block-1].size[0], /NOEDIT)

size_y = cw_field(sizeBase, title = "            y:", $
                   value = tree_ptr[current_block-1].size[1], /NOEDIT)

refine = cw_field(refineBase, title = "refinement level:", $
                   value = tree_ptr[current_block-1].lrefine, /NOEDIT)

nodetype = cw_field(refineBase, title = "nodetype:        ", $
                    value = tree_ptr[current_block-1].nodeType, /NOEDIT)

vars = widget_droplist(varSelectBase, title = 'Variables: ', uvalue = 'var', $
                       value = varnames)

data = widget_table(tableBase, $
                    ALIGNMENT=0, $
                    COLUMN_LABELS=string(indgen(params_ptr.nxb)), $
                    FORMAT='(g13.6)', $
                    COLUMN_WIDTH=105, $
                    X_SCROLL_SIZE=8, $
                    Y_SCROLL_SIZE=8, $
                    ROW_LABELS=string(reverse(indgen(params_ptr.nyb),1)), $
                    VALUE=reverse(reform(data_ptr[current_var,current_block-1,*,*]),2), $
                    /RESIZEABLE_COLUMNS)


parent = widget_button(parentBase, value = "parent", uvalue = "parent")
if (tree_ptr[current_block-1].gid[4] LT 0) then $
  widget_control, parent, sensitive=0

xlneigh = widget_button(neighBase, value = "-X neighbor", uvalue = "xlneigh")
xrneigh = widget_button(neighBase, value = "+X neighbor", uvalue = "xrneigh")
ylneigh = widget_button(neighBase, value = "-Y neighbor", uvalue = "ylneigh")
yrneigh = widget_button(neighBase, value = "+Y neighbor", uvalue = "yrneigh")

if (tree_ptr[current_block-1].gid[0] LT 0) then $
  widget_control, xlneigh, sensitive=0
if (tree_ptr[current_block-1].gid[1] LT 0) then $
  widget_control, xrneigh, sensitive=0
if (tree_ptr[current_block-1].gid[2] LT 0) then $
  widget_control, ylneigh, sensitive=0
if (tree_ptr[current_block-1].gid[3] LT 0) then $
  widget_control, yrneigh, sensitive=0

child1 = widget_button(childBaseB, value = "-X/-Y child", uvalue = "child1")
child2 = widget_button(childBaseB, value = "+X/-Y child", uvalue = "child2")
child3 = widget_button(childBaseA, value = "-X/+Y child", uvalue = "child3")
child4 = widget_button(childBaseA, value = "+X/+Y child", uvalue = "child4")

if (tree_ptr[current_block-1].gid[5] LT 0) then $
  widget_control, child1, sensitive=0
if (tree_ptr[current_block-1].gid[6] LT 0) then $
  widget_control, child2, sensitive=0
if (tree_ptr[current_block-1].gid[7] LT 0) then $
  widget_control, child3, sensitive=0
if (tree_ptr[current_block-1].gid[8] LT 0) then $
  widget_control, child4, sensitive=0

window, XSIZE=600, YSIZE=600
grid_visualize, BLOCK=current_block

jump = widget_button(zoomABase, value = "Jump to block: ", uvalue = "jump")
blockJump = cw_field(zoomABase, title = "", value=current_block)

xblock = cw_field(zoomBBase, title = "find block, x:", $
                  value=tree_ptr[current_block-1].coord[0])
yblock = cw_field(zoomBBase, title = "            y:", $
                  value=tree_ptr[current_block-1].coord[1])
find = widget_button(zoomBBase, value = "Find", uvalue = "find", xsize=50)

exit = widget_button(zoomCBase, value = " Exit ", uvalue = "exit")

widget_control, mainBase, /realize


info = {mainBase:mainBase, $
        blockBase:blockBase, $
        dataBase:dataBase, $
        coordsBase:coordsBase, $
        sizeBase:sizeBase, $
        refineBase:refineBase, $
        varSelectBase:varSelectBase, $
        coord_x:coord_x, $
        coord_y:coord_y, $
        size_x:size_x, $
        size_y:size_y, $
        refine:refine, $
        nodetype:nodetype, $
        vars:vars, $
        data:data, $
        parent:parent, $
        xlneigh:xlneigh, $
        xrneigh:xrneigh, $
        ylneigh:ylneigh, $
        yrneigh:yrneigh, $
        child1:child1, $
        child2:child2, $
        child3:child3, $
        child4:child4, $
        blockJump:blockJump, $
        jump:jump, $
        xblock:xblock, $
        yblock:yblock, $
        find:find, $
        exit:exit}

widget_control, mainBase, set_uvalue=info, /no_copy
xmanager, 'data_browse', mainBase


end
