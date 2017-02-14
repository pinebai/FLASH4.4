;-----------------------------------------------------------------------------
; xplot1d_amr.pro
;
;   plot output from flash -- for 1d data
;
;   this routine is typically driven by xflash
;
;   general flow:
;       -- read in data 
;       -- convert data from block structure to 2d array
;       -- scale data to colormap and display
;       -- plot block boundaries if desired
;
;-----------------------------------------------------------------------------

pro xplot1d_amr, FILE_INFO = fileInfo, $
                     VARIABLE_INFO = variable, $
                     OPTIONS = options, $
                     OUTPUT = output, $
                     ZOOM = zoom, $
                     COLORMAP = colormap, $
                     WIDGET_INFORMATION = info, $
                     PARTICLE_WIDGET = particle_widget, $
                     LABEL_OPT = label_info, $
                     KNOWN_VARIABLES = varnames, $
                     PROBLEM_INFO = problem, $
                     TREE_PTR = tree_ptr, $
                     PARAMS_PTR = params_ptr, $
                     DATA_PTR = data_ptr, $
                     PLOT_COUNT = plot_count, $
                     DEBUG = debug

common save, filename, pblm_ctr_lim, orientation, ndim, geometry

IF (Keyword_Set(DEBUG)) THEN print, '*** in xplot1d_amr'

xflash_dir = get_xflash_path()
IF Keyword_Set(debug) THEN print, 'looking for ', xflash_dir + 'flash_colors.tbl'
loadct, 3, file = xflash_dir + 'flash_colors.tbl', /SILENT

; eventually, the color to use should be picked from the pull down
; menu on the main widget -- this needs to be changed depending 
; on the resolution

use_color = color('red')

; set the plot title into titleName
set_plot_title,TITLE=titleName,VARIABLE_INFO=variable,OPTIONS=options



set_plot, 'X'

!p.background = color('white')
!p.color      = color('black')

;set the number of plots on a page
set_number_plots,PLOT_COUNT=plot_count,NX=nx,NY=ny,MULTI_PLOTS=multi_plots_per_page



!p.multi      = [0,NX,NY]

if (multi_plots_per_page EQ 1) then begin
    !x.omargin = [5,5]
    !y.omargin = [10,10]
    !x.margin = [10,3]
    !y.margin = [5,5]
endif else begin
    !x.omargin = [5,5]
    !y.omargin = [10,3]
    !x.margin = [10,3]
    !y.margin = [4,2]
endelse


IF Keyword_Set(debug) THEN print, 'multi_plots_per_page = ', multi_plots_per_page

IF  (output.type EQ 1) then begin
    !p.charsize   = 1.0
endif else BEGIN 
    if output.hsize GT 1200 then begin
        !p.charsize  = 3.00
        !p.charthick = 2.00
    endif else if output.hsize GT 600 then begin
        !p.charsize  = 1.25
        !p.charthick = 1.00
    endif else begin
        !p.charsize  = 0.90
        !p.charthick = 1.00
    endelse 
ENDELSE 


; -1 is default system
; 0 is 
; 1 is truetype fonts
; !p.font=-1


; FLASH files end with a 4-digit file number that is incremented as
; files are produced.  To get the base, strip the end off of the
; filename.  For now, ignore the possiblity of a .gz ending.
fileBase = strmid(fileInfo.prototype,0,strlen(fileInfo.prototype)-4)

pathEnd = strpos(fileBase, '/', /REVERSE_SEARCH)
fileBaseClean = strcompress(strmid(fileBase, pathEnd+1, $
                                   strlen(fileBase)-pathEnd))


;==============================================================================
; loop over the files
;==============================================================================
for ifile = fileInfo.startSuffix, fileInfo.endSuffix, fileInfo.step do begin


; set the current plot # and the x and y position (in plotspace)
; of the current plot on the page -- 1 based indexing here
    plot_number = !p.multi[0]
    if (plot_number EQ 0) then plot_number = NX*NY
    plot_number = NX*NY - plot_number 
    
    xplot_pos = ( (plot_number) mod NX ) + 1
    yplot_pos = (plot_number - xplot_pos + 1)/NX + 1
    
    IF Keyword_Set(debug) THEN print, '## plot info:', plot_number, xplot_pos, yplot_pos

; ---- determine the filename -------------------------------------------------
    filename = fileBase + string(ifile, format = '(i4.4)')

    IF Keyword_Set(debug) THEN print, 'working on ' + filename

; later on, we'll figure out whether the file has been read already,
; or if we need to unzip it.  For now read the file
    itype = determine_file_type(filename)
    
    if (itype EQ -1) then begin
        result = dialog_message('Error: file does not exist', $
                                /error, title = 'ERROR!!')
        return
    endif
    
    if (n_elements(info) NE 0) then widget_control, info.status, $
      set_value = 'status: reading from file'

; --- read the data -----------------------------------------------------------
    read_amr, filename, TREE=tree, DATA=unk, PARAMETERS=params, PARTICLES=particles, $
              GEOMETRY=geometry
    
    tree_ptr = tree
    data_ptr = unk
    params_ptr = params

    if (n_elements(info) NE 0) then $
      widget_control, info.status, $
                      set_value = 'status: done reading'
    
    
; dump a few statistics on the screen
    IF Keyword_Set(debug) THEN print, '   number of blocks:           ', params.totBlocks
    

    if params.time LT 1.e-21 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e24) + ' ys'
    endif else if params.time LT 1.e-18 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e21) + ' zs'
    endif else if params.time LT 1.e-15 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e18) + ' as'
    endif else if params.time LT 1.e-12 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e15) + ' fs'
    endif else if params.time LT 1.e-9 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e12) + ' ps'
    endif else if params.time LT 1.e-6 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e9) + ' ns'
    endif else if params.time LT 1.e-3 then begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time*1e6) + ' !4l!3s'

    endif else if params.time GT 1.0e15 then begin
        timeout  = 'time = ' + $
          string(format = '(f7.3)',params.time/3.1557e16) + ' Gyr'

    endif else if params.time GT 1.0e12 then begin
        timeout  = 'time = ' + $
          string(format = '(f7.3)',params.time/3.1557e13) + ' Myr'

    endif else if params.time GT 99999.9 then begin
        timeout  = 'time =' + $
          string(format = '(f15.1)',params.time) + ' s'

    endif else if params.time GT 999.9 then begin
        timeout  = 'time = ' + $
          string(format = '(f7.1)',params.time) + ' s'

    endif else begin
        timeout  = 'time = ' + $
                   string(format = '(f7.3)',params.time) + ' s'
    endelse



; ---- create the output filename ---------------------------------------------
    
    outfile = fileBaseClean
    
    if options.blocks then outfile = outfile + 'blk_'
    
    outfile = outfile + strcompress(variable.name, /REMOVE_ALL) + $
              string(ifile, format = '(i4.4)')
    
    case output.type of
        1: outfile = outfile + '.ps'
        2: begin
            
; newer version of IDL no longer support GIFs
            if (!VERSION.RELEASE LE 5.3) then begin
                outfile = outfile + '.gif'
            endif else begin
                outfile = outfile + '.png'
            endelse
        end
        else:
    endcase
    
    
; create the x window if necessary
    if ifile EQ fileInfo.startSuffix then begin
        if NOT output.type then begin
            window, xsize = output.hsize, $
                    ysize = output.vsize, $
                    title = 'xflash1D plot -- ' + filename 
        endif
    endif
    

; temporary storage for data manipulation
    top_blocks = (size(where(tree[*].nodeType EQ 1)))[1]
    aray   = fltarr(params.nxb*top_blocks)
    bray   = fltarr(params.nxb*top_blocks)
    
    
; grid centered variables
    ipos = 0

    for block = 0, params.totBlocks-1 do begin
        
        if (tree[block].nodeType EQ 1) then begin
            
            xmin = tree[block].bndBox[0,0]

            x_block = (tree[block].size[0]/params.nxb)* $
                      (findgen(params.nxb) + 0.5) + xmin

            aray[ipos:ipos+params.nxb-1] = x_block


            if (varnames(var_index(variable.name)) EQ 'snd_spd') then begin

                bray[ipos:ipos+params.nxb-1] = $
                  reform(sqrt(unk[var_index('gamc'),block,*,0,0] * $
                              unk[var_index('pres'),block,*,0,0] / $
                              unk[var_index('dens'),block,*,0,0]))
                

            endif else if (varnames(var_index(variable.name)) EQ $
                           'mach') then begin

                bray[ipos:ipos+params.nxb-1] = $
                  reform(abs(unk[var_index('velx'),block,*,0,0]) / $
                         sqrt(unk[var_index('gamc'),block,*,0,0] * $
                              unk[var_index('pres'),block,*,0,0] / $
                              unk[var_index('dens'),block,*,0,0]))
                

            endif else if (varnames(var_index(variable.name)) EQ $
                           'int_ener') then begin

                bray[ipos:ipos+params.nxb-1] = $
                  reform(unk[var_index('ener'),block,*,0,0] - $
                         .5*unk[var_index('velx'),block,*,0,0]^2)
                

            endif else if (varnames(var_index(variable.name)) EQ $
                           'ekin/eint') then begin

                bray[ipos:ipos+params.nxb-1] = $
                  reform(.5*unk[var_index('velx'),block,*,0,0]^2 / $
                         (unk[var_index('ener'),block,*,0,0]  $
                          - .5*unk[var_index('velx'),block,*,0,0]^2))


            endif else if (varnames(var_index(variable.name)) EQ $
                           'tot_vel') then begin

                bray[ipos:ipos+params.nxb-1] = $
                  abs(unk[var_index('velx'),block,*,0,0])


            endif else begin

                bray[ipos:ipos+params.nxb-1] = $
                  reform(unk(var_index(variable.name),block,*,0,0))
                

            endelse
            
            ipos = ipos + params.nxb
            
; end of nodetype test
        endif
        
;end of loop over blocks
    endfor


; absolute value the y data if requested
    if options.abs EQ 1 then bray = abs(bray)


; sort the arrays so if the points are connected, it makes sense       
    iind = sort(aray)
    x    = aray(iind)
    y    = bray(iind)
    

; print the data to a file if requested
;    if options.afile then begin
    
;        openw,2,outfile + '_raw.dat'
;        ac = '  node       x' + '      ' + titleName
;        printf,2,ac
;        print,ac
;        print,params.nxb
;        print,top_blocks
;        print,params.nxb*top_blocks
;
;        for i=0,params.nxb*top_blocks-1 do printf,2,FORMAT='(i4,100E13.5)',i,x[i],y[i]
;        close,2

;    endif

    
; state the min and max of the variable being plotted
    IF Keyword_set(debug) THEN print, varnames[var_index(variable.name)] + ' extrema: ', min(bray),max(bray)

    if (variable.auto EQ 1) then begin
        variable.min = min(bray)
        variable.max = max(bray)
    endif


; set the minimum by making sure we are between the minimum of the
; domain, and the minimum coord of the block in the max x direction
    if (zoom.xmin GE 0.d0) then begin
        xmin_tot = min(tree[*].bndBox[0,0])
        xmin = (zoom.xmin > xmin_tot) < max(tree[*].bndBox[0,0])
    endif else begin
        xmin = min(tree[*].bndBox[0,0])
        xmin_tot = xmin
    endelse
    
    if (zoom.xmax GE 0.d0) then begin
        xmax_tot = max(tree[*].bndBox[1,0])        
        xmax = (zoom.xmax < xmax_tot) > xmin*1.00001
    endif else begin
        xmax = max(tree[*].bndBox[1,0])
        xmax_tot = xmax
    endelse

    if (variable.auto EQ 1) then begin

; the y axis limits set automatically

        ymin = min(bray)
        ymax = max(bray)

        dy_ = ymax - ymin

        if dy_ NE 0. then begin
            if options.log then begin
                lymin = alog10(ymin)
                lymax = alog10(ymax)
                ldy   = lymax - lymin
                ymin  = 10.^(lymin-0.05*ldy)
                ymax  = 10.^(lymax+0.05*ldy)
            endif else begin
                ymin = ymin - 0.05*dy_
                ymax = ymax + 0.05*dy_
            endelse
        endif

        if ymin EQ ymax then begin
            ymin = ymin - 0.05*dy_
            ymax = ymax + 0.05*dy_
        endif

        if ymin EQ ymax then begin
            ymin = 0.95*ymin
            ymax = 1.05*ymax
        endif

        if ymin EQ ymax then begin
            ymin = -0.1
            ymax = +0.1
        endif

        if options.log then begin
            lymin = alog10(ymin)
            lymax = alog10(ymax)
            ldy   = lymax - lymin
            ymin  = 10.^(lymin-0.05*ldy)
            ymax  = 10.^(lymax+0.05*ldy)
        endif

    endif else begin

; the y axis limits are specified by the user

        ymin = variable.min
        ymax = variable.max

    endelse
; set the axis names 
    xstring = 'Distance (cm)'
    ystring = titleName
  
    
; initialize postscript if need be
    if (n_elements(info) NE 0) then $
      widget_control, info.status, $
                      set_value = 'status: initialize device'
    
    
; postscript portrait initialization
    if output.type EQ 1 then begin

        current_device = !d.name
        set_plot, 'PS'

        if (plot_number EQ 0) then begin
            xsize = 7.5
            ysize = 7.5
            device, file = outfile, xsize=xsize, ysize=ysize, retain=2, $
                    xoff=0.5, yoff=1.75, /inch, /color, bits_per_pixel = 8
        endif
    endif
    


; set the flags if log values are to be used
    xlogon = 0
    ylogon = 0
;    if options.xlog then xlogon=1
    if options.log then ylogon=1
    
    
; adjust the line thickness for postscript plots
    xthick1 = 1.0
    xthick2 = 2.0
    if output.type EQ 1 then xthick2 = 6.0
    
    
; straight plot
;    if options.oplot EQ 0 then begin

; first draw the axes
    IF Keyword_Set(debug) THEN print, 'looking to make the axes'

    if (multi_plots_per_page EQ 0) then begin
        plot, [0], [0], $
              xtitle = xstring, ytitle = ystring, $
              xrange = [xmin,xmax], yrange = [ymin,ymax], $
              xlog = xlogon, ylog = ylogon, $
              xstyle = 1, ystyle = 1, $
              color = black
    endif else begin
        plot, [0], [0], $
              xtitle = xstring, ytitle = ystring, $
              xrange = [xmin,xmax], yrange = [ymin,ymax], $
              xlog = xlogon, ylog = ylogon, $
              xstyle = 1, ystyle = 1, $
              color = black
    endelse

    oplot, x, y, $
           linestyle = 0, thick = xthick2, color = use_color
    
; plot the blocks
    if options.blocks then begin 

        usersym, [0, 0], [-4, 4]

        for i = 0, params.totBlocks-1 do begin

            if (tree[i].nodeType EQ 1) then begin

; use the bndBox to get the x position of the block boundary, and
; interpolate to find the data value at that point
                y_interp = interpol(y,x,[tree[i].bndBox[0,0]])

                oplot, [tree[i].bndBox[0,0]], [y_interp], psym=8, color=color('ltblue'), thick=xthick2
                
            endif

        endfor

    endif
    
    
; floating label

    if (label_info.enabled) then begin

        new_ymin = float(ymin)
        new_ymax = float(ymax)

        ypos = label_info.posy*(new_ymax-new_ymin)+new_ymin

        if (multi_plots_per_page EQ 0) then begin
            xyouts, label_info.posx*(xmax-xmin)+xmin, ypos, $
                    label_info.label, charsize = label_info.size, charthick = label_info.thick, $
                    color = label_info.color
        endif else begin

            labels = strarr(NX*NY)
            count2 = 0
            for count=0, 8 do begin
                if (label_info.label[count] NE '') then begin
                    labels[count2] = label_info.label[count]
                    if labels[count2] EQ '-1' then labels[count2] = ''
                    count2 = count2+1
                endif
            endfor

            xyouts, label_info.posx*(xmax-xmin)+xmin, ypos, $
                    labels[plot_number], alignment = .5, charsize = label_info.size, $
                    charthick = label_info.thick, color = label_info.color
        endelse

    endif


    if (particle_widget.enabled) then begin

      partvelvec, particles, PARTICLE_WIDGET = particle_widget, /OVER, $
        COLOR=color('verygray'), COUNTER=counter, PARTICLE_TAG=particleTag, $
        OLD_PARTS=old_parts, DIR=xflash_dir, PART_TEMP = parts_temp
    endif



;  finally add some information to the bottom of the plot
    gridtxt = 'number of blocks = ' + string(format = '(i7)', params.totBlocks)
    gridtxt2 = 'AMR levels = ' + string(format = '(i5)', max(tree[*].lrefine))

    if (multi_plots_per_page EQ 0) then begin

        if options.annotate EQ 1 then begin
            if output.type EQ 1 then begin
                xyouts, .05, .05, timeout, color = color('black'), /normal   

                xyouts, .05, .02, gridtxt + ', ' + gridtxt2, $
                        color = color('black'), /normal

                xyouts, .05, -.05, filename, color = color('dkgray'), $
                        /normal, charsize = .8
                
            endif else begin
                xyouts, .05, .07, timeout,  color = color('black'), /normal   
                xyouts, .05, .04, gridtxt,  color = color('black'), /normal  
                xyouts, .05, .01, gridtxt2, color = color('black'), /normal  
                
            endelse
        endif

    endif
    
; empty the buffer
    empty

    if (plot_number EQ NX*NY-1 OR ifile EQ fileInfo.endSuffix) then begin

        case output.type of
            1: begin
                device, /close
                set_plot, current_device
            end
            2: color_bitmap, outfile
            else:
        endcase
    endif
    
    if (n_elements(info) NE 0) then $
      widget_control, info.status, $
                      set_value = 'status: awaiting orders . . .'
    
endfor

end






