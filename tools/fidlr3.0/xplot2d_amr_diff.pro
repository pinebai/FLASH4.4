;-----------------------------------------------------------------------------
; xplot2d_amr_diff.pro
;
;    This function will plot the diff of a variable in 2 checkpoint
;    files.
;
;    general flow:
;       -- read in 2 data files 
;       -- convert data from block structure to 2d array
;       -- subtract 2 data files 
;       -- scale data to colormap and display
;       -- plot block boundaries if desired
;       -- plot velocity vectors, if desired
;
;-----------------------------------------------------------------------------

pro xplot2d_amr_diff, FILENAME1 = filename1, $
                      FILENAME2 = filename2, $
                      READ = iread, $
                      VARIABLE_INFO1 = variable1, $
                      VARIABLE_INFO2 = variable2, $
                      OPTIONS = options, $
                      OUTPUT = output, $
                      VECTOR = vector, $
                      ZOOM = zoom, $
                      COLORMAP = colormap, $
                      WIDGET_INFORMATION = info, $
                      CONTOUR_OPT = contours, $
                      PARTICLE_WIDGET = particle_widget, $
                      PROBLEM_INFO = problem, $
                      KNOWN_VARIABLES1 = varnames1, $
                      TREE1_PTR = tree1_ptr, $
                      PARAMS1_PTR = params1_ptr, $
                      DATA1_PTR = data1_ptr, $
                      KNOWN_VARIABLES2 = varnames2, $
                      TREE2_PTR = tree2_ptr, $
                      PARAMS2_PTR = params2_ptr, $
                      DATA2_PTR = data2_ptr, $
                      PLOT_COUNT = plot_count, $
                      OUTFILE_NAME = outfilename, $
                      CREDIT = credit, $
                      DEBUG = debug, $
                      NOINTERP = nointerp

; common in the variables from the read-in routine


; create a common to store the max variable quantity.  During the loop
; over the files, it is stored here
common maxvar, max_var

common save_state, params1, tree1, unk1, params2, tree2, unk2

common variables, varnames, native_num

; load the colortable -- xflash uses special colortables, the first 12
; colors are reserved and constant across the different tables
xflash_dir = get_xflash_path()

; get the colortable
iclrmap = color_index(colormap, MIN_VALUE=colorMin, MAX_VALUE=colorMax)
loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT

; set the plot title into titleName
set_plot_title,TITLE=titleName,VARIABLE_INFO=variable1,OPTIONS=options

set_plot, 'X'

!p.background = color('white')
!p.color      = color('black')

; set the number of plots on a page
set_number_plots,PLOT_COUNT=plot_count,NX=nx,NY=ny,MULTI_PLOTS=multi_plots_per_page,DEBUG=debug

!p.multi      = [0,NX,NY]

if (output.type EQ 1) then begin
    if (multi_plots_per_page) then begin
        !p.charsize   = 0.75
    endif else begin
        !p.charsize   = 1.0
    endelse

endif else begin
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
endelse

; FLASH files end with a 4-digit file number that is incremented as
; files are produced.  To get the base, strip the end off of the
; filename.  For now, ignore the possiblity of a .gz ending.
;fileBase = strmid(file1Info.prototype,0,strlen(file1Info.prototype)-4)

;pathEnd = strpos(fileBase, '/', /REVERSE_SEARCH)
;fileBaseClean = strcompress(strmid(fileBase, pathEnd+1, $
;                                   strlen(fileBase)-pathEnd))


    
; diff hack: forcing plot to middle of plotspace; assuming only one plot
xplot_pos = 1
yplot_pos = 1
    

; ---- determine the filename -------------------------------------------------
;    filename = fileBase + string(ifile, format = '(i4.4)')

;filename1 = file1Info.prototype
;filename2 = file2Info.prototype
IF Keyword_set(debug) THEN print, 'working on ' + filename1 + ' and ' + filename2

; only read the data if it hasn't changed since the last time
if (iread EQ 1) then begin

; later on, we'll figure out whether the file has been read already,
; or if we need to unzip it.  For now read the file
    itype1 = determine_file_type(filename1)
    itype2 = determine_file_type(filename2)
        
    if (itype1 EQ -1 OR itype2 EQ -1) then begin
        result = dialog_message('Error: file does not exist', $
                                /error, title = 'ERROR!!')
        return
    endif
    

    if (n_elements(info) NE 0) then widget_control, info.status, $
          set_value = 'status: reading from file'

; --- read the data -----------------------------------------------------------
    read_amr, filename1, TREE=tree1, DATA=unk1, PARAMETERS=params1, STORED_VARS=varnames1
    read_amr, filename2, TREE=tree2, DATA=unk2, PARAMETERS=params2, STORED_VARS=varnames2
        
    varnames = varnames1
    tree1_ptr = tree1
    data1_ptr = unk1
    params1_ptr = params1

    tree2_ptr = tree2
    data2_ptr = unk2
    params2_ptr = params2

    if (n_elements(info) NE 0) then $
      widget_control, info.status, $
      set_value = 'status: done reading'
endif else begin
    IF (Keyword_set(debug)) THEN print, '^^^^ skipping the read'
endelse

; dump a few statistics on the screen
IF Keyword_set(debug) THEN print, '   number of blocks:           ', params1.totBlocks

lrefine_max = max(tree1[*].lrefine)
top_level = (size(where(tree1[*].lrefine EQ lrefine_max)))[1]
IF Keyword_set(debug) THEN print, '   number of top level zones: ', top_level*params1.nxb*params1.nyb

uniform = 0l
uniform = params1.ntopx*params1.nxb*2^(lrefine_max-1)* $
  params1.ntopy*params1.nyb*2^(lrefine_max-1)

IF Keyword_set(debug) THEN print, '   fraction of uniform grid:   ', $
  float(top_level*params1.nxb*params1.nyb)/float(uniform)

if params1.time LT 1.e-9 then begin
    timeout  = 'time = ' + $
      string(format = '(f9.3)',params1.time*1e12) + ' ps'

endif else if params1.time LT 1.e-6 then begin
    timeout  = 'time = ' + $
      string(format = '(f9.3)',params1.time*1e9) + ' ns'

endif else if params1.time LT 1.e-3 then begin
    timeout  = 'time = ' + $
      string(format = '(f9.3)',params1.time*1e6) + ' !4l!3s'

endif else if params1.time LT 3600. then begin
    timeout  = 'time = ' + $
      string(format = '(f9.3)',params1.time) + ' s'

endif else if params1.time LT 86400. then begin
    timeout  = 'time = ' + $
      string(format = '(f9.3)',params1.time/3600.) + ' hours'

endif else if params1.time GT 1.0e15 then begin
    timeout  = 'time = ' + $
      string(format = '(f7.3)',params.time/3.1557e16) + ' Gyr'

endif else if params1.time GT 1.0e12 then begin
    timeout  = 'time = ' + $
      string(format = '(f7.3)',params.time/3.1557e13) + ' Myr'

endif else begin
    timeout  = 'time = ' + $
      string(format = '(f9.3)',params1.time/86400.) + ' days'

endelse
    

; ---- create the output filename ---------------------------------------------

outfile = filename1 + '_'

if options.blocks then outfile = outfile + 'blk_'

if n_elements(outfilename) EQ 0 then begin

    outfile = outfile + strcompress(variable1.name, /REMOVE_ALL) ;+ $
      ;string(ifile, format = '(i4.4)')
   
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
        
endif else begin
    outfile = outfilename
endelse
  
if (output.type eq 1) then print, 'outputting to ', outfile

;------------------------------------------------------------------------------
; load the variable to be plotted into temporary storage
;------------------------------------------------------------------------------
IF Keyword_set(debug) THEN print, 'variable1 name = ', variable1.name
varnames = varnames1
temp_arr1 = create_variable(unk1, variable1.name,debug)
varnames = varnames2
temp_arr2 = create_variable(unk2, variable2.name,debug)

;if options.abs EQ 1 then temp_arr = abs(temp_arr)

IF Keyword_set(debug) THEN print, variable1.name + ' extrema: ', min(temp_arr1), max(temp_arr1)


;------------------------------------------------------------------------------
; find the size of the plotting area
;------------------------------------------------------------------------------

; by default, set the plot extrema to the domain limits.  If this was
; overridden by the zoom boxes, then make sure that the value given
; for the limit is within the computational domain.
    
xmin_tot = min(tree1[*].bndBox[0,0])
xmax_tot = max(tree1[*].bndBox[1,0])        

ymin_tot = min(tree1[*].bndBox[0,1])
ymax_tot = max(tree1[*].bndBox[1,1])

; right now, a value of -1.0 means use the default.  This should be
; changed at some point to allow any valid negative coordinate.

if (zoom.xmin NE -1.0d0) then begin
    xmin = (zoom.xmin > xmin_tot) < max(tree1[*].bndBox[0,0])
    xmin = (xmin)[0]
endif else begin
    xmin = xmin_tot
endelse

if (zoom.ymin NE -1.0d0) then begin
    ymin = (zoom.ymin > ymin_tot) < max(tree1[*].bndBox[0,1])
    ymin = (ymin)[0]
endif else begin
    ymin = ymin_tot
endelse


if (zoom.xmax GE 0.d0) then begin
    xmax = (zoom.xmax < xmax_tot) > xmin*1.00001
    xmax = (xmax)[0]
endif else begin
    xmax = xmax_tot
endelse

if (zoom.ymax GE 0.d0) then begin
    ymax = (zoom.ymax < ymax_tot) > ymin*1.00001
    ymax = (ymax)[0]
endif else begin
    ymax = ymax_tot
endelse

; compute the aspect ratio, and define the plot size
; determine the window size and 
dy = ymax - ymin
dx = xmax - xmin

case problem.orientation of
    0: dataAspectRatio = dy/dx
    1: dataAspectRatio = dx/dy
    2: dataAspectRatio = dy/dx
endcase

; compute the fraction of the total domain that is being plotted
xfrac = dx/(xmax_tot - xmin_tot)
yfrac = dy/(ymax_tot - ymin_tot)

; ---- initialize the device --------------------------------------------------

plot_number = 0 

; choose orientation based on dataAspectRatio
if (n_elements(info) NE 0) then $
  widget_control, info.status, set_value = 'status: initialize device'

if output.type EQ 1 then begin
    current_device = !d.name
    set_plot, 'PS'

    if dataAspectRatio GE 1. AND NX LE NY then begin
        iorient = 0
    endif else begin
        iorient = 1
    endelse
    
    case iorient of

; portrait orientation
        0: begin

            xsize = 7.5
            ysize = 9.5
            
            if (plot_number EQ 0) then begin
                device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
                  XOFF = 0.5, YOFF = 0.75, /INCH, /COLOR, $
                  BITS_PER_PIXEL = 8
            endif

        end

; landscape orientation
        1: begin

            xsize = 10.
            ysize = 6.

            if (plot_number EQ 0) then begin
                device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
                  XOFF = 1.25, YOFF = 10.5, /INCH, /COLOR, $
                  BITS_PER_PIXEL = 8, /LANDSCAPE
            endif
            
        end
    endcase

    deviceAspect = (ysize)/(xsize)


endif else begin
    deviceAspect = (float(output.vsize))/ $
      (float(output.hsize))
endelse 

; ---- determine the bounds for the plot --------------------------------------

; set the normal coordinates of the portion of the display/page you
; wish to use -- leave a little margin so the plots don't run to the
; edge of the page
page_nx1 = 0.05
page_nx2 = 0.95

page_ny1 = 0.05
page_ny2 = 0.90
    

; compute the bounding normal coordinates for this subplot
; the current subplot occupies the rectangle with opposite
; corners (nx1, ny1) to (nx2, ny2).  (note, 0,0 is in the 
; lower left corner in normal coords).

dpagex = page_nx2 - page_nx1
dpagey = page_ny2 - page_ny1

nx1 = (dpagex/float(NX)) * (xplot_pos - 1 + 0.15) + page_nx1
nx2 = (dpagex/float(NX)) * (xplot_pos - 0.05)  + page_nx1

ny1 = page_ny2 - (dpagey/float(NY)) * (yplot_pos - 0.15) 
ny2 = page_ny2 - (dpagey/float(NY)) * (yplot_pos - 1 + 0.05) 

; compute the aspect ratio of the domain and window (paper) together
aspect = dataAspectRatio/deviceAspect

; aspect compares the aspect ratio of the data to that of the device,
; aspect = (dy_data/dy_device)/(dx_data/dx_device).  If aspect > 1, 
; it means that the y-coordinate is the one that is going to set the
; overall scaling.  If aspect < 1, then the x-coordinate will set the
; scaling.  Consider each case separately below.

if (aspect GE 1.) then begin 

; the y size of the data sets the scaling of the plot

; set the initial values of the y min and max normal coordinates.  We
; leave some room for the axis labels
    py1 = ny1
    py2 = ny2

; compute the x size, using the aspect ratio
    dpy = py2 - py1
    dpx = dpy/aspect < 0.7* (dpagex/float(NX))


; recompute dpy, in case dpx was adjusted
    dpy = aspect*dpx

; compute the plot coordinates
    px1 = nx1
    px2 = px1 + dpx

    py_center = 0.5*(ny1 + ny2)
    py1 = py_center - .5*dpy
    py2 = py_center + .5*dpy

; set the plot and legend bounding box -- the legend will be vertical
    plot_pos   = [px1, py1, px2, py2]

; the legend goes in the remaining space in the x direction -- figure
; out how much room is left there
    dx_legend = (dpagex/float(NX)) * (xplot_pos) +page_nx1 - px2
    lwidth = dx_legend/4 < 0.25*(px2 - px1)
    lcenter = 0.5* ((dpagex/float(NX)) * (xplot_pos) + page_nx1 + px2)
    legend_pos = [lcenter-0.5*lwidth, py1, lcenter+0.5*lwidth, py2]
    
endif else begin 
        
; the x size of the data sets the scaling of the plot


; set the initial x min and max normal coordiantes
    px1 = nx1
    px2 = nx2

    dpx = px2 - px1
    dpy = aspect*dpx < 0.7* (dpagey/float(NY))

; recompute dpx, in case dpy was adjusted
    dpx = dpy/aspect

; recompute the plot coordinates
    px_center = 0.5*(nx1 + nx2)
    px1 = px_center - .5*dpx 
    px2 = px_center + .5*dpx

    py2 = ny2 
    py1 = py2 - dpy

; set the plot and legend bounding box -- the legend will be horizontal
    plot_pos   = [px1, py1, px2, py2]

; the legend goes in the remaining space in the lower y direction --
; figure out how much room is left there
    dy_legend = py1 - (page_ny2 - (dpagey/float(NY)) * (yplot_pos))
    lheight = dy_legend/4 < 0.25*(py2 - py1)
    lcenter = 0.5*(py1 + (page_ny2 - (dpagey/float(NY)) * (yplot_pos)))
    legend_pos = [px1, lcenter-0.5*lheight, px2, lcenter+0.5*lheight]
endelse
    

; create the plot window if necessary
;if (ifile EQ file1Info.startSuffix) then begin

if (output.type NE 1) then $
  window, 0, XSIZE = output.hsize, YSIZE = output.vsize, $
  TITLE = 'AMR plot'

IF Keyword_set(debug) THEN print, 'main plot window = ', !D.WINDOW
;endif


;------------------------------------------------------------------------------
; determine whether to sub-sample the data
;------------------------------------------------------------------------------

; assume 1024 pixels on screen and 2048 on paper (~ 300dpi

sample = 0

max_refine = max(tree1[*].lrefine)

; compute the number of pixels in a uniform grid of our sub-domain
max_pixels = params1.ntopx*params1.nxb*2.^(max_refine-1)*xfrac > $
  params1.ntopy*params1.nyb*2.^(max_refine-1)*yfrac

if output.type EQ 1 then begin
    ideal = 2048
     
    if max_pixels GT ideal then $
      sample = fix(alog(max_pixels/ideal) / alog(2.)) > 0

endif else begin
    ideal = 1024

    if max_pixels GT ideal then $
      sample = fix(alog(max_pixels/ideal) / alog(2.)) > 0

endelse

; make sure we have atleast one point in every block
;    if (2^sample GT (params.nxb < params.nyb)) then begin
;        sample = fix(alog(params.nxb < params.nyb)/alog(2.))
;    endif

IF Keyword_set(debug) THEN print, 'sub-sampling by ', sample, ' levels of refinement'


;------------------------------------------------------------------------------
; scale the data
;------------------------------------------------------------------------------
if (n_elements(info) NE 0) then $
  widget_control, info.status, SET_VALUE = 'status: scaling data'


time_start = systime(/seconds)

IF Keyword_set(debug) THEN print, max(temp_arr1)

; first put it on a uniform grid
IF Keyword_set(debug) THEN openw, 1, 'nsdata1'

nsdata1 = merge_amr(temp_arr1, TREE = tree1, PARAMETERS = params1, $
                   XMERGE = x, YMERGE = y, SAMPLE = sample, $
                   XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DOUBLE=1)

IF Keyword_set(debug) THEN openw, 2, 'nsdata2'
nsdata2 = merge_amr(temp_arr2, TREE = tree2, PARAMETERS = params2, $
                   XMERGE = x, YMERGE = y, SAMPLE = sample, $
                   XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DOUBLE=1)

IF Keyword_set(debug) THEN printf, 1, nsdata1
close, 1
IF Keyword_set(debug) THEN printf, 2, nsdata2
close, 2
nsdata = nsdata1 - nsdata2

IF Keyword_set(debug) THEN print, nsdata
; if we are automagically scaling, reset the data range
if variable1.auto EQ 1 then begin
    variable1.min = min(nsdata)
    variable1.max = max(nsdata)
endif


IF Keyword_set(debug) THEN print, 'variable1 plot range: ', variable1.min, variable1.max

if options.log EQ 1 then nsdata = alog10(temporary(nsdata))

; if we are max'ing the variable, treat it specially
if options.max EQ 0 then begin
    case options.log of
        0: nsdata = scale_color(temporary(nsdata), $
                                VARMAX = variable1.max, $
                                VARMIN = variable1.min, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)

        1: nsdata = scale_color(temporary(nsdata), $
                                VARMAX = variable1.max, $
                                VARMIN = variable1.min, /log, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)
    endcase

    time_end = systime(/seconds)

    IF Keyword_set(debug) THEN print, '**** time to merge and scale = ', time_end - time_start
endif else begin

       
; take the maximum
    if ifile EQ file1Info.startSuffix then begin
        max_var = nsdata
        IF Keyword_set(debug) THEN print, 'max = ', max(max_var), min(max_var)
    endif else begin
        max_var = nsdata > max_var
    endelse

; scale it
    case options.log of
        0: nsdata = scale_color(max_var, $
                                VARMAX = variable1.max, $
                                VARMIN = variable1.min, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)
            
        1: nsdata = scale_color(max_var, $
                                VARMAX = variable1.max, $
                                VARMIN = variable1.min, /log, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)
    endcase
endelse
       
if problem.orientation EQ 1 then nsdata = transpose(temporary(nsdata))


;-----------------------------------------------------------------------------
; create the axes -- just get the data coordinates defined
;-----------------------------------------------------------------------------

; create the title -- if we are printing it
if (multi_plots_per_page EQ 0) then begin
    
    if (options.annotate EQ 1) then begin
        title = titleName
    endif else begin
        title = titleName ;' '
    endelse

endif else begin
        
    title = timeout
        
endelse

case problem.orientation of
    0: begin
        plot, [xmin, xmax], [ymin, ymax], POS = plot_pos, $
          XSTYLE = 5, YSTYLE = 5, TITLE = title
        oplot, [xmin, xmax], [ymax, ymin]
    end

    1: begin
; transposed orientation (x on the vertical axis)
        plot, [ymin, ymax], [xmin, xmax], POS = plot_pos, $
          XSTYLE = 5, YSTYLE = 5, TITLE = title
        oplot, [ymin, ymax], [xmax, xmin]
    end
    
    2: begin
        plot, [xmin, xmax], [ymin, ymax], POS = plot_pos, $
          XSTYLE = 5, YSTYLE = 5, YRANGE = [ymax,ymin], TITLE = title
        oplot, [xmin, xmax], [ymax, ymin]
    end
endcase

; if there are multiple plots on a page -- throw up the title
;if (multi_plots_per_page EQ 1 AND plot_number EQ 0) then begin
;    xyouts, .5, 0.5*(1.0 + page_ny2), titleName, $
;      color = color('black'), /normal, alignment = 0.5, charsize=2
;endif


;------------------------------------------------------------------------------
; plot the data
;------------------------------------------------------------------------------
if (n_elements(info) NE 0) then $
  widget_control, info.status, $
  set_value = 'status: plotting data'

; convert the limits of the domain into normal coordinates
case problem.orientation of
    0: begin
        lower = convert_coord([xmin, ymin], /DATA, /TO_NORMAL)
        upper = convert_coord([xmax, ymax], /DATA, /TO_NORMAL)
    end
    1: begin
        lower = convert_coord([ymin, xmin], /DATA, /TO_NORMAL)
        upper = convert_coord([ymax, xmax], /DATA, /TO_NORMAL)
    end
    2: begin
        lower = convert_coord([xmin, ymax], /DATA, /TO_NORMAL)
        upper = convert_coord([xmax, ymin], /DATA, /TO_NORMAL)
        nsdata = rotate(temporary(nsdata),7)
    end
endcase

sz = size(nsdata)

tvimage, nsdata, position = [lower[0], lower[1], upper[0], upper[1]], $
  /OVERPLOT, /HALF_HALF, NOINTERPOLATION=nointerp

;------------------------------------------------------------------------------
; plot any user-defined contours
;------------------------------------------------------------------------------

IF Keyword_set(debug) THEN print, 'contour 0: ', contours[0]

for i = 0, (size(contours))[1]-1 do begin
        
    if (contours[i].enabled) then begin
        
        if (n_elements(info) NE 0) then $
          widget_control, info.status, $
          set_value = 'status: plotting contour #' + string(i,'(i2)')

        sctr = merge_amr(reform(unk(contours[i].var,*,*,*,*)), $
                         TREE = tree, PARAMETERS = params, $
                         XMERGE = x, YMERGE = y, sample = sample, $
                         XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])
            
        IF Keyword_set(debug) THEN print, 'sctr: min, max: ', min(sctr), max(sctr)
        if problem.orientation EQ 1 then begin
            sctr = transpose(temporary(sctr))
            ty = y
            y = x
            x = ty
            undefine, ty
        endif

        contour, sctr, x, y, levels = [contours[i].value], $
          c_color = [contours[i].color], /OVERPLOT
        
        undefine, sctr
    endif
        
endfor

;------------------------------------------------------------------------------
; plot the vectors
;------------------------------------------------------------------------------

if vector.enabled AND $
  var_index(vector.xcomp) GE 0 AND $
  var_index(vector.ycomp) GE 0 then begin

  velx = drop_1st_dimension(unk[var_index(vector.xcomp),*,*,*,*])
  vely = drop_1st_dimension(unk[var_index(vector.ycomp),*,*,*,*])
    
  if (n_elements(info) NE 0) then $
    widget_control, info.status, $
    set_value = 'status: plotting vectors'
        
  case problem.orientation of
    0: begin
      velx = merge_amr(velx, TREE = tree, PARAMETERS = params, $
                       sample = sample, $
                       XMERGE = x, YMERGE = y, $
                       XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])
            
      vely = merge_amr(vely, TREE = tree, PARAMETERS = params, $
                       sample = sample, $
                       XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])

    end
    1: begin
      vely = merge_amr(velx, TREE = tree, PARAMETERS = params, $
                       sample = sample, $
                       XMERGE = x, YMERGE = y, $
                       XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])

      vely = transpose(temporary(vely))

      velx = merge_amr(vely, TREE = tree, PARAMETERS = params, $
                       sample = sample, $
                       XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])

      velx = transpose(temporary(velx))
    end
    2: begin
      velx = merge_amr(velx, TREE = tree, PARAMETERS = params, $
                       sample = sample, $
                       XMERGE = x, YMERGE = y, $
                       XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])

      vely = merge_amr(vely, TREE = tree, PARAMETERS = params, $
                       sample = sample, $
                       XRANGE = [xmin,xmax], YRANGE = [ymin,ymax])

    end
  endcase

  if (problem.orientation EQ 1) then begin
    ty = y
    y = x
    x = ty
    undefine, ty
  endif
        
  xpts = (size(x))[1]
  ypts = (size(y))[1]

  posx = x # replicate(1,ypts)
  posy = replicate(1,xpts) # y

  vectorplot, velx, vely, posx, posy, /over, color = color('verygray'), $
    minmag = vector.min_vector, maxmag = vector.max_vector, $
    typvel = vector.typical_vector, $
    xskip = vector.xskip, yskip = vector.yskip, legend = [.70,.03], $
    outline = vector.outline

  undefine, velx
  undefine, vely

endif

if (particle_widget.enabled) then begin
  partvelvec, particles, PARTICLE_WIDGET=particle_widget, $
    /OVER, COLOR = color('verygray')
endif

if options.blocks then draw_blocks, color('ltblue'), TREE = tree1, PARAMETERS = params1, $
  ORIENTATION = problem.orientation

if (options.blocks AND var_index('NEVER') GT 0) then begin

    tstep_set = courant()

    case problem.orientation of
        0: begin
            plots, [tree[tstep_set].bndBox[0,0], $
                    tree[tstep_set].bndBox[1,0]], $
              [tree[tstep_set].bndBox[0,1], $
               tree[tstep_set].bndBox[1,1]], $
              color = color('red'), noclip = 0
            plots, [tree[tstep_set].bndBox[0,0], $
                    tree[tstep_set].bndBox[1,0]], $
              [tree[tstep_set].bndBox[1,1], $
               tree[tstep_set].bndBox[0,1]], $
              color = color('red'), noclip = 0
        end
        1: begin
            plots, [tree[tstep_set].bndBox[0,1], $
                    tree[tstep_set].bndBox[1,1]], $
              [tree[tstep_set].bndBox[0,0], $
               tree[tstep_set].bndBox[1,0]], $
              color = color('red'), noclip = 0
            plots, [tree[tstep_set].bndBox[1,1], $
                    tree[tstep_set].bndBox[0,1]], $
              [tree[tstep_set].bndBox[0,0], $
               tree[tstep_set].bndBox[1,0]], $
              color = color('red'), noclip = 0
        end
        2: begin
            plots, [tree[tstep_set].bndBox[0,0], $
                    tree[tstep_set].bndBox[1,0]], $
              [tree[tstep_set].bndBox[0,1], $
               tree[tstep_set].bndBox[1,1]], $
              color = color('red'), noclip = 0
            plots, [tree[tstep_set].bndBox[0,0], $
                    tree[tstep_set].bndBox[1,0]], $
              [tree[tstep_set].bndBox[1,1], $
               tree[tstep_set].bndBox[0,1]], $
              color = color('red'), noclip = 0
        end
    endcase
endif

;------------------
; restore the axes
;------------------

case problem.orientation of 
    0: begin
        if (problem.name EQ 'xray_cyl') then begin
            case options.tick of 
                0: begin        ; no tick marks
                end
                1: begin        ; display tick marks
                    axis, xmin, ymin, xaxis = 0, xstyle = 1, $
                      xtitle = 'r (cm)'
                    axis, xmin, ymin, yaxis = 0, ystyle = 1, $
                      ytitle = 'z (cm)'
                    axis, xmin, ymax, xaxis = 1, xstyle = 1, $
                      xtickformat = 'nolabel'
                    axis, xmax, ymin, yaxis = 1, ystyle = 1, $
                      ytickformat = 'nolabel'
                end
            endcase
        endif else if (problem.name EQ 'rt') then begin
            case options.tick of 
                0: begin        ; no tick marks
                end
                1: begin        ; display tick marks
                    axis, xmin, ymax, xaxis = 0, xstyle = 1, $
                      xtitle = 'x (cm)'
                    axis, xmin, ymax, yaxis = 0, ystyle = 1, $
                      ytitle = 'y (cm)'
                    axis, xmin, ymin, xaxis = 1, xstyle = 1, $
                      xtickformat = 'nolabel'
                    axis, xmax, ymax, yaxis = 1, ystyle = 1, $
                      ytickformat = 'nolabel'
                end
            endcase
        endif else begin
            case options.tick of 
                0: begin        ; no tick marks
                end
                1: begin        ; display tick marks
                    axis, xmin, ymin, xaxis = 0, xstyle = 1, $
                      xtitle = 'x (cm)'
                    axis, xmin, ymin, yaxis = 0, ystyle = 1, $
                      ytitle = 'y (cm)'
                    axis, xmin, ymax, xaxis = 1, xstyle = 1, $
                      xtickformat = 'nolabel'
                    axis, xmax, ymin, yaxis = 1, ystyle = 1, $
                      ytickformat = 'nolabel'
                end
            endcase
        endelse
    end
    1: begin
        case options.tick of 
            0: begin            ; no tick marks
            end
            1: begin            ; display tick marks
                axis, ymin, xmin, xaxis = 0, xstyle = 1, $
                  xtitle = 'y (cm)'
                axis, ymin, xmin, yaxis = 0, ystyle = 1, ytitle = 'x (cm)'
                axis, ymin, xmax, xaxis = 1, xstyle = 1, $
                  xtickformat = 'nolabel'
                axis, ymax, xmin, yaxis = 1, ystyle = 1, $
                  ytickformat = 'nolabel'
            end
        endcase
    end
    2: begin
        case options.tick of 
            0: begin            ; no tick marks
            end
            1: begin            ; display tick marks
                axis, xmin, ymax, xaxis = 0, xstyle = 1, xtitle = 'x (cm)'
                axis, xmin, ymax, yaxis = 0, ystyle = 1, ytitle = 'y (cm)'
                axis, xmin, ymin, xaxis = 1, xstyle = 1, $
                  xtickformat = 'nolabel'
                axis, xmax, ymax, yaxis = 1, ystyle = 1, $
                  ytickformat = 'nolabel'
            end
        endcase
    end
endcase
    
; save the plot window status so we can play with it later
p_hydro = !p & x_hydro = !x & y_hydro = !y


;-----------------------
; add a colorbar legend if desired
;-----------------------
if options.colorbar EQ 1 then begin    
    if aspect LT 1 then begin
        
        case options.log of
            
            0: colorbar2, float(variable1.min[0]), float(variable1.max[0]), $
              legend_pos, $
              COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5
            
            1: colorbar2, alog10(variable1.min[0]), alog10(variable1.max[0]), $
              legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, $
              CHARSIZE = 1.5
        endcase
        
    endif else begin

        case options.log of
            
            0: vcolorbar, float(variable1.min[0]), float(variable1.max[0]), $
              legend_pos, $
              COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5
            
            1: vcolorbar, alog10(variable1.min[0]), alog10(variable1.max[0]), $
              legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, $
              CHARSIZE = 1.5
            
        endcase
    endelse
endif
    
; restore the plot to the hydro plot so we can play with it
!p = p_hydro & !x = x_hydro & !y = y_hydro


; --------------------------------------------------------
;  finally add some information to the bottom of the plot
; --------------------------------------------------------
   
gridtxt = 'number of blocks = ' + string(format = '(i7)', params1.totBlocks)
gridtxt2 = 'AMR levels = ' + string(format = '(i5)', max(tree1[*].lrefine))

; if desired, print the time, filename, and credits
    
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

    if (n_elements(credit) NE 0) then begin
        
        if output.type EQ 1 then begin
            xyouts, .45, .05, credit, color = color('black'), /normal   
        endif else begin
            xyouts, .45, .05, credit,  color = color('black'), /normal   
        endelse
    endif

endif

; empty the buffer
empty

;if (plot_number EQ NX*NY-1 OR ifile EQ file1Info.endSuffix) then begin

case output.type of
    1: begin
        device, /close
        set_plot, current_device
    end
    2: color_bitmap, outfile
    else:
endcase
;endif

if (n_elements(info) NE 0) then $
  widget_control, info.status, $
  set_value = 'status: awaiting orders . . .'

;endfor


end


