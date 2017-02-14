;-----------------------------------------------------------------------------
; xplot2d_amr.pro
;
;    xplot2d_amr serves as the backend for the xflash3 widget.  This 
;    function will plot 2 or 3-d FLASH data, using the parameters 
;    feed in from xflash_new.
;
;    general flow:
;       -- read in data 
;       -- convert data from block structure to 2d array
;       -- scale data to colormap and display
;       -- plot block boundaries if desired
;       -- plot velocity vectors, if desired
;
;-----------------------------------------------------------------------------

pro xplot2d_amr, FILE_INFO = fileInfo, $
                 READ = iread, $
                 VARIABLE_INFO = variable, $
                 OPTIONS = options, $
                 OUTPUT = output, $
                 VECTOR = vector, $
                 ZOOM = zoom, $
                 COLORMAP = colormap, $
                 WIDGET_INFORMATION = info, $
                 CONTOUR_OPT = contours, $
                 PARTICLE_WIDGET = particle_widget, $
                 LABEL_OPT = label_info, $
                 PROBLEM_INFO = problem, $
                 KNOWN_VARIABLES = varnames, $
                 TREE_PTR = tree_ptr, $
                 PARAMS_PTR = params_ptr, $
                 DATA_PTR = data_ptr, $
                 PLOT_COUNT = plot_count, $
                 OUTFILE_NAME = outfilename, $
                 CREDIT = credit, $
                 DEBUG = debug
; common in the variables from the read-in routine

; create a common to store the max variable quantity.  During the loop
; over the files, it is stored here
common maxvar, max_var

common save_state, params, tree, unk, particles
common save, filename, pblm_ctr_lim, orientation, ndim, geometry

common plot_dat_local, xsize_old, ysize_old

if N_ELEMENTS(xsize_old) EQ 0 then xsize_old=0
if N_ELEMENTS(ysize_old) EQ 0 then ysize_old=0

; load the colortable -- xflash uses special colortables, the first 12
; colors are reserved and constant across the different tables
xflash_dir = get_xflash_path()

; get the colortable
iclrmap = color_index(colormap, MIN_VALUE=colorMin, MAX_VALUE=colorMax)
loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT

; set the plot title into titleName
set_plot_title,TITLE=titleName,VARIABLE_INFO=variable,OPTIONS=options

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
fileBase = strmid(fileInfo.prototype,0,strlen(fileInfo.prototype)-4)

pathEnd = strpos(fileBase, '/', /REVERSE_SEARCH)
fileBaseClean = strcompress(strmid(fileBase, pathEnd+1, strlen(fileBase)-pathEnd))

; for particle trajectories
if (particle_widget.traj) then begin
  tmp = (fileInfo.endSuffix - fileInfo.startSuffix)/fileInfo.step
  print, "in xplot2d_amr, tmp is: ", tmp
  if (particle_widget.traj_color) then begin
    old_parts = fltarr(3, tmp+fileInfo.step)
  endif else begin
    old_parts = fltarr(2, tmp+fileInfo.step)
  endelse
endif

;==============================================================================
; loop over the files
;==============================================================================
counter = -1
for ifile = fileInfo.startSuffix, fileInfo.endSuffix, fileInfo.step do begin

  counter++

  ; set the current plot # and the x and y position (in plotspace)
  ; of the current plot on the page -- 1 based indexing here
  plot_number = !p.multi[0]
  if (plot_number EQ 0) then plot_number = NX*NY
  plot_number = NX*NY - plot_number 

  xplot_pos = ( (plot_number) mod NX ) + 1
  yplot_pos = (plot_number - xplot_pos + 1)/NX + 1

  if Keyword_Set(debug) then print, '## plot info:', plot_number, xplot_pos, yplot_pos

  ; ---- determine the filename -------------------------------------------------
  filename = fileBase + string(ifile, format = '(i4.4)')

  if Keyword_Set(debug) then print, 'working on ' + filename

  ; later on, we'll figure out whether the file has been read already,
  ; or if we need to unzip it.  For now read the file
  itype = determine_file_type(filename)

  if (itype EQ -1) then begin
    result = dialog_message('Error: file does not exist', $
                            /error, title = 'ERROR!!')
    return
  endif

  ; only read the data if it hasn't changed since the last time
  if (iread EQ 1) then begin
    if (n_elements(info) NE 0) then widget_control, info.status, $
      set_value = 'status: reading from file'

    ; --- read the data -----------------------------------------------------------
    read_amr, filename, TREE=tree, DATA=unk, PARAMETERS=params, PARTICLES=particles, $
              GEOMETRY=geometry
    tree_ptr = tree
    data_ptr = unk
    params_ptr = params

    if (n_elements(info) NE 0) then widget_control, info.status, $
      set_value = 'status: done reading'

  endif else begin
    if Keyword_set(debug) then print, '^^^^ skipping the read'
  endelse


  ; dump a few statistics on the screen
  if Keyword_Set(debug) then print, '   number of blocks:           ', params.totBlocks

  lrefine_max = max(tree[*].lrefine)
  top_level = (size(where(tree[*].lrefine EQ lrefine_max)))[1]
  if Keyword_Set(debug) then  print, '   number of top level zones: ', top_level*params.nxb*params.nyb

  uniform = 0l
  uniform = params.ntopx*params.nxb*2^(lrefine_max-1)* $
            params.ntopy*params.nyb*2^(lrefine_max-1)

  if Keyword_Set(debug) then     print, '   fraction of uniform grid:   ', $
    float(top_level*params.nxb*params.nyb)/float(uniform)

  if params.time LT 1.e-9 then begin
    timeout  = 'time = ' + $
               string(format = '(f9.5)',params.time*1e12) + ' ps'

  endif else if params.time LT 1.e-6 then begin
    timeout  = 'time = ' + $
               string(format = '(f9.5)',params.time*1e9) + ' ns'

  endif else if params.time LT 1.e-3 then begin
    timeout  = 'time = ' + $
               string(format = '(f9.5)',params.time*1e6) + ' !4l!3s'

  endif else if params.time LT 3600. then begin
    timeout  = 'time = ' + $
               string(format = '(f9.5)',params.time) + ' s'

  endif else if params.time LT 86400. then begin
    timeout  = 'time = ' + $
      string(format = '(f9.5)',params.time/3600.) + ' hours'
    
  endif else if params.time GT 1.0e15 then begin
    timeout  = 'time = ' + $
      string(format = '(f7.3)',params.time/3.1557e16) + ' Gyr'

  endif else if params.time GT 1.0e12 then begin
      timeout  = 'time = ' + $
        string(format = '(f7.3)',params.time/3.1557e13) + ' Myr'
  endif else begin
    timeout  = 'time = ' + $
               string(format = '(f9.5)',params.time/86400.) + ' days'

  endelse

  ; ---- create the output filename ---------------------------------------------

  outfile = fileBaseClean

  if options.blocks then outfile = outfile + 'blk_'

  if n_elements(outfilename) EQ 0 then begin

    outfile = outfile + strcompress(variable.name, /REMOVE_ALL) + $
              string(ifile, format = '(i4.4)')

    case output.type of
      1: outfile = outfile + '.ps'
      2: begin
        ; newer version of IDL no longer support Gifs
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

  print, 'outputting to ', outfile

  ;------------------------------------------------------------------------------
  ; load the variable to be plotted into temporary storage
  ;------------------------------------------------------------------------------
  if Keyword_Set(debug) then    print, 'variable name = ', variable.name
  temp_arr = create_variable(unk, variable.name,debug)

  if options.abs EQ 1 then temp_arr = abs(temp_arr)

  if Keyword_Set(debug) then print, variable.name + ' extrema: ', min(temp_arr), max(temp_arr)

  ;------------------------------------------------------------------------------
  ; find the size of the plotting area
  ;------------------------------------------------------------------------------

  if (params.geometry EQ "CARTESIAN" OR $
      params.geometry EQ "CYLINDRICAL") then begin

    ; by default, set the plot extrema to the domain limits.  If this was
    ; overridden by the zoom boxes, then make sure that the value given
    ; for the limit is within the computational domain.

    xmin_tot = min(tree[*].bndBox[0,0])
    xmax_tot = max(tree[*].bndBox[1,0])        

    ymin_tot = min(tree[*].bndBox[0,1])
    ymax_tot = max(tree[*].bndBox[1,1])

    ; right now, a value of -1.0 means use the default.  This should be
    ; changed at some point to allow any valid negative coordinate.

    if (zoom.xmin NE -1.d0) then begin
      xmin = (zoom.xmin > xmin_tot) < max(tree[*].bndBox[0,0])
      xmin = (xmin)[0]
    endif else begin
      xmin = xmin_tot
    endelse

    if (zoom.ymin NE -1.0d0) then begin
      ymin = (zoom.ymin > ymin_tot) < max(tree[*].bndBox[0,1])
      ymin = (ymin)[0]
    endif else begin
      ymin = ymin_tot
    endelse


    if (zoom.xmax NE -1.d0) then begin
      xmax = (zoom.xmax < xmax_tot) > xmin*1.00001
      xmax = (xmax)[0]
    endif else begin
      xmax = xmax_tot
    endelse

    if (zoom.ymax NE -1.d0) then begin
      ymax = (zoom.ymax < ymax_tot) > ymin*1.00001
      ymax = (ymax)[0]
    endif else begin
      ymax = ymax_tot
    endelse

  endif else begin

    theta_min = min(tree[*].bndBox[0,1])
    theta_max = max(tree[*].bndBox[1,1])

    ; temporary hack
    if zoom.ymin EQ -1 then zoom.ymin = theta_min
    if zoom.ymax EQ -1 then zoom.ymax = theta_max

    theta_min = theta_min > zoom.ymin
    theta_max = theta_max < zoom.ymax

    dtheta = theta_max - theta_min

    ; for now, assume that we start at r = 0
    r_min = min(tree[*].bndBox[0,0])
    r_max = max(tree[*].bndBox[1,0])

    if zoom.xmin EQ -1 then zoom.xmin = r_min
    if zoom.xmax EQ -1 then zoom.xmax = r_max

    r_min = zoom.xmin > r_min
    r_max = zoom.xmax < r_max

    if Keyword_Set(debug) then  print, 'theta range: ', theta_min, theta_max
    if Keyword_Set(debug) then  print, 'r range: ', r_min, r_max

    ; figure out how big of a uniformly gridded data array we need.  For
    ; simplicity, we will always put the theta = 0 vertically.  Depending
    ; on how many quadrants were span, our x and y extrema change

    if (dtheta LT !pi/2. - 1.e-4) then begin
      if Keyword_Set(debug) then  print, dtheta, !pi/2.

      ; we are less than one quadrant

      if (theta_max LE !pi/2) then begin
        xmin = 0.
        ymin = 0.

        xmax = r_max*sin(theta_max)
        ymax = r_max*cos(theta_min)

      endif else if (theta_max LT !pi) then begin
        xmin = 0.
        ymax = 0.

        xmax = r_max*sin(theta_min)
        ymin = r_max*cos(theta_max)

      endif else begin
        print, 'ERROR: angle range not yet suported'
      endelse
    endif else begin

      ; we span more than one quadrant
      if ( theta_max GT !pi ) then begin

        xmin = -r_max
        xmax =  r_max

        ymin = -r_max
        ymax =  r_max

      endif else begin

        xmin = 0.d0
        xmax = r_max

        ymin = r_max*cos(theta_max)
        ymax = r_max*cos(theta_min)

      endelse
    endelse
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

  if (params.geometry EQ "CARTESIAN" OR $
      params.geometry EQ "CYLINDRICAL") then begin

    ; compute the fraction of the total domain that is being plotted
    xfrac = dx/(xmax_tot - xmin_tot)
    yfrac = dy/(ymax_tot - ymin_tot)
  endif


  ; ---- initialize the device --------------------------------------------------

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
        ysize = 10.0
        if (plot_number EQ 0) then begin
          device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
                  XOFF = 0.5, YOFF = 0.5, /INCH, /COLOR, $
                  BITS_PER_PIXEL = 8, LANDSCAPE=0
        endif
      end

      ; landscape orientation
      1: begin
        xsize = 10.
        ysize = 7.5

        if (plot_number EQ 0) then begin
          device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
                  XOFF = 0.5, YOFF = 10.5, /INCH, /COLOR, $
                  BITS_PER_PIXEL = 8, /LANDSCAPE
        endif
      end
    endcase

    deviceAspect = (ysize)/(xsize)

  endif else begin
    deviceAspect = (float(output.vsize))/ $
                   (float(output.hsize))
  endelse 

  ; compute the aspect ratio of the domain and window (paper) together
  aspect = dataAspectRatio/deviceAspect

  ; ---- determine the bounds for the plot --------------------------------------

  ; set the normal coordinates of the portion of the display/page you
  ; wish to use -- leave a little margin so the plots don't run to the
  ; edge of the page
  if (multi_plots_per_page EQ 0) then begin
    page_nx1 = 0.05
    page_nx2 = 0.95

    page_ny1 = 0.05
    page_ny2 = 0.90
  endif else begin

    if aspect GT 1 then begin
      page_nx1 = 0.05
      page_nx2 = 0.85

      page_ny1 = 0.05
      page_ny2 = 0.90
    endif else begin
      page_nx1 = 0.05
      page_nx2 = 0.95

      page_ny1 = 0.15
      page_ny2 = 0.90
    endelse
  endelse

  ; compute the bounding normal coordinates for this subplot
  ; the current subplot occupies the rectangle with opposite
  ; corners (nx1, ny1) to (nx2, ny2).  (note, 0,0 is in the 
  ; lower left corner in normal coords).

  dpagex = page_nx2 - page_nx1
  dpagey = page_ny2 - page_ny1

  if multi_plots_per_page EQ 0 then begin
    nx1 = (dpagex/float(NX)) * (xplot_pos - 1 + 0.15) + page_nx1
    nx2 = (dpagex/float(NX)) * (xplot_pos - 0.05)  + page_nx1

    ny1 = page_ny2 - (dpagey/float(NY)) * (yplot_pos - 0.15) 
    ny2 = page_ny2 - (dpagey/float(NY)) * (yplot_pos - 1 + 0.05) 
  endif else begin
    nx1 = (dpagex/float(NX)) * (xplot_pos - 1 + 0.15) + page_nx1
    nx2 = (dpagex/float(NX)) * (xplot_pos - 0.05)  + page_nx1

    ny1 = page_ny2 - (dpagey/float(NY)) * (yplot_pos - 0.15) 
    ny2 = page_ny2 - (dpagey/float(NY)) * (yplot_pos - 1 + 0.05) 
  endelse

  ; aspect compares the aspect ratio of the data to that of the device,
  ; aspect = (dy_data/dy_device)/(dx_data/dx_device).  If aspect > 1, 
  ; it means that the y-coordinate is the one that is going to set the
  ; overall scaling.  If aspect < 1, then the x-coordinate will set the
  ; scaling.  Consider each case separately below.

  if (multi_plots_per_page EQ 0) then begin
    cb_factor = 0.7
  endif else begin
    cb_factor = 0.90
  endelse

  if (aspect GE 1.) then begin 

    ; the y size of the data sets the scaling of the plot

    ; set the initial values of the y min and max normal coordinates.  We
    ; leave some room for the axis labels
    py1 = ny1
    py2 = ny2

    ; compute the x size, using the aspect ratio
    dpy = py2 - py1
    dpx = dpy/aspect < cb_factor * (dpagex/float(NX))

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
    dpy = aspect*dpx < cb_factor * (dpagey/float(NY))

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
  if (ifile EQ fileInfo.startSuffix) then begin

    if (output.type NE 1) then begin
      if (output.hsize NE xsize_old) OR (output.vsize NE ysize_old) then begin
        window, 0, XSIZE = output.hsize, YSIZE = output.vsize, $
        TITLE = 'AMR plot'

        xsize_old=output.hsize
        ysize_old=output.vsize
      endif
    endif

    if Keyword_Set(debug) then print, 'main plot window = ', !D.WINDOW

  endif

  ;-----------------------------------------------------------------------------
  ; create the axes -- just get the data coordinates defined
  ;-----------------------------------------------------------------------------

  ; create the title -- if we are printing it
  if (multi_plots_per_page EQ 0) then begin

    if (particle_widget.data_enabled) then begin
      title = 'Particles by Temp'
    endif else if (options.annotate EQ 1) then begin
      title = titleName
    endif else begin
      title = ' '
    endelse

  endif else begin
    title = timeout
  endelse


  if Keyword_Set(debug) then print, xmin,xmax,ymin,ymax

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
  if (multi_plots_per_page EQ 1 AND plot_number EQ 0) then begin
    xyouts, .5, 0.5*(1.0 + page_ny2), titleName, $
            color = color('black'), /normal, alignment = 0.5, charsize=2
  endif

  if options.procdist eq '' then begin
    ;------------------------------------------------------------------------------
    ; determine whether to sub-sample the data
    ;------------------------------------------------------------------------------

    if (params.geometry EQ "CARTESIAN" OR $
        params.geometry EQ "CYLINDRICAL") then begin

      ; assume 1024 pixels on screen and 2048 on paper (~ 300dpi

      sample = 0

      max_refine = max(tree[*].lrefine)

      ; compute the number of pixels in a uniform grid of our sub-domain
      max_pixels = params.ntopx*params.nxb*2.^(max_refine-1)*xfrac > $
                   params.ntopy*params.nyb*2.^(max_refine-1)*yfrac

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
      if Keyword_Set(debug) then print, 'sub-sampling by ', sample, ' levels of refinement'

    endif

    ;------------------------------------------------------------------------------
    ; scale the data
    ;------------------------------------------------------------------------------
    if (n_elements(info) NE 0) then $
      widget_control, info.status, SET_VALUE = 'status: scaling data'

    time_start = systime(/seconds)

    if Keyword_Set(debug) then help, temp_arr
    if Keyword_Set(debug) then print, max(temp_arr)

    if (params.geometry EQ "CARTESIAN" OR $
        params.geometry EQ "CYLINDRICAL") then begin

      ; first put it on a uniform grid
      nsdata = merge_amr(temp_arr, TREE = tree, PARAMETERS = params, $
                         XMERGE = x, YMERGE = y, SAMPLE = sample, $
                         XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)

      ; if we are automatically scaling, reset the data range
      if variable.auto EQ 1 then begin
        variable.min = min(nsdata)
        variable.max = max(nsdata)
      endif

    endif else begin
      nsdata = merge_polar(temp_arr, TREE = tree, PARAMETERS = params, $
                           XMERGE = x, YMERGE = y, SAMPLE = sample, $
                           RRANGE = [r_min,r_max], TRANGE = [theta_min,theta_max])

      index_valid_data = where(nsdata NE -1.e33)

      if variable.auto EQ 1 then begin
        variable.min = min(nsdata[index_valid_data])
        variable.max = max(nsdata[index_valid_data])
      endif

    endelse

    if Keyword_Set(debug) then print, 'variable plot range: ', variable.min, variable.max

    if options.log EQ 1 then nsdata = alog10(temporary(nsdata))

    ; if we are max'ing the variable, treat it specially
    if options.max EQ 0 then begin
      case options.log of
        0: nsdata = scale_color(temporary(nsdata), $
                                VARMAX = variable.max, $
                                VARMIN = variable.min, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)

        1: nsdata = scale_color(temporary(nsdata), $
                                VARMAX = variable.max, $
                                VARMIN = variable.min, /log, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)
      endcase

      time_end = systime(/seconds)

      if Keyword_Set(debug) then print, '**** time to merge and scale = ', time_end - time_start

    endif else begin

      ; take the maximum
      if ifile EQ fileInfo.startSuffix then begin
        max_var = nsdata
        if Keyword_Set(debug) then print, 'max = ', max(max_var), min(max_var)
      endif else begin
        max_var = nsdata > max_var
      endelse

      ; scale it
      case options.log of
        0: nsdata = scale_color(max_var, $
                                VARMAX = variable.max, $
                                VARMIN = variable.min, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)

        1: nsdata = scale_color(max_var, $
                                VARMAX = variable.max, $
                                VARMIN = variable.min, /log, $
                                COLORMAP_MIN = colorMin, $
                                COLORMAP_MAX = colorMax)
      endcase
    endelse 

    if problem.orientation EQ 1 then nsdata = transpose(temporary(nsdata))

    ;-----------------------------------------------------------------------------
    ; setup particle stuff
    ;-----------------------------------------------------------------------------

    ; the goal of this is to plot particles by coordinates colored according to
    ; temperature (hopefully more potential color values to come).  this can be 
    ; done over an existing hydro plot, or over a dummy colormap (set to the minimum
    ; value of the current colormap)
    ; DEV but you need to check what data exists in the particles
    ; structure DOH

    if (particle_widget.data_enabled) then begin

      parts_temp = particles.ptemp
      case options.log of
        0: parts_temp = scale_color(parts_temp, $
                                    VARMAX = variable.max, $
                                    VARMIN = variable.min, $
                                    COLORMAP_MIN = colorMin, $
                                    COLORMAP_MAX = colorMax)
        1: parts_temp = scale_color(alog10(parts_temp), $
                                    VARMAX = variable.max, $
                                    VARMIN = variable.min, /log, $
                                    COLORMAP_MIN = colorMin, $
                                    COLORMAP_MAX = colorMax)
      endcase

      if (particle_widget.dummy) then begin
        nsdata[*,*] = colorMin
      endif 

    endif else begin
      parts_temp = 1.
    endelse


    ;------------------------------------------------------------------------------
    ; plot the data
    ;------------------------------------------------------------------------------
    if (n_elements(info) NE 0) then widget_control, info.status, $
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
      /OVERPLOT, /HALF_HALF

    ;------------------------------------------------------------------------------
    ; plot any user-defined contours
    ;------------------------------------------------------------------------------

    if Keyword_Set(debug) then print, 'contour 0: ', contours[0]

    for i = 0, (size(contours))[1]-1 do begin

      if (contours[i].enabled) then begin

        if (n_elements(info) NE 0) then widget_control, info.status, $
          set_value = 'status: plotting contour #' + string(i,'(i2)')

        if (params.geometry EQ "CARTESIAN" OR $
          params.geometry EQ "CYLINDRICAL") then begin
          sctr = merge_amr(reform(unk(contours[i].var,*,*,*,*)), $
                           TREE = tree, PARAMETERS = params, $
                           XMERGE = x, YMERGE = y, sample = sample, $
                           XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)
        endif else begin

          sctr = merge_polar(reform(unk(contours[i].var,*,*,*,*)), $
                             TREE = tree, PARAMETERS = params, $
                             XMERGE = x, YMERGE = y, sample = sample, $
                             RRANGE = [r_min,r_max], TRANGE = [theta_min,theta_max])
        endelse

        if Keyword_Set(debug) then print, 'sctr: min, max: ', min(sctr), max(sctr)
        if problem.orientation EQ 1 then begin
          sctr = transpose(temporary(sctr))
          tempy = y
          y = x
          x = tempy
          undefine, tempy
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

      if (n_elements(info) NE 0) then widget_control, info.status, $
        set_value = 'status: plotting vectors'

      case problem.orientation of
        0: begin
          if (params.geometry EQ "CARTESIAN" OR $
            params.geometry EQ "CYLINDRICAL") then begin

            velx = merge_amr(velx, TREE = tree, PARAMETERS = params, $
                             sample = sample, $
                             XMERGE = x, YMERGE = y, $
                             XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)

            vely = merge_amr(vely, TREE = tree, PARAMETERS = params, $
                             sample = sample, $
                             XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)
          endif else begin
            velr = merge_polar(velx, TREE = tree, PARAMETERS = params, $
                               sample = sample, $
                               XMERGE = x, YMERGE = y, $
                               RRANGE = [r_min,r_max], TRANGE = [theta_min,theta_max])

            velt = merge_polar(vely, TREE = tree, PARAMETERS = params, $
                               sample = sample, $
                               RRANGE = [r_min,r_max], TRANGE = [theta_min,theta_max])
          endelse

        end

        1: begin
          if (params.geometry EQ "CARTESIAN" OR $
            params.geometry EQ "CYLINDRICAL") then begin
            vely = merge_amr(velx, TREE = tree, PARAMETERS = params, $
                             sample = sample, $
                             XMERGE = x, YMERGE = y, $
                             XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)

            vely = transpose(temporary(vely))

            velx = merge_amr(vely, TREE = tree, PARAMETERS = params, $
                             sample = sample, $
                             XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)

            velx = transpose(temporary(velx))
          endif else begin
            print, 'not sure what orientation = 1 means in polar coords'
          endelse
        end

        2: begin
          if (params.geometry EQ "CARTESIAN" OR $
            params.geometry EQ "CYLINDRICAL") then begin
            velx = merge_amr(velx, TREE = tree, PARAMETERS = params, $
                             sample = sample, $
                             XMERGE = x, YMERGE = y, $
                             XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)
                        
            vely = merge_amr(vely, TREE = tree, PARAMETERS = params, $
                             sample = sample, $
                             XRANGE = [xmin,xmax], YRANGE = [ymin,ymax], DEBUG=debug)
          endif else begin
            print, 'not sure what orientation = 2 means in polar coords'
          endelse
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

      ; if we are in spherical coords, velr is actually the radial velocity,
      ; and velt is actually the theta velocity.  We need to make these into
      ; x and y velocities for the velocity vector plotter.

      ; the main thing that we need is an array of the theta for each
      ; element of velr and velt
      if (params.geometry EQ "SPHERICAL") then begin
        theta = atan(posx, posy)

        velx = velr*sin(theta) + velt*cos(theta)
        vely = velr*cos(theta) - velt*sin(theta)
      endif

      vectorplot, velx, vely, posx, posy, /over, color = color('verygray'), $
        minmag = vector.min_vector, maxmag = vector.max_vector, $
        typvel = vector.typical_vector, $
        xskip = vector.xskip, yskip = vector.yskip, legend = [.70,.03], $
        outline = vector.outline

      undefine, velx
      undefine, vely

    endif  ; vector.enabled == true

    if (particle_widget.enabled) then begin
      if (particle_widget.traj) then begin
        if counter EQ 0 then begin
          tags = particles.tag
          numParticles = n_elements(tags) ; one tag for each particle means
                                          ; number of tags = number of particles
          soughtParticle = particle_widget.partnum

          if (soughtParticle LE 0) then begin
            msg1 = "Particle number must be 1 or greater"
            result = dialog_message(msg1)
            break
          endif else if (soughtParticle GT numParticles) then begin
            msg1 = strjoin(["Chosen particle number ", strtrim(soughtParticle,1), " exceeds"])
            msg2 = strjoin(["the total number of particles ", strtrim(numParticles,1), "."])
            msg3 = ""
            msg4 = "Use 'Particle Options' button to adjust."
            result = dialog_message ([msg1, msg2, msg3, msg4])
            break
          endif else begin
            soughtParticle = particle_widget.partnum
            particleTag = tags[soughtParticle-1]
          endelse
        endif
      endif

      partvelvec, particles, PARTICLE_WIDGET = particle_widget, /OVER, $
        COLOR=color('verygray'), COUNTER=counter, PARTICLE_TAG=particleTag, $
        OLD_PARTS=old_parts, DIR=xflash_dir, PART_TEMP = parts_temp
    endif

    if (particle_widget.enabled AND particle_widget.traj_color) then begin
      iclrmap = color_index(colormap, MIN_VALUE=colorMin, MAX_VALUE=colorMax)
      loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT
    endif

    if (params.geometry EQ "CARTESIAN"  OR $
      params.geometry EQ "CYLINDRICAL") then begin
      if options.blocks then begin
        if Keyword_Set(debug) then print, 'drawing just blocks'
        draw_blocks, color('black'), TREE = tree, PARAMETERS = params, $
          ORIENTATION = problem.orientation

      endif
    endif else begin
      if Keyword_Set(debug) then  print, 'polar?'
      draw_blocks_polar, color('ltblue'), TREE = tree, PARAMETERS = params, $
        ORIENTATION = problem.orientation
    endelse

  endif else begin

    ; just show the distribution of blocks on processors
    if (params.geometry EQ "CARTESIAN" OR $
      params.geometry EQ "CYLINDRICAL") then begin

      if Keyword_set(debug) then print, 'procdist show!'
      show_proc_distribution, LEVEL=options.procdist, TREE = tree, PARAMETERS = params, $
      ORIENTATION = problem.orientation, SHOW_NUMBERS = options.showprocnums

      draw_blocks, color('black'), LEVEL=options.procdist, TREE = tree, PARAMETERS = params, $
        ORIENTATION = problem.orientation

    endif

  endelse

  if ((options.procdist OR options.blocks) AND var_index('NEVER') GT 0) then begin

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
          0: begin    ; no tick marks
          end
          1: begin    ; display tick marks
            axis, xmin, ymin, xaxis = 0, xstyle = 1, xtitle = 'r (cm)'
            axis, xmin, ymin, yaxis = 0, ystyle = 1, ytitle = 'z (cm)'
            axis, xmin, ymax, xaxis = 1, xstyle = 1, xtickformat = 'nolabel'
            axis, xmax, ymin, yaxis = 1, ystyle = 1, ytickformat = 'nolabel'
          end
        endcase
      endif else if (problem.name EQ 'rt') then begin
        case options.tick of 
          0: begin    ; no tick marks
          end
          1: begin    ; display tick marks
            axis, xmin, ymax, xaxis = 0, xstyle = 1, xtitle = 'x (cm)'
            axis, xmin, ymax, yaxis = 0, ystyle = 1, ytitle = 'y (cm)'
            axis, xmin, ymin, xaxis = 1, xstyle = 1, xtickformat = 'nolabel'
            axis, xmax, ymax, yaxis = 1, ystyle = 1, ytickformat = 'nolabel'
          end
        endcase
      endif else begin
        case options.tick of 
          0: begin    ; no tick marks
          end
          1: begin    ; display tick marks
            axis, xmin, ymin, xaxis = 0, xstyle = 1, xtitle = 'x (cm)'
            axis, xmin, ymin, yaxis = 0, ystyle = 1, ytitle = 'y (cm)'
            axis, xmin, ymax, xaxis = 1, xstyle = 1, xtickformat = 'nolabel'
            axis, xmax, ymin, yaxis = 1, ystyle = 1, ytickformat = 'nolabel'
          end
        endcase
      endelse
    end

    1: begin
      case options.tick of 
        0: begin        ; no tick marks
        end
        1: begin        ; display tick marks
          axis, ymin, xmin, xaxis = 0, xstyle = 1, xtitle = 'y (cm)'
          axis, ymin, xmin, yaxis = 0, ystyle = 1, ytitle = 'x (cm)'
          axis, ymin, xmax, xaxis = 1, xstyle = 1, xtickformat = 'nolabel'
          axis, ymax, xmin, yaxis = 1, ystyle = 1, ytickformat = 'nolabel'
        end
      endcase
    end

    2: begin
      case options.tick of 
        0: begin        ; no tick marks
        end
        1: begin        ; display tick marks
          axis, xmin, ymax, xaxis = 0, xstyle = 1, xtitle = 'x (cm)'
          axis, xmin, ymax, yaxis = 0, ystyle = 1, ytitle = 'y (cm)'
          axis, xmin, ymin, xaxis = 1, xstyle = 1, xtickformat = 'nolabel'
          axis, xmax, ymax, yaxis = 1, ystyle = 1, ytickformat = 'nolabel'
        end
      endcase
    end
  endcase

  ; save the plot window status so we can play with it later
  p_hydro = !p & x_hydro = !x & y_hydro = !y


  ;-----------------------
  ; add the floating label if desired
  ;-----------------------

  if label_info.enabled then begin
    if multi_plots_per_page EQ 0 then begin
      case problem.orientation of
        0: begin
          xyouts, label_info.posx*(xmax-xmin)+xmin, label_info.posy*(ymax-ymin)+ymin,$
            label_info.label[0], alignment = .5, charsize = label_info.size, $
            charthick = label_info.thick, color = label_info.color
        end
        1: begin
          xyouts, label_info.posx*(ymax-ymin)+ymin, label_info.posy*(xmax-xmin)+xmin,$
            label_info.label[0], alignment = .5, charsize = label_info.size, $
            charthick = label_info.thick, color = label_info.color
        end
        2: begin
          xyouts, label_info.posx*(xmax-xmin)+xmin, label_info.posy*(ymin-ymax)+ymax,$
            label_info.label[0], alignment = .5, charsize = label_info.size, $
            charthick = label_info.thick, color = label_info.color
        end
      endcase

    endif else begin

      labels = strarr(NX*NY)
      count2=0
      for count=0, 8 do begin
        if (label_info.label[count] NE '') then begin
          labels[count2] = label_info.label[count]
          if labels[count2] EQ '-1' then labels[count2]=''
          count2 = count2+1
        endif
      endfor

      case problem.orientation of
        0: begin
          xyouts, label_info.posx*(xmax-xmin)+xmin, label_info.posy*(ymax-ymin)+ymin,$
            labels[plot_number], alignment = .5, charsize = label_info.size, $
            charthick = label_info.thick, color = label_info.color
        end
        1: begin
          xyouts, label_info.posx*(ymax-ymin)+ymin, label_info.posy*(xmax-xmin)+xmin,$
            labels[plot_number], alignment = .5, charsize = label_info.size, $
            charthick = label_info.thick, color = label_info.color
        end
        2: begin
          xyouts, label_info.posx*(xmax-xmin)+xmin, label_info.posy*(ymin-ymax)+ymax,$
            labels[plot_number], alignment = .5, charsize = label_info.size, $
            charthick = label_info.thick, color = label_info.color
        end
      endcase
    endelse
  endif


  ;-----------------------
  ; add a colorbar legend if desired
  ;-----------------------
  if options.colorbar EQ 1 AND multi_plots_per_page EQ 0 then begin    
    if aspect LT 1 then begin
      if (particle_widget.traj_color AND particle_widget.enabled) then begin

        legend_pos = [px1, lcenter-.65*lheight, px2, lcenter-0.15*lheight]
        iclrmap = color_index('Grayscale', MIN_VALUE=colorMin, MAX_VALUE=colorMax)
        loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT

        case particle_widget.log of
          0: colorbar2, particle_widget.color_min, particle_widget.color_max, $
            legend_pos, COLORMIN=colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

          1: colorbar2, alog10(particle_widget.color_min), alog10(particle_widget.color_max), $
            legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5
        endcase

        legend_pos = [px1, lcenter+0.5*lheight, px2, lcenter+lheight]
        iclrmap = color_index(colormap, MIN_VALUE=colorMin, MAX_VALUE=colorMax)
        loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT

      endif

      case options.log of
        0: colorbar2, float(variable.min[0]), float(variable.max[0]), $
          legend_pos, $
          COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

        1: colorbar2, alog10(variable.min[0]), alog10(variable.max[0]), $
          legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, $
          CHARSIZE = 1.5
      endcase

    endif else begin

      if (particle_widget.traj_color AND particle_widget.enabled) then begin

        legend_pos = [lcenter+.15*lwidth, py1, lcenter+.65*lwidth, py2]
        iclrmap = color_index('Grayscale', MIN_VALUE=colorMin, $
                              MAX_VALUE=colorMax)
        loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT

        case particle_widget.log of
          0: vcolorbar, particle_widget.color_min, $
            particle_widget.color_max, $
            legend_pos, COLORMIN=colorMin, COLORMAX = colorMax, $
            CHARSIZE = 1.5

          1: vcolorbar, alog10(particle_widget.color_min), $
            alog10(particle_widget.color_max), $
            legend_pos, COLORMIN=colorMin, COLORMAX = colorMax, $
            CHARSIZE = 1.5
        endcase

        legend_pos = [lcenter-1.25*lwidth, py1, lcenter-.75*lwidth, py2]
        iclrmap = color_index(colormap, MIN_VALUE=colorMin, MAX_VALUE=colorMax)
        loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT
      endif

      case options.log of

        0: vcolorbar, float(variable.min[0]), float(variable.max[0]), $
          legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

        1: vcolorbar, alog10(variable.min[0]), alog10(variable.max[0]), $
          legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

      endcase
    endelse
  endif

  ; restore the plot to the hydro plot so we can play with it
  !p = p_hydro & !x = x_hydro & !y = y_hydro

  ; --------------------------------------------------------
  ;  finally add some information to the bottom of the plot
  ; --------------------------------------------------------

  gridtxt = 'number of blocks = ' + string(format = '(i7)', params.totBlocks)
  gridtxt2 = 'AMR levels = ' + string(format = '(i5)', max(tree[*].lrefine))

  ; if desired, print the time, filename, and credits

  if (multi_plots_per_page EQ 0) then begin
    if options.annotate EQ 1 then begin
      if output.type EQ 1 then begin
        xyouts, .05, .05, timeout, color = color('black'), /normal   

        xyouts, .05, .02, gridtxt + ', ' + gridtxt2, color = color('black'), /normal

        xyouts, .05, -.05, filename, color = color('dkgray'), /normal, charsize = .8

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

  ; multiple colorbars (for particle trajectories) are not yet available with 
  ; multiple plots.  kr 08/04
  if (plot_number EQ NX*NY-1 OR ifile EQ fileInfo.endSuffix) then begin
    if options.colorbar EQ 1 AND multi_plots_per_page EQ 1 then begin    
      if aspect LT 1 then begin
        case options.log of
          0: colorbar2, float(variable.min[0]), float(variable.max[0]), $
            [0.10, 0.05, 0.90, 0.1], COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

          1: colorbar2, alog10(variable.min[0]), alog10(variable.max[0]), $
            [0.10, 0.05, 0.90, 0.1], COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5
        endcase

      endif else begin
        case options.log of
          0: vcolorbar, float(variable.min[0]), float(variable.max[0]), $
            [0.9, 0.10 ,0.95 ,0.90], COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

          1: vcolorbar, alog10(variable.min[0]), alog10(variable.max[0]), $
            [0.9, 0.10, 0.95, 0.90], COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

        endcase
      endelse
    endif

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
    widget_control, info.status, set_value = 'status: awaiting orders . . .'

endfor

end
