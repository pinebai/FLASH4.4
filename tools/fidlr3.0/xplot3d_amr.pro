;-----------------------------------------------------------------------------
; Xplot3d_amr.pro -- MZ 3-11-00
;
;   plot output from flash -- for 3d data slices
;
;   this routine is typically driven by xflash3
;
;   general flow:
;       -- read in data 
;       -- convert data from block structure to 2d slice 
;       -- scale data to colormap and display
;       -- plot block boundaries if desired
;       -- plot velocity vectors, if desired
;
;----------------------------------------------------------------------------
pro xplot3d_amr, FILE_INFO = fileInfo, $
                 VARIABLE_INFO = variable, $
                 OPTIONS = options, $
                 CORNERS = corners, $
                 OUTPUT = output, $
                 PROBLEM_INFO = problem, $
                 LABEL_OPT = label_info, $
                 SLICE_DIR = islice_dir, $
                 ZOOM = zoom, $
                 COLORMAP = colormap, $
                 WIDGET_INFORMATION = info, $
                 DATA_PTR = data_ptr, $
                 PARAMS_PTR = params_ptr, $
                 XMERGE = xmerge, $
                 YMERGE = ymerge, $
                 ZMERGE = zmerge, $
                 DEBUG = debug


common system, xflash_dir

; create a common to store the max variable quantity.  During the loop
; over the files, it is stored here
common maxvar, max_var

common save, filename, pblm_ctr_lim, orientation, ndim, geometry

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

!p.multi      = [0,1,1]


if (output.type EQ 1) then begin
  !p.charsize   = 1.0
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
fileBaseClean = strcompress(strmid(fileBase, pathEnd+1, $
                                   strlen(fileBase)-pathEnd))


;==============================================================================
; loop over the files
;==============================================================================
for ifile = fileInfo.startSuffix, fileInfo.endSuffix, fileInfo.step do begin

  ; ---- determine the filename -------------------------------------------------
  filename = fileBase + string(ifile, format = '(i4.4)')


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

  ; read in only the variable to be plotted -- save memory
  read_amr, filename, VAR_NAME=variable.name, $
    TREE=tree, DATA=unk, PARAMETERS=params, $
    GEOMETRY=geometry

  params_ptr = params

  if (n_elements(info) NE 0) then widget_control, info.status, $
    set_value = 'status: done reading'

  ; dump a few statistics on the screen
  if Keyword_Set(debug) then print, '   number of blocks:           ', params.totBlocks

  lrefine_max = max(tree[*].lrefine)
  top_level = (size(where(tree[*].lrefine EQ lrefine_max)))[1]
  if Keyword_Set(debug) then print, '   number of top level zones: ', top_level*params.nxb*params.nyb*params.nzb

  uniform = params.ntopx*params.nxb*2^(lrefine_max-1)* $
    params.ntopy*params.nyb*2^(lrefine_max-1)* $
    params.ntopz*params.nzb*2^(lrefine_max-1)

  if Keyword_Set(debug) then print, '   fraction of uniform grid:   ', $
    float(top_level*params.nxb*params.nyb*params.nzb)/float(uniform)

  if params.time LT 1.e-9 then begin
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
    
    
  ;-----------------------------------------------------------------------------
  ; load the variable to be plotted into temporary storage
  ;-----------------------------------------------------------------------------
  if options.procdist eq '' then begin
    temp_arr = drop_1st_dimension(unk(0,*,*,*,*))

    if options.abs EQ 1 then temp_arr = abs(temp_arr)

    if Keyword_Set(debug) then print, variable.name + ' extrema: ', min(temp_arr), max(temp_arr)
        
    ; if we are automagically scaling, reset the data range
    if variable.auto EQ 1 then begin
      variable.min = min(temp_arr)
      variable.max = max(temp_arr)
    endif
  endif 

  ;-----------------------------------------------------------------------------
  ; find the size of the plotting area
  ;-----------------------------------------------------------------------------

  ; by default, set the plot extrema to the domain limits.  If this was
  ; overridden by the zoom boxes, then make sure that the value given
  ; for the limit is within the computational domain.
  xmin_tot = min(tree[*].bndBox[0,0])
  xmax_tot = max(tree[*].bndBox[1,0])        

  ymin_tot = min(tree[*].bndBox[0,1])
  ymax_tot = max(tree[*].bndBox[1,1])

  zmin_tot = min(tree[*].bndBox[0,2])
  zmax_tot = max(tree[*].bndBox[1,2])


  ; right now, a value of -1.0 means use the default.  This should be
  ; changed at some point to allow any valid negative coordinate.

  ; mimimum coordinates
  if (zoom.xmin NE -1.0d0) then begin
    xmin = (zoom.xmin > xmin_tot) < xmax_tot
  endif else begin
    xmin = xmin_tot
  endelse

  if (zoom.ymin NE -1.0d0) then begin
    ymin = (zoom.ymin > ymin_tot) < ymax_tot
  endif else begin
    ymin = ymin_tot
  endelse

  if (zoom.zmin NE -1.0d0) then begin
    zmin = (zoom.zmin > zmin_tot) < zmax_tot
  endif else begin
    zmin = zmin_tot
  endelse

  ; maximum coordinates
  if (zoom.xmax NE -1.0d0) then begin
    xmax = (zoom.xmax < xmax_tot) > xmin*1.00001
  endif else begin
    xmax = xmax_tot
  endelse

  if (islice_dir EQ 2) then xmax = xmin

  if (zoom.ymax NE -1.0d0) then begin
    ymax = (zoom.ymax < ymax_tot) > ymin*1.00001
  endif else begin
    ymax = ymax_tot
  endelse

  if (islice_dir EQ 1) then ymax = ymin
    
  if (zoom.zmax NE -1.0d0) then begin
    zmax = (zoom.zmax < zmax_tot) > zmin*1.00001
  endif else begin
    zmax = zmax_tot
  endelse

  if (islice_dir EQ 0) then zmax = zmin

  ; compute the aspect ratio, and define the plot size
  ; determine the window size and 
  case islice_dir of
    0: begin  ; x-y plane
      dvertical   = ymax - ymin
      dhorizontal = xmax - xmin

      vfrac = dvertical/(ymax_tot - ymin_tot)
      hfrac = dhorizontal/(xmax_tot - xmin_tot)
    end
    1: begin  ; x-z plane
      dvertical   = zmax - zmin
      dhorizontal = xmax - xmin

      vfrac = dvertical/(zmax_tot - zmin_tot)
      hfrac = dhorizontal/(xmax_tot - xmin_tot)
    end
    2: begin  ; y-z plane
      dvertical   = zmax - zmin
      dhorizontal = ymax - ymin

      vfrac = dvertical/(zmax_tot - zmin_tot)
      hfrac = dhorizontal/(ymax_tot - ymin_tot)
    end
  endcase

  case problem.orientation of
    0: dataAspectRatio = dvertical/dhorizontal
    1: dataAspectRatio = dhorizontal/dvertical
    2: dataAspectRatio = dvertical/dhorizontal
  endcase

  ; ---- initialize the device --------------------------------------------------

  if (n_elements(info) NE 0) then $
    widget_control, info.status, set_value = 'status: initialize device'


  if output.type EQ 1 then begin
    current_device = !d.name
    set_plot, 'PS'

    if dataAspectRatio GE 1. then begin
      iorient = 0
    endif else begin
      iorient = 1
    endelse

    case iorient of

    ; portrait orientation
    0: begin
      xsize = 7.5
      ysize = 9.5

      device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
        XOFF = 0.5, YOFF = 0.75, /INCH, /COLOR, $
        BITS_PER_PIXEL = 8
    end

    ; landscape orientation
    1: begin
      xsize = 10.
      ysize = 6.

      device, FILE = outfile, XSIZE = xsize, YSIZE = ysize, $
        XOFF = 1.25, YOFF = 10.5, /INCH, /COLOR, $
        BITS_PER_PIXEL = 8, /LANDSCAPE
      end
    endcase

    deviceAspect = ysize/xsize
        
  endif else begin
    deviceAspect = output.vsize/output.hsize
  endelse

  aspect = dataAspectRatio/deviceAspect

  ; ---- determine the bounds for the plot --------------------------------------

  ; set the normal coordinates of the portion of the display/page you
  ; wish to use -- leave a little margin so the plots don't run to the
  ; edge of the page
  page_nx1 = 0.05
  page_nx2 = 0.95

  page_ny1 = 0.05
  page_ny2 = 0.90

  ; if we are considering multiple plots / page, we'd set some
  ; boundaries here -- just take the full page for now

  dpagex = page_nx2 - page_nx1
  dpagey = page_ny2 - page_ny1

  nx1 = page_nx1
  nx2 = page_nx2

  ny1 = page_ny1
  ny2 = page_ny2

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
    dpx = dpy/aspect < .70* dpagex

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
    dx_legend = dpagex + page_nx1 - px2
    lwidth = dx_legend/4 < 0.25*(px2 - px1)
    lcenter = 0.5* (dpagex + page_nx1 + px2)
    legend_pos = [lcenter-0.5*lwidth, py1, lcenter+0.5*lwidth, py2]

  endif else begin 

    ; the x size of the data sets the scaling of the plot

    ; set the initial x min and max normal coordiantes
    px1 = nx1
    px2 = nx2

    dpx = px2 - px1
    dpy = aspect*dpx < .7* dpagey

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
    dy_legend = py1 - (page_ny2 - dpagey)
    lheight = dy_legend/4 < 0.25*(py2 - py1)
    lcenter = 0.5*(py1 + (page_ny2 - dpagey))
    legend_pos = [px1, lcenter-0.5*lheight, px2, lcenter+0.5*lheight]
  endelse


  ; create the plot window if necessary
  if (ifile EQ fileInfo.startSuffix) then begin
    if (output.type NE 1) then $
      window, XSIZE = output.hsize, YSIZE = output.vsize, $
      TITLE = 'AMR plot'
  endif

  ;-----------------------------------------------------------------------------
  ; create the axes -- just get the data coordinates defined
  ;-----------------------------------------------------------------------------

  ; create the title -- if we are printing it
  if (options.annotate EQ 1) then begin
    title = titleName
  endif else begin
    title = ' '
  endelse

  case islice_dir of
    0: begin  ; x-y
      case problem.orientation of
        0: begin
          plot, [xmin, xmax], [ymin, ymax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [xmin, xmax], [ymax, ymin]
        end

        1: begin
          plot, [ymin, ymax], [xmin, xmax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [ymin, ymax], [xmax, xmin]
        end

        2: begin
          plot, [xmin, xmax], [ymin, ymax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, YRANGE=[ymax,ymin], TITLE = title
          oplot, [xmin, xmax], [ymax, ymin]
        end
      endcase
    end

    1: begin  ; x-z
      case problem.orientation of
        0: begin
          plot, [xmin, xmax], [zmin, zmax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [xmin, xmax], [zmax, zmin]
        end

        1: begin
          plot, [zmin, zmax], [xmin, xmax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [zmin, zmax], [xmax, xmin]
        end

        2: begin
          plot, [xmin, xmax], [zmin, zmax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [xmin, xmax], [zmax, zmin]
        end
      endcase
    end

    2: begin  ; y-z
      case problem.orientation of
        0: begin
          plot, [ymin, ymax], [zmin, zmax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [ymin, ymax], [zmax, zmin]
        end

        1: begin
          plot, [zmin, zmax], [ymin, ymax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [zmin, zmax], [ymax, ymin]
        end

        2: begin
          plot, [ymin, ymax], [zmin, zmax], POS = plot_pos, $
            XSTYLE = 5, YSTYLE = 5, TITLE = title
          oplot, [ymin, ymax], [zmax, zmin]
        end
      endcase
    end
  endcase

  if options.procdist eq '' then begin

    ;------------------------------------------------------------------------------
    ; determine whether to sub-sample the data
    ;------------------------------------------------------------------------------

    ; assume 1024 pixels on screen and 4096 on paper

    sample = 0

    max_refine = max(tree[*].lrefine)

    ; compute the number of pixels in a uniform grid of our sub-domain

    case islice_dir of 
      0: max_pixels = params.ntopx*params.nxb*2.^(max_refine-1)*hfrac > $
        params.ntopy*params.nyb*2.^(max_refine-1)*vfrac

      1: max_pixels = params.ntopx*params.nxb*2.^(max_refine-1)*hfrac > $
        params.ntopz*params.nzb*2.^(max_refine-1)*vfrac

      2: max_pixels = params.ntopy*params.nyb*2.^(max_refine-1)*hfrac > $
        params.ntopz*params.nzb*2.^(max_refine-1)*vfrac
    endcase

    if output.type EQ 1 then begin
      ideal = 4096

      if max_pixels GT ideal then $
        sample = fix(alog(max_pixels/ideal) / alog(2.)) > 0

    endif else begin
      ideal = 1024

      if max_pixels GT ideal then $
        sample = fix(alog(max_pixels/ideal) / alog(2.)) > 0

    endelse

    ; make sure we have atleast one point in every block
    if (2^sample GT (params.nxb < params.nyb < params.nzb)) then begin
      sample = fix(alog(params.nxb < params.nyb < params.nzb)/alog(2.))
    endif


    ;------------------------------------------------------------------------------
    ; scale the data
    ;------------------------------------------------------------------------------
    if (n_elements(info) NE 0) then $
      widget_control, info.status, SET_VALUE = 'status: scaling data'

    ; first put it on a uniform grid
    if Keyword_Set(debug) then begin 
      print, 'ranges:'
      print, xmin, xmax
      print, ymin, ymax
      print, zmin, zmax
     endif   

    ; if we are max'ing the variable, treat it specially
    if options.max EQ 0 then begin

      ; merge it -- corners are handled by merge_amr automagically
      nsdata = merge_amr(temp_arr, TREE=tree, PARAMETERS=params, $
                         XMERGE=x, YMERGE=y, ZMERGE=z, $
                         SAMPLE = sample, $
                         XRANGE = [xmin,xmax], $
                         YRANGE = [ymin,ymax], $
                         ZRANGE = [zmin,zmax])

      nsdata = reform(temporary(nsdata))

      if options.log EQ 1 then nsdata = alog10(temporary(nsdata))

      data_ptr = nsdata
      xmerge = x
      ymerge = y
      zmerge = z

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
    endif else begin

      ; merge it manually
      stempvar = merge_amr(temp_arr, TREE=tree, PARAMETERS=params, $
                           XMERGE = x, YMERGE = y, ZMERGE = z, $
                           SAMPLE = sample, $
                           XRANGE = [xmin,xmax], $
                           YRANGE = [ymin,ymax], $
                           ZRANGE = [zmin,zmax])

      stempvar = reform(temporary(stempvar))

      if options.log EQ 1 then stempvar = alog10(temporary(stempvar))

      ; take the maximum
      if ifile EQ fileInfo.startSuffix then begin
        max_var = stempvar
        if Keyword_Set(debug) then print, 'max = ', max(max_var), min(max_var)
      endif else begin
        max_var = stempvar > max_var
      endelse

      data_ptr = max_var
      xmerge = x
      ymerge = y
      zmerge = z

      undefine, stempvar

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

        
    ;------------------------------------------------------------------------------
    ; plot the data
    ;------------------------------------------------------------------------------
    if (n_elements(info) NE 0) then widget_control, info.status, $
      set_value = 'status: plotting data'

    ; convert the limits of the domain into normal coordinates
    case islice_dir of
      0: begin
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
      end

      1: begin
        case problem.orientation of
          0: begin
            lower = convert_coord([xmin, ymin], /DATA, /TO_NORMAL)
            upper = convert_coord([xmax, ymax], /DATA, /TO_NORMAL)
          end
          1: begin
            lower = convert_coord([zmin, xmin], /DATA, /TO_NORMAL)
            upper = convert_coord([zmax, xmax], /DATA, /TO_NORMAL)
          end
          2: begin
            lower = convert_coord([xmin, ymax], /DATA, /TO_NORMAL)
            upper = convert_coord([xmax, ymin], /DATA, /TO_NORMAL)
            nsdata = rotate(temporary(nsdata),7)
          end
        endcase
      end

      2: begin
        case problem.orientation of
          0: begin
            lower = convert_coord([ymin, zmin], /DATA, /TO_NORMAL)
            upper = convert_coord([ymax, zmax], /DATA, /TO_NORMAL)
          end
          1: begin
            lower = convert_coord([zmin, ymin], /DATA, /TO_NORMAL)
            upper = convert_coord([zmax, ymax], /DATA, /TO_NORMAL)
          end
          2: begin
            lower = convert_coord([ymin, zmax], /DATA, /TO_NORMAL)
            upper = convert_coord([ymax, zmin], /DATA, /TO_NORMAL)
            nsdata = rotate(temporary(nsdata),7)
          end
        endcase
      end
    endcase

    sz = size(nsdata)

    tvimage, nsdata, pos = [lower[0], lower[1], upper[0], upper[1]], /overplot

    if (params.geometry EQ "CARTESIAN") then begin
      if options.blocks then begin
        draw3d_blocks, color('black'), TREE = tree, PARAMETERS = params, $
          ORIENTATION = problem.orientation, $
          XRANGE = [xmin,xmax], $
          YRANGE = [ymin,ymax], $
          ZRANGE = [zmin,zmax], $
          SLICE_DIR = islice_dir
      endif
    endif
  endif else begin 
    ; just show the distribution of blocks on processors
    if (params.geometry EQ "CARTESIAN") then begin
      print, 'procdist show!'
      show_3d_proc_distribution, LEVEL=options.procdist, TREE = tree, PARAMETERS = params, $
        ORIENTATION = problem.orientation, $
        XRANGE = [xmin,xmax], $
        YRANGE = [ymin,ymax], $
        ZRANGE = [zmin,zmax], $
        SLICE_DIR = islice_dir, $
        SHOW_NUMBERS=options.showprocnums

      draw3d_blocks, color('black'), LEVEL=options.procdist, $
        TREE = tree, PARAMETERS = params, $
        ORIENTATION = problem.orientation, $
        XRANGE = [xmin,xmax], $
        YRANGE = [ymin,ymax], $
        ZRANGE = [zmin,zmax], $
        SLICE_DIR = islice_dir
    endif
  endelse 

  ;------------------
  ; restore the axes
  ;------------------
  case islice_dir of 
    0: begin
      case problem.orientation of 
        0: begin
          if (problem.name EQ 'xray_cyl') then begin
            case options.tick of 
              0: begin ; no tick marks
              end
              1: begin ; display tick marks
                axis, xmin, ymin, xaxis = 0, XSTYLE = 1, $
                  xtitle = 'r (cm)'
                axis, xmin, ymin, yaxis = 0, YSTYLE = 1, $
                  ytitle = 'z (cm)'
                axis, xmin, ymax, xaxis = 1, XSTYLE = 1, $
                  xtickformat = 'nolabel'
                axis, xmax, ymin, yaxis = 1, YSTYLE = 1, $
                  ytickformat = 'nolabel'
              end
            endcase
          endif else begin
            case options.tick of 
              0: begin ; no tick marks
              end
              1: begin ; display tick marks
                axis, xmin, ymin, xaxis = 0, XSTYLE = 1, $
                  xtitle = 'x (cm)'
                axis, xmin, ymin, yaxis = 0, YSTYLE = 1, $
                  ytitle = 'y (cm)'
                axis, xmin, ymax, xaxis = 1, XSTYLE = 1, $
                  xtickformat = 'nolabel'
                axis, xmax, ymin, yaxis = 1, YSTYLE = 1, $
                  ytickformat = 'nolabel'
              end
            endcase
          endelse
        end

        1: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, ymin, xmin, xaxis = 0, XSTYLE = 1, xtitle = 'y (cm)'
              axis, ymin, xmin, yaxis = 0, YSTYLE = 1, ytitle = 'x (cm)'
              axis, ymin, xmax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, ymax, xmin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
        end

        2: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, xmin, ymax, xaxis = 0, XSTYLE = 1, xtitle = 'x (cm)'
              axis, xmin, ymax, yaxis = 0, YSTYLE = 1, ytitle = 'y (cm)'
              axis, xmin, ymin, xaxis = 1, XSTYLE = 1, $
                xtickformat = 'nolabel'
              axis, xmax, ymax, yaxis = 1, YSTYLE = 1, $
                ytickformat = 'nolabel'
            end
          endcase
        end
      endcase
    end

    1: begin
      case problem.orientation of 
        0: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, xmin, zmin, xaxis = 0, XSTYLE = 1, xtitle = 'x (cm)'
              axis, xmin, zmin, yaxis = 0, YSTYLE = 1, ytitle = 'z (cm)'
              axis, xmin, zmax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, xmax, zmin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
        end
        1: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, zmin, xmin, xaxis = 0, XSTYLE = 1, xtitle = 'z (cm)'
              axis, zmin, xmin, yaxis = 0, YSTYLE = 1, ytitle = 'x (cm)'
              axis, zmin, xmax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, zmax, xmin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
        end
        2: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, xmin, zmin, xaxis = 0, XSTYLE = 1, xtitle = 'x (cm)'
              axis, xmin, zmin, yaxis = 0, YSTYLE = 1, ytitle = 'z (cm)'
              axis, xmin, zmax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, xmax, zmin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
        end
      endcase
    end

    2: begin
      case problem.orientation of 
        0: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, ymin, zmin, xaxis = 0, XSTYLE = 1, xtitle = 'y (cm)'
              axis, ymin, zmin, yaxis = 0, YSTYLE = 1, ytitle = 'z (cm)'
              axis, ymin, zmax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, ymax, zmin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
        end
        1: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, zmin, ymin, xaxis = 0, XSTYLE = 1, xtitle = 'z (cm)'
              axis, zmin, ymin, yaxis = 0, YSTYLE = 1, ytitle = 'y (cm)'
              axis, zmin, ymax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, zmax, ymin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
        end
        2: begin
          case options.tick of 
            0: begin ; no tick marks
            end
            1: begin ; display tick marks
              axis, ymin, zmin, xaxis = 0, XSTYLE = 1, xtitle = 'y (cm)'
              axis, ymin, zmin, yaxis = 0, YSTYLE = 1, ytitle = 'z (cm)'
              axis, ymin, zmax, xaxis = 1, XSTYLE = 1, xtickformat = 'nolabel'
              axis, ymax, zmin, yaxis = 1, YSTYLE = 1, ytickformat = 'nolabel'
            end
          endcase
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
    case islice_dir of
      0: begin
        case problem.orientation of
          0: begin
            xyouts, label_info.posx*(xmax-xmin)+xmin, $
              label_info.posy*(ymax-ymin)+ymin, label_info.label, $
              alignment = .5, charsize = label_info.size, $
                                charthick = label_info.thick, color = label_info.color
          end
          1: begin
            xyouts, label_info.posx*(ymax-ymin)+ymin, $
              label_info.posy*(xmax-xmin)+xmin, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color=label_info.color
          end
          2: begin
            xyouts, label_info.posx*(xmax-xmin)+xmin, $
              label_info.posy*(ymin-ymax)+ymax, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color=label_info.color
          end
        endcase
      end
      1: begin
        case problem.orientation of
          0: begin
            xyouts, label_info.posx*(xmax-xmin)+xmin, $
              label_info.posy*(zmax-zmin)+zmin, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color=label_info.color
          end
          1: begin
            xyouts, label_info.posx*(zmax-zmin)+zmin, $
              label_info.posy*(xmax-xmin)+xmin, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color=label_info.color
          end
          2: begin
            xyouts, label_info.posx*(xmax-xmin)+xmin, $
              label_info.posy*(zmin-zmax)+zmax, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color=label_info.color
          end
        endcase
      end
      2: begin
        case problem.orientation of
          0: begin
            xyouts, label_info.posx*(ymax-ymin)+ymin, $
              label_info.posy*(zmax-zmin)+zmin, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color = label_info.color
          end
          1: begin
            xyouts, label_info.posx*(zmax-zmin)+zmin, $
              label_info.posy*(ymax-ymin)+ymin, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color=label_info.color
          end
          2: begin
            xyouts, label_info.posx*(ymax-ymin)+ymin, $
              label_info.posy*(zmin-zmax)+zmax, label_info.label, $
              alignment = .5, charsize = label_info.size, $
              charthick = label_info.thick, color = label_info.color
          end
        endcase
      end
    endcase
  endif


  ;-----------------------
  ; add a colorbar legend if desired
  ;-----------------------
  if options.colorbar EQ 1 then begin
    if aspect LT 1 then begin
      case options.log of
        0: colorbar2, float(variable.min[0]), float(variable.max[0]), legend_pos, $
          COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

        1: colorbar2, alog10(variable.min[0]), alog10(variable.max[0]), $
          legend_pos, COLORMIN = colorMin, COLORMAX = colorMax, $
          CHARSIZE = 1.5
      endcase

    endif else begin
      case options.log of
        0: vcolorbar, float(variable.min[0]), float(variable.max[0]), legend_pos, $
          COLORMIN = colorMin, COLORMAX = colorMax, CHARSIZE = 1.5

        1: vcolorbar, alog10(variable.min[0]), alog10(variable.max[0]), $
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

  gridtxt = 'number of blocks = ' + string(format = '(i7)', params.totBlocks) 
  gridtxt2 = 'AMR levels = ' + string(format = '(i5)', max(tree[*].lrefine))

  ; if desired, print the time, filename, and credits
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

  ; empty the buffer
  empty

  case output.type of
    1: begin
      device, /close
      set_plot, current_device
    end
    2: color_bitmap, outfile
    else:
  endcase

  if (n_elements(info) NE 0) then widget_control, info.status, $
    set_value = 'status: awaiting orders . . .'
    
endfor

end
