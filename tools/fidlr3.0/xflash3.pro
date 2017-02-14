;=============================================================================
; xflash  -- general FLASH data plotting interface
;
;   plot 1, 2, or 3-d data with velocity vectors, contours, and
;   particle data overlayed.
;
; type xflash3 at the idl prompt to run
; 
; modify xflash_defaults.pro to add additional problems
;
;   xflash is an interface to the collection of idl routines which
;   work with flash data.  xflash collects the plot options and passes
;   them on to xplot_amr, which does the actual plotting.  
;
;=============================================================================

; event handler
pro xflash_event, ev

; common all the controls into the widget

common variables, varnames, native_num

common save_parts, particle_widget
common save_label, label_info

common save, filename, pblm_ctr_lim, orientation, ndim, geometry

common save_color_info, current_color_id, current_color_name
common save_default_info, current_default_id, current_default_name
common save_contour_info, contours
common save_vector_info, vector
common save_histogram_info, nbins, max_scalep
common save_number_info, current_number_id, current_number_name

common save_ptr, tree_ptr, data_ptr, params_ptr, xmerge_ptr, ymerge_ptr, zmerge_ptr

common save_path, current_path
common save_filename, filename_old, filename_last_read, filename_last_read_id, $
  old_start_sfx, old_end_sfx

common save_slice_plane, last_slice_plane

common save_filetype, filetype

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; convert the debug flag into a plain parameters
debug = info.debug

; take action pased on which widget was changed
case uval of
    
  ;-------------------------------------------------------------------------
  ; file actions
  ;-------------------------------------------------------------------------
  'fopen': begin
    filters = ['*hdf*', '*_ncmpi_*']
    filename = dialog_pickfile(FILTER=filters, /MUST_EXIST, $
                               PATH=current_path, DIALOG_PARENT=info.mainBase)

    ; if no file was selected then skip every thing else
    if (filename NE "") then begin

      if (strpos(filename, '.dat') < 0) then begin
        ; file is not a .dat file

        ; update the current path for next time
        pathEnd = strpos(filename, '/', /REVERSE_SEARCH)
        current_path = strcompress(strmid(filename, 0, pathEnd))

        baseNmEnd = strpos(filename, '_', /REVERSE_SEARCH)
        filenameBase = strcompress(strmid(filename, 0, baseNmEnd+1))
        print, 'filenameBase is: ', filenameBase
                
        oldBaseNmEnd = strpos(filename_old, '_', /REVERSE_SEARCH)
        oldFilenameBase = strcompress(strmid(filename_old, 0, oldBaseNmEnd+1))
        print, 'oldFilenameBase is: ', oldFilenameBase

        same_base = 0
        if (oldFilenameBase EQ filenameBase) then same_base = 1

        itype = determine_file_type(filename) ; 'itype'= 1 if hdf5, 2 if netcdf, or -1 if unrecognized

        if (itype EQ -1) then begin
          result = dialog_message(['ERROR: file does not appear to be', $
                                   'a valid FLASH data file.'], $
                                   title='ERROR', $
                                   dialog_parent=info.mainBase)
        endif else begin
          widget_control, info.finfo, sensitive=1
                    
          pathEnd = strpos(filename, '/', /reverse_search)

          ; ---- update the prototype file label ----------------------------------------
          widget_control, info.fprototype, $
            set_value='Prototype File: '+ $
            strcompress(strmid(filename,pathEnd+1, strlen(filename)-pathEnd))

          ; ---- get the dimensionality -------------------------------------------------
          ndim = determine_file_dimensionality(filename)
          print, 'ndim is: ', ndim

          ; ---- setup some pointers to pass the data -----------------------------------
          ; if (n_elements(tree_ptr) GT 0) then ptr_free, tree_ptr
          ; if (n_elements(params_ptr) GT 0) then ptr_free, params_ptr
          ; if (n_elements(data_ptr) GT 0) then ptr_free, data_ptr

          ; is the file double (2) or single (1) precision?
          precision = determine_file_precision(filename)
          print, 'precision is: ', precision

          ; call a standard function to set up the tree 
          ; note that original flash call was
          ;   gid:lonarr(2*ndim+1+2^ndim), hence definitions below
          nfaces = 2*ndim
          nchild = 2^ndim
          define_tree,precision,nfaces,nchild,ndim,TREE=tree
          tree_ptr = ptr_new(tree)
          params_ptr = ptr_new(params)
                    
          data_ptr = ptr_new(0.0)
          xmerge_ptr = ptr_new(0.0)
          ymerge_ptr = ptr_new(0.0)
          zmerge_ptr = ptr_new(0.0)

          ; ---- get the suffix ---------------------------------------------------------
          baseNmEnd = strpos(filename, '_', /REVERSE_SEARCH)
                    
          suffix = strcompress(strmid(filename,baseNmEnd+1,4))

          ; ---- sensitize the different plot options -----------------------------------
          widget_control, info.fstartSuffix, sensitive=1
          widget_control, info.fstartSuffix, SET_VALUE=suffix

          ; set end suffix to beginning suffix by default
          widget_control, info.fendSuffix, sensitive=1
          widget_control, info.fendSuffix, SET_VALUE=suffix

          widget_control, info.fstep, sensitive=1

          widget_control, info.out_buttons, sensitive=1
          widget_control, info.hsize, sensitive=1
          widget_control, info.vsize, sensitive=1

          ; make the option buttons clickable...
          widget_control, info.opt_buttons, sensitive=1

          widget_control, info.zoomLabel, sensitive=1

          ; I decided a better default was to start with the 'auto'
          ; button turned on and these windows turned off -nttaylor
          ; widget_control, info.varMin, sensitive=1
          ; widget_control, info.varMax, sensitive=1

          widget_control, info.autoData, sensitive=1

          ; ---- populate the variable drop box -----------------------------------------
          ; get the names of the unknowns stored in the file
          varnames = get_var_list(filename)
          print, 'varnames is: ', varnames

          if n_elements(varnames) EQ 1 then begin
                        
            if (varnames[0] EQ -1) then begin

              ; files from FLASH1 will not have the unknown names stored
              result = dialog_message(['Warning: file does not have' ,  $
                                       'variable name information.', $
                                       'Unable to proceed . . .'], $
                                       title = 'Warning', $
                                       dialog_parent = info.mainBase)
            endif

          endif

          if (same_base EQ 0) then begin

            ; ---- check for particles ----------------------------------------------------
            numParticles = get_particle_number(filename)

            ; populate parallel block distribution level selection dropdown list
            get_lrefine_max_min, filename, ndim, MAX=maxLevels, MIN=minLevels
            levels = indgen(maxLevels - minLevels+1) + 1
            levels = string(levels)
            print, 'maxLevels is: ', maxLevels
            print, 'minLevels is: ', minLevels
            print, 'levels is: ', levels

            ; create structure to hold data from
            ; particle widget (cf. xparticle.pro)
            particle_widget = {enabled:0, $
                               plot_vel:0, $
                               sym_size: 1.0, $
                               show_tag: 0, $
                               typical_velocity:10., $
                               data_enabled:0, $
                               dummy:0, $
                               traj:0,$
                               traj_color:0,$
                               color_max:1.e10, $
                               color_min: 1.e7, $
                               log:1,$
                               partnum:1880}

            widget_control, info.levelDroplist, /destroy
            widget_control, info.showprocnums, /destroy

            info.levelDroplist = widget_droplist(info.blockDistBase, $
                                                 title = 'Level: ', $
                                                 uvalue = 'blockDistLevel', $ 
                                                 value = levels) 

            info.showProcNums = cw_bgroup(info.blockDistBase, ['Show Proc Numbers'], $
                                          /nonexclusive, $
                                          set_value=[0], $
                                          uvalue = 'showProcNums')

            widget_control, info.levelDroplist, sensitive=0
            widget_control, info.showProcNums, sensitive=0


            ; ---- determine which derived variables are possible -------------------------
            native_num = (size(varnames))[1]

            ; total velocity
            totVelDefined = 0

            case ndim of
              1: begin
                if (var_index('velx') NE -1) then begin
                  varnames = add_var(temporary(varnames), 'tot_vel')
                  totVelDefined = 1
                endif
              end
                            
              2: begin
                if (var_index('velx') NE -1 AND $
                  var_index('vely') NE -1) then begin
                  varnames = add_var(temporary(varnames), 'tot_vel')
                  totVelDefined = 1
                endif
              end
                            
              3: begin
                if (var_index('velx') NE -1 AND $
                  var_index('vely') NE -1 AND $
                  var_index('velz') NE -1) then begin
                  varnames = add_var(temporary(varnames), 'tot_vel')
                  totVelDefined = 1
                endif
              end
            endcase

                        
            ; sound speed
            if var_index('pres') NE -1 AND $
               var_index('dens') NE -1 AND $
               var_index('gamc') NE -1 then $
               varnames = add_var(temporary(varnames), 'snd_spd')
                        
            ; mach number
            if (var_index('pres') NE -1 AND $
               var_index('dens') NE -1 AND $
               var_index('gamc') NE -1 AND $
               totVelDefined EQ 1) then $
               varnames = add_var(temporary(varnames), 'mach')

            ; internal energy and ratio of kinetic to internal
            if (var_index('ener') NE -1 AND $
               totVelDefined EQ 1) then begin
                 varnames = add_var(temporary(varnames), 'int_ener')
                 varnames = add_var(temporary(varnames), 'ekin/eint')
            endif

            ; gravitational energy / internal energy
            if (var_index('gpot') NE -1 AND $
                var_index('ener') NE -1 AND $
                var_index('velx') NE -1 AND $
                var_index('vely') NE -1) then begin
                  varnames = add_var(temporary(varnames), 'egrav/eint')

            endif

            widget_control, info.vars, /destroy
                        
            ; make sure to store the new instance of the droplist back to the info
            ; structure so it is available to others
            info.vars = widget_droplist(info.varBase, $
                                        title = 'Mesh Variables: ', $
                                        uvalue = 'var', value = varnames)


            ; ---- load the generic defaults ----------------------------------------------
            xflash_defaults, 'Generic', $
                             CTR_LIM=ctr_lim, ISWP=iswp, $
                             CONTOURS=contours_def, VECTOR=vector_def, $
                             NBINS = nbins, HIST_SCALE = max_scale

            ; save the contours
            contours = contours_def
            vector = vector_def

            ; we need to save ctr_lim, so other events can access it
            pblm_ctr_lim = ctr_lim
            orientation = iswp


            ; get the geometry
            geometry = determine_geometry(filename)

            if (geometry eq "") then begin
               continue = dialog_message("Geometry string is blank in your file.  You can pick a geometry or abort.  Do you want to pick a geometry?", /question, dialog_parent=info.mainBase)
               if (continue eq "Yes") then begin
                 geometry = xpick_geometry(info.mainBase)
               endif else begin
                 widget_control, info.mainBase, /destroy
               endelse 
            endif 


            ; ---- activate the options based on the # of dims ----------------------------
            widget_control, info.zoomLabel, sensitive = 1
                        
            widget_control, info.xmin, sensitive = 1
            widget_control, info.xmax, sensitive = 1

            case ndim of
              1: begin

                widget_control, info.ymin, sensitive = 0
                widget_control, info.ymax, sensitive = 0

                widget_control, info.zmin, sensitive = 0
                widget_control, info.zmax, sensitive = 0

                widget_control, info.zoombox, sensitive = 0
                widget_control, info.reset, sensitive = 0

                widget_control, info.slicePlane, sensitive = 0
                                
                widget_control, info.contourOpts, sensitive = 0
                widget_control, info.vectorOpts, sensitive = 0
                widget_control, info.particleOpts, sensitive = 0

                ; DEV histogram is not working disable buttons throughout
                widget_control, info.histogramOpts, sensitive = 0
                widget_control, info.labelOpts, sensitive = 1

                print, 'numParticles is: ', numParticles
                if (numParticles GT 0) then begin
                  ; turn on the 'particle options' button
                  widget_control, info.particleOpts, sensitive = 1
                  particle_vars = get_particle_vars(filename)
                endif

              end

              2: begin

                widget_control, info.ymin, sensitive = 1
                widget_control, info.ymax, sensitive = 1

                widget_control, info.zmin, sensitive = 0
                widget_control, info.zmax, sensitive = 0

                widget_control, info.enableBlockDist, sensitive=1
                                
                widget_control, info.zoombox, sensitive = 1
                widget_control, info.reset, sensitive = 1

                widget_control, info.slicePlane, sensitive = 0
                                
                widget_control, info.contourOpts, sensitive = 1
                widget_control, info.vectorOpts, sensitive = 1

                ; DEV histogram is not working, disable button
                widget_control, info.histogramOpts, sensitive = 0
                widget_control, info.labelOpts, sensitive = 1
                                
                print, 'numParticles is: ', numParticles
                if (numParticles GT 0) then begin
                  ; turn on the 'particle options' button
                  widget_control, info.particleOpts, sensitive = 1
                  particle_vars = get_particle_vars(filename)
                endif
              end

              3: begin
                                
                widget_control, info.ymin, sensitive = 1
                widget_control, info.ymax, sensitive = 1
                                
                ; by default, we are in the x-y plane, so only one z value is needed
                widget_control, info.zmin, sensitive = 1
                widget_control, info.zmax, sensitive = 0

                widget_control, info.enableBlockDist, sensitive=1
                                
                widget_control, info.zoombox, sensitive = 0
                widget_control, info.reset, sensitive = 0

                widget_control, info.slicePlane, sensitive = 1
                last_slice_plane = -1
                widget_control, info.contourOpts, sensitive = 0
                widget_control, info.vectorOpts, sensitive = 0
                widget_control, info.particleOpts, sensitive = 0
                ; DEV histogram is not working, disable button
                widget_control, info.histogramOpts, sensitive = 0
                widget_control, info.labelOpts, sensitive = 1
                                
              end
                            
            endcase

          endif ; (same_base EQ 0)

          ; ---- floating label info ----------------------------------------------------
          label = strarr(9)
          label[0] = '(a)'

          label_info = {enabled:0, $
                        label:label, $
                        size:1.5, $
                        thick:1.,$
                        multi:0, $
                        multi_opt:0, $
                        color:color('black'), $
                        posx:.9, $
                        posy:.92 }


          ; ---- activate the plot button -----------------------------------------------
          widget_control, info.plot, sensitive = 1
          ; DEV histogram is not working, disable button
          widget_control, info.histogram, sensitive = 0

          ; ---- store the new value of info.vars ---------------------------------------
          widget_control, ev.top, set_uvalue = info, /no_copy

        endelse

        filetype = 1  ; hdf5 and netcdf files are filetype 1
        filename_old = filename

      endif else begin ; (strpos(filename, '.dat') is greater than 0)

        ; we have a .dat file
        pathEnd = strpos(filename, '/', /REVERSE_SEARCH)

        widget_control, info.fprototype, $
                        set_value='Prototype File: '+ $
                        strcompress(strmid(filename,pathEnd+1, $
                                    strlen(filename)-pathEnd))

        varnames = get_dat_vars(filename)

        widget_control, info.vars, /destroy
                
        ; make sure to store the new instance of the droplist back to the info
        ; structure so it is available to others
        info.vars = widget_droplist(info.varBase, $
                                    title = 'Variables: ', $
                                    uvalue = 'var', value = varnames)
                
        ; activate the necessary widgets
        widget_control, info.fstartSuffix, sensitive=0
        widget_control, info.fendSuffix, sensitive=0
        widget_control, info.fstep, sensitive=0
                
        widget_control, info.out_buttons, sensitive=1
        widget_control, info.hsize, sensitive=1
        widget_control, info.vsize, sensitive=1
                
        widget_control, info.opt_buttons, sensitive=1
                
        widget_control, info.zoomLabel, sensitive=0
                
        widget_control, info.varMin, sensitive=1
        widget_control, info.varMax, sensitive=1
                
        widget_control, info.autoData, sensitive=1

        widget_control, info.xmin, sensitive = 0
        widget_control, info.xmax, sensitive = 0

        widget_control, info.ymin, sensitive = 0
        widget_control, info.ymax, sensitive = 0
                
        widget_control, info.zmin, sensitive = 0
        widget_control, info.zmax, sensitive = 0
                
        widget_control, info.zoombox, sensitive = 0
        widget_control, info.reset, sensitive = 0

        widget_control, info.slicePlane, sensitive = 0

        widget_control, info.contourOpts, sensitive = 0
        widget_control, info.vectorOpts, sensitive = 0
        widget_control, info.histogramOpts, sensitive = 0
                
        widget_control, info.plot, sensitive = 1


        xflash_defaults, 'Generic', CTR_LIM=ctr_lim                        
                
        ; we need to save ctr_lim, so other events can access it
        pblm_ctr_lim = ctr_lim


        ; ---- store the new value of info.vars ---------------------------------------
        widget_control, ev.top, set_uvalue = info, /no_copy


        ; dat files are filetype 2
        filetype = 2
      endelse

    endif else begin
      filename = filename_old
    endelse

  end ; fopen


  'finfo': begin
    xfile_info, filename
  end


  ;-----------------------------------------------------------------------------
  ; default options
  ;-----------------------------------------------------------------------------
  'def_item': begin

    widget_control, ev.id, Get_Value=defaultProblemName

    ; to fake the appearance of a checkmark next to the selected menu
    ; item, the first two characters of the problem name are spaces, or 
    ; an x followed by a space.  Strip these off when calling
    ; xflash_defaults

    ; load the defaults for the selected problem
    xflash_defaults, strmid(defaultProblemName,2), $
                     CTR_LIM=ctr_lim, DEFAULT_NAMES=problems, ISWP=iswp, $
                     CONTOURS=contours_def, VECTOR=vector_def

    contours = contours_def
    vector = vector_def

    ; remove the x off the old problem name and put it on the new one
    widget_control, current_default_id, GET_VALUE=oldName
    widget_control, current_default_id, SET_VALUE='  ' + strmid(oldName,2)

    widget_control, ev.id, SET_VALUE='x ' + strmid(defaultProblemName,2)

    ; store the index         
    current_default_id = ev.id
    current_default_name = strmid(defaultProblemName,2)

    ; store the data ranges globally
    pblm_ctr_lim = ctr_lim
    orientation = iswp

    ; get the variable number
    ivar = widget_info(info.vars, /droplist_select)

    ; update the data ranges
    widget_control, info.varMin, set_value = pblm_ctr_lim[ivar,0]
    widget_control, info.varMax, set_value = pblm_ctr_lim[ivar,1]

    ; reset the zoom options
    widget_control, info.xmin, set_value = -1
    widget_control, info.xmax, set_value = -1

    widget_control, info.ymin, set_value = -1        
    widget_control, info.ymax, set_value = -1

    widget_control, info.zmin, set_value = -1
    widget_control, info.zmax, set_value = -1

  end ; 'def_item'


  ;-----------------------------------------------------------------------------
  ; color options
  ;-----------------------------------------------------------------------------
  'colorItem': begin

    widget_control, ev.id, GET_VALUE=colorName

    ; I could not figure out how to have IDL put a checkmark next to the
    ; color name that is currently selected off the menu, so we fake this
    ; below

    ; remove the 'x' off the old color
    widget_control, current_color_id, GET_VALUE=oldName
    widget_control, current_color_id, SET_VALUE='  ' + strmid(oldName,2)

    ; put an 'x' on the new one to indicate that it is selected
    widget_control, ev.id, SET_VALUE = 'x ' + strmid(colorName,2)

    ; store the new id and name so we can use them next time
    current_color_id = ev.id
    current_color_name = strmid(colorName,2)

  end ; 'colorItem'

  ;-----------------------------------------------------------------------------
  ; xy count options
  ;-----------------------------------------------------------------------------
  'xycountItem': begin

    widget_control, ev.id, GET_VALUE=numberName

    widget_control, current_number_id, GET_VALUE=oldName
    widget_control, current_number_id, SET_VALUE='  ' + strmid(oldName,2)

    widget_control, ev.id, SET_VALUE = 'x ' + strmid(numberName,2)

    current_number_id = ev.id
    current_number_name = strmid(numberName,2)

  end ; 'xycountItem'

  ;-----------------------------------------------------------------------------
  ; variable changed 
  ;-----------------------------------------------------------------------------
  'var': begin

    ; get the new variable number
    ivar = widget_info(info.vars, /droplist_select)

    ; put the new contour limits in
    widget_control, info.varMin, set_value = pblm_ctr_lim[ivar,0]
    widget_control, info.varMax, set_value = pblm_ctr_lim[ivar,1]

  end ; 'var'

  ;-----------------------------------------------------------------------------
  ; output changed 
  ;-----------------------------------------------------------------------------
  'output': begin

    widget_control, info.out_buttons, get_value = ipost

    if ipost EQ 1 then begin
      widget_control, info.hsize, sensitive = 0
      widget_control, info.vsize, sensitive = 0
    endif else begin
      widget_control, info.hsize, sensitive = 1
      widget_control, info.vsize, sensitive = 1
    endelse

  end

  ;-----------------------------------------------------------------------------
  ; enable the block distribution plot selected
  ;-----------------------------------------------------------------------------
  'enableBlockDist': begin
    widget_control, info.enableBlockDist, get_value = ienable

    if ienable[0] EQ 1 then begin
      widget_control, info.levelDroplist, sensitive=1
      widget_control, info.showProcNums, sensitive=1
    endif else begin
      widget_control, info.levelDroplist, sensitive=0
      widget_control, info.showProcNums, sensitive=0
    endelse
  end

  ;-----------------------------------------------------------------------------
  ; automatic data range limits selected
  ;-----------------------------------------------------------------------------
  'auto': begin

    ; get the value of the auto button
    widget_control, info.autoData, get_value = iauto

    if iauto[0] EQ 1 then begin
      widget_control, info.varMin, sensitive=0
      widget_control, info.varMax, sensitive=0    
    endif else begin
      widget_control, info.varMin, sensitive=1
      widget_control, info.varMax, sensitive=1    
    endelse
  end

  ;-----------------------------------------------------------------------------
  ; slice plane changed (3d)
  ;-----------------------------------------------------------------------------
  'slice': begin
    widget_control, info.slicePlane, get_value = islice_dir

    case islice_dir of
      0: begin
        widget_control, info.xmax, sensitive=1
        widget_control, info.ymax, sensitive=1
        widget_control, info.zmax, sensitive=0
      end

      1: begin
        widget_control, info.xmax, sensitive=1
        widget_control, info.ymax, sensitive=0
        widget_control, info.zmax, sensitive=1
      end

      2: begin
        widget_control, info.xmax, sensitive=0
        widget_control, info.ymax, sensitive=1
        widget_control, info.zmax, sensitive=1
      end
    endcase

    if (islice_dir EQ last_slice_plane) then begin
      widget_control, info.slice1d, SENSITIVE=1
    end else begin
      widget_control, info.slice1d, SENSITIVE=0
    end
  end


  ;-----------------------------------------------------------------------------
  ; contour options pressed
  ;-----------------------------------------------------------------------------
  'ctropt': begin
        
    xcontour, VARIABLES=varnames

  end


  ;-----------------------------------------------------------------------------
  ; vector options pressed
  ;-----------------------------------------------------------------------------
  'vecopt': begin
    xvector, VARIABLES=varnames, NATIVE_NUMBER = native_num
    if Keyword_Set(debug) then print, 'x,y comp = ', vector.xcomp, vector.ycomp
  end

  ;-----------------------------------------------------------------------------
  ; floating label pressed
  ;-----------------------------------------------------------------------------

  'labelopt': begin
    xlabel
  end

  ;-----------------------------------------------------------------------------
  ; histogram options pressed
  ;-----------------------------------------------------------------------------
  'histopt': begin
    xhist
  end


  ;-----------------------------------------------------------------------------
  ; plot pressed
  ;-----------------------------------------------------------------------------
  'plot': begin

    ; file
    if (filetype EQ 1) then begin  ; this is an hdf5 or netcdf file

      widget_control, info.fstartSuffix, get_value = startSuffix
      startSuffix = fix(startSuffix[0])

      widget_control, info.fendSuffix, get_value = endSuffix
      endSuffix = fix(endSuffix[0])

      if (endSuffix LT startSuffix) then endSuffix = startSuffix
            
      widget_control, info.fstep, get_value = step
      step = fix(step[0])
      if (step EQ 0) then step = 1

      fileInfo = {prototype:filename, $
                  startSuffix:startSuffix, $
                  endSuffix:endSuffix, $
                  step:step}


      ; figure out if we need to read the data once again.  If only one
      ; file is specified by the fileInfo structure (i.e. startSuffix =
      ; endSuffix), and it is the same as the last file read. then we
      ; don't need to read the data again.

      baseNmEnd = strpos(fileInfo.prototype, '_', /REVERSE_SEARCH)
      filenameBase = strcompress(strmid(fileInfo.prototype, 0, baseNmEnd+1))
            
      filename_read = filenameBase + string(fileInfo.startSuffix, format = '(i4.4)')

      SPAWN, 'strings ' + filename + ' | head -8 | grep -v \[a-z\] | grep -v \[A-Z\] | tail -1' , filename_read_id

      ; supposedly needed on SuSE9.1 with IDL 6.0; have no way to test it

      filename_read_id = filename_read_id[0]

      ; if we are doing a loop over several files, or if we did a loop last
      ; time, reread the file this time through.
      if (fileInfo.startSuffix EQ fileInfo.endSuffix AND $
        old_start_sfx EQ old_end_sfx) then begin
                
        ; create the filename that we are examining

        if (filename_read EQ filename_last_read) then begin
          if ( filename_read_id EQ filename_last_read_id ) then begin
            iread = 0
          endif else begin
            iread = 1
          endelse
        endif else begin
          iread = 1
        endelse
      endif else begin
        iread = 1
      endelse
            
      filename_last_read = filename_read
      filename_last_read_id = filename_read_id

      old_start_sfx = fileInfo.startSuffix
      old_end_sfx = fileInfo.endSuffix

      ; output
      widget_control, info.out_buttons, get_value = outputType

      widget_control, info.vsize, get_value = vpixels
      widget_control, info.hsize, get_value = hpixels

      ; create a structure to hold the output information
      output = {type:outputType,         $ ; flag for output medium
                hsize:float(hpixels[0]), $ ; # of pixels for graphic in x
                vsize:float(vpixels[0])} ; # of pixels for graphic in y


      ; variable
      ivar = widget_info(info.vars, /droplist_select)

      ; plotting options
      widget_control, info.opt_buttons, get_value = opts

      widget_control, info.enableBlockDist, get_value = ienable

      blockDistLevel = 0
      showProcNums = 0
      if ienable[0] EQ 1 then begin
        get_lrefine_max_min, filename, ndim, MAX=maxLevels, MIN=minLevels
        levels = indgen(maxLevels - minLevels+1) + 1
        levels = string(levels)
        ilevel = widget_info(info.levelDroplist, /droplist_select)
        widget_control, info.showProcNums, get_value = ishowProcNums
        showProcNums = ishowProcNums[0]
        blockDistLevel = levels[long(ilevel[0])]
      endif

      ; data range
      widget_control, info.varMin, get_value = dataMin
      widget_control, info.varMax, get_value = dataMax
            
      widget_control, info.autoData, get_value = auto

      pblm_ctr_lim[ivar,0] = dataMin
      pblm_ctr_lim[ivar,1] = dataMax
            
      variable = {name:varnames(ivar),  $
                  min:dataMin, $
                  max:dataMax, $
                  auto:auto[0]}


      ; create a structure to hold the options
      options = {blocks:opts[3],   $ ; draw block boundaries
                 procdist:blockDistLevel, $ ; color blocks based on processor  
                 showprocnums:showProcNums, $ ; show the processor numbers on the blocks
                 log:opts[0],      $ ; take the log of data
                 annotate:opts[4], $ ; plot time and credits
                 colorbar:opts[5], $ ; plot colorbar
                 max:opts[2],      $ ; plot max of variable
                 abs:opts[1],      $ ; take abs value of var
                 tick:opts[6]} ; display tick marks


      ; zoom options
      widget_control, info.slicePlane, get_value = plane
      last_slice_plane = plane

      widget_control, info.xmin, get_value = xmin
      widget_control, info.xmax, get_value = xmax
            
      widget_control, info.ymin, get_value = ymin
      widget_control, info.ymax, get_value = ymax
            
      widget_control, info.zmin, get_value = zmin
      widget_control, info.zmax, get_value = zmax

      ; create a structure to hold the zoom information
      zoom = {plane:plane, $
              xmin:xmin,   $
              xmax:xmax,   $
              ymin:ymin,   $
              ymax:ymax,   $
              zmin:zmin,   $
              zmax:zmax}

      problemInfo = {orientation:orientation, $
                     name:current_default_name}

      ; start plotting, based on dimensionality
      if (ndim EQ 1) then begin
        xplot1d_amr, FILE_INFO = fileInfo, $
                     VARIABLE_INFO = variable, $
                     OPTIONS = options, $
                     OUTPUT = output, $
                     ZOOM = zoom, $
                     COLORMAP = current_color_name, $
                     WIDGET_INFORMATION = info, $
                     PARTICLE_WIDGET = particle_widget, $
                     PROBLEM_INFO = problemInfo, $
                     KNOWN_VARIABLES = varnames, $
                     LABEL_OPT = label_info, $
                     TREE_PTR = tree_ptr, $
                     PARAMS_PTR = params_ptr, $
                     DATA_PTR = data_ptr, $
                     PLOT_COUNT = current_number_name, $
                     DEBUG = debug
                
        if (current_number_name EQ '1') then begin
          widget_control, info.query, SENSITIVE=1
        endif else begin
          widget_control, info.query, SENSITIVE=0
        endelse

        widget_control, info.slice1d, SENSITIVE=0

      endif else if (ndim EQ 2) then begin
        xplot2d_amr, FILE_INFO = fileInfo, $
                     READ = iread, $
                     VARIABLE_INFO = variable, $
                     OPTIONS = options, $
                     OUTPUT = output, $
                     VECTOR = vector, $
                     ZOOM = zoom, $
                     COLORMAP = current_color_name, $
                     WIDGET_INFORMATION = info, $
                     CONTOUR_OPT = contours, $
                     PARTICLE_WIDGET = particle_widget, $
                     LABEL_OPT = label_info, $
                     PROBLEM_INFO = problemInfo, $
                     KNOWN_VARIABLES = varnames, $
                     TREE_PTR = tree_ptr, $
                     PARAMS_PTR = params_ptr, $
                     DATA_PTR = data_ptr, $
                     PLOT_COUNT = current_number_name, $
                     DEBUG = debug

        if (current_number_name EQ '1') then begin
          widget_control, info.query, SENSITIVE=1
          widget_control, info.slice1d, SENSITIVE=1
        endif else begin
          widget_control, info.query, SENSITIVE=0
          widget_control, info.slice1d, SENSITIVE=0
        endelse

      endif else if (ndim EQ 3) then begin
        xplot3d_amr, FILE_INFO = fileInfo, $
                     VARIABLE_INFO = variable, $
                     OPTIONS = options, $
                     OUTPUT = output, $
                     PROBLEM_INFO = problemInfo, $
                     LABEL_OPT = label_info, $
                     SLICE_DIR = plane, $ 
                     ZOOM = zoom, $
                     COLORMAP = current_color_name, $
                     WIDGET_INFORMATION = info, $
                     DATA_PTR = data_ptr, $
                     PARAMS_PTR = params_ptr, $
                     XMERGE = xmerge_ptr, $
                     YMERGE = ymerge_ptr, $
                     ZMERGE = zmerge_ptr, $
                     DEBUG = debug

        if (current_number_name EQ '1') then begin
          widget_control, info.slice1d, SENSITIVE=1
        endif else begin
          widget_control, info.slice1d, SENSITIVE=0
        endelse

        widget_control, info.query, SENSITIVE=0
      endif 

    endif else begin ; 'filetype' was not equal to 1

      if Keyword_Set(debug) then print, 'looking to plot .dat'
      ; plot the dat file
      ; output
      widget_control, info.out_buttons, get_value = outputType

      widget_control, info.vsize, get_value = vpixels
      widget_control, info.hsize, get_value = hpixels

      ; create a structure to hold the output information
      output = {type:outputType,         $ ; flag for output medium
                hsize:float(hpixels[0]), $ ; # of pixels for graphic in x
                vsize:float(vpixels[0])} ; # of pixels for graphic in y

      if Keyword_Set(debug) then print, 'formed output'

      ; variable
      ivar = widget_info(info.vars, /droplist_select)

      if Keyword_Set(debug) then print, 'got variable'

      ; plotting options
      widget_control, info.opt_buttons, get_value = opts

      ; data range
      widget_control, info.varMin, get_value = dataMin
      widget_control, info.varMax, get_value = dataMax
            
      widget_control, info.autoData, get_value = auto

      if Keyword_Set(debug) then print, 'got var info'
      variable = {name:varnames(ivar),  $
                  min:dataMin, $
                  max:dataMax, $
                  auto:auto[0]}

      if Keyword_Set(debug) then print, 'formatted variable'

      ; create a structure to hold the options
      options = {blocks:opts[3],   $ ; draw block boundaries
                 log:opts[0],      $ ; take the log of data
                 annotate:opts[4], $ ; plot time and credits
                 colorbar:opts[5], $ ; plot colorbar
                 max:opts[2],      $ ; plot max of variable
                 abs:opts[1],      $ ; take abs value of var
                 tick:opts[6]} ; display tick marks

      if Keyword_Set(debug) then print, 'about to plot'

      plot_dat, filename, $
                OUTPUT=output, $
                VARIABLE=variable, $
                OPTIONS=options

      if Keyword_Set(debug) then print, 'done plotting'
            
    endelse
        
  end

  ;----------------------------------------------------------------------------
  ; particle options pressed
  ;----------------------------------------------------------------------------
  'partopt': begin
    xparticle
  end

  ;-----------------------------------------------------------------------------
  ; histogram pressed
  ;-----------------------------------------------------------------------------
  'histogram': begin

    ; file
    widget_control, info.fstartSuffix, get_value = startSuffix
    startSuffix = fix(startSuffix[0])

    widget_control, info.fendSuffix, get_value = endSuffix
    endSuffix = fix(endSuffix[0])
    if (endSuffix LT startSuffix) then endSuffix = startSuffix

    widget_control, info.fstep, get_value = step
    step = fix(step[0])
    if (step EQ 0) then step = 1

    fileInfo = {prototype:filename, $
                startSuffix:startSuffix, $
                endSuffix:endSuffix, $
                step:step}

    ; output
    widget_control, info.out_buttons, get_value = outputType

    widget_control, info.vsize, get_value = vpixels
    widget_control, info.hsize, get_value = hpixels

    ; create a structure to hold the output information
    output = {type:outputType,         $ ; flag for output medium
              hsize:float(hpixels[0]), $ ; # of pixels for graphic in x
              vsize:float(vpixels[0])} ; # of pixels for graphic in y

    ; variable
    ivar = widget_info(info.vars, /droplist_select)

    ; data range
    widget_control, info.varMin, get_value = dataMin
    widget_control, info.varMax, get_value = dataMax

    widget_control, info.autoData, get_value = auto

    pblm_ctr_lim[ivar,0] = dataMin
    pblm_ctr_lim[ivar,1] = dataMax

    variable = {name:varnames(ivar),  $
                min:dataMin[0], $
                max:dataMax[0], $
                auto:auto[0]}

    widget_control, info.opt_buttons, get_value = opts

    options = {blocks:opts[3],   $ ; draw block boundaries
               log:opts[0],      $ ; take the log of data
               annotate:opts[4], $ ; plot time and credits
               colorbar:opts[5], $ ; plot colorbar
               max:opts[2],      $ ; plot max of variable
               abs:opts[1],      $ ; take abs value of var
               tick:opts[6]} ; display tick marks


    hist_driver, FILE_INFO = fileInfo, $
                 VARIABLE_INFO = variable, $
                 OPTIONS = options, $
                 NBINS = nbins, $
                 HIST_SCALE = max_scale, $
                 OUTPUT = output, $
                 KNOWN_VARIABLES = varnames        

  end


  ;-----------------------------------------------------------------------------
  ; query pressed
  ;-----------------------------------------------------------------------------
  'query': begin

    widget_control, info.status, $
                    set_value = 'click in the domain to get info'

    case orientation of
      0: cursor, x_query, y_query
      1: cursor, y_query, x_query
      2: cursor, x_query, y_query
    endcase

    query, X_CURS = x_query, Y_CURS = y_query, WIDGET_INFORMATION = info, $
           TREE_PTR = tree_ptr, DATA_PTR = data_ptr, PARAMS_PTR = params_ptr, $
           ORIENTATION = orientation, VAR_NAMES = varnames

    widget_control, info.status, set_value = 'status: awaiting orders'

  end

  ;-----------------------------------------------------------------------------
  ; zoombox pressed
  ;-----------------------------------------------------------------------------
  'zoombox': begin

    widget_control, info.status, $
                    set_value = 'left click = move; middle click = resize; right click = done'

    box_cursor, x0, y0, nx, ny
        
    p0_zoom = convert_coord([x0], [y0], /device, /to_data)
    p1_zoom = convert_coord([x0+nx], [y0+ny], /device, /to_data)

    widget_control, info.xmin, set_value = p0_zoom[0]
    widget_control, info.xmax, set_value = p1_zoom[0]

    widget_control, info.ymin, set_value = p0_zoom[1]
    widget_control, info.ymax, set_value = p1_zoom[1]

    widget_control, info.zmin, set_value = -1
    widget_control, info.zmax, set_value = -1

    widget_control, info.status, $
                    set_value = 'status: awaiting orders'

  end


  ;-----------------------------------------------------------------------------
  ; reset pressed
  ;-----------------------------------------------------------------------------
  'reset': begin

    widget_control, info.xmin, set_value = -1
    widget_control, info.xmax, set_value = -1

    widget_control, info.ymin, set_value = -1
    widget_control, info.ymax, set_value = -1

    widget_control, info.zmin, set_value = -1
    widget_control, info.zmax, set_value = -1

  end

  ;-----------------------------------------------------------------------------
  ; slice pressed
  ;-----------------------------------------------------------------------------
  '1dslice': begin

    widget_control, info.status, $
                    set_value = 'click for 1-d slice: left (vertical); right (horizontal)'

    if Keyword_Set(debug) then print, 'going to get the cursor position'

    iout = 0

    case orientation of

      0: begin
        cursor, x_query, y_query

        button = !MOUSE.BUTTON

        if (button EQ 1) then begin
          dir = 1
        endif else if (button EQ 4) then begin
          dir = 0
        endif else begin
          widget_control, info.status, set_value = 'Error: invalid direction'
          iout = 1
        endelse
                
      end
            
      1: begin
        cursor, y_query, x_query, /DOWN, /DATA

        if Keyword_Set(debug) then print, 'x,y = ', x_query, y_query
        button = !MOUSE.BUTTON
        if (button EQ 1) then begin
          dir = 0
        endif else if (button EQ 4) then begin
          dir = 1
        endif else begin
          widget_control, info.status, set_value = 'Error: invalid direction'
          iout = 1
        endelse

      end
      2: begin
        cursor, x_query, y_query

        button = !MOUSE.BUTTON

        if (button EQ 1) then begin
          dir = 1
        endif else if (button EQ 4) then begin
          dir = 0
        endif else begin
          widget_control, info.status, set_value = 'Error: invalid direction'
          iout = 1
        endelse

      end

    endcase 

    ; make sure we are in the domain

    ; get the limits of the current plot domain
    plt_x_min = !x.crange[0]
    plt_x_max = !x.crange[1]

    plt_y_min = !y.crange[0]
    plt_y_max = !y.crange[1]

    if ((orientation EQ 0 AND ((x_query LT plt_x_min) OR $
                               (x_query GT plt_x_max) OR $
                               (y_query LT plt_y_min) OR $
                               (y_query GT plt_y_max))) OR $
        (orientation EQ 1 AND ((y_query LT plt_x_min) OR $
                               (y_query GT plt_x_max) OR $
                               (x_query LT plt_y_min) OR $
                               (x_query GT plt_y_max))) OR  $
        (orientation EQ 2 AND ((x_query LT plt_x_min) OR $
                               (x_query GT plt_x_max) OR $
                               (y_query GT plt_y_min) OR $
                               (y_query LT plt_y_max)))) $
    then begin

      widget_control, info.status, set_value = 'status: outside of domain'
      iout = 1

    endif

    if (iout EQ 0) then begin

      ivar = widget_info(info.vars, /droplist_select)

      if (ndim EQ 3) then begin
        widget_control, info.slicePlane, get_value = plane

        case plane of
          0: begin    ; x-y plane
            case dir of ; along x direction
              0: begin
                dirTitle = "x"
                widget_control, info.xmin, get_value = coord_min
                widget_control, info.xmax, get_value = coord_max
              end
              1: begin ; along y direction
                dirTitle = "y"
                widget_control, info.ymin, get_value = coord_min
                widget_control, info.ymax, get_value = coord_max
              end
            endcase
          end
                    
          1: begin    ; x-z plane
            case dir of ; along x direction
              0: begin
                dirTitle = "x"
                widget_control, info.xmin, get_value = coord_min
                widget_control, info.xmax, get_value = coord_max
              end
              1: begin ; along z direction
                dirTitle = "z"
                widget_control, info.zmin, get_value = coord_min
                widget_control, info.zmax, get_value = coord_max
              end
            endcase
          end
                    
          2: begin    ; y-z plane
            case dir of ; along y direction
              0: begin
                dirTitle = "y"
                widget_control, info.ymin, get_value = coord_min
                widget_control, info.ymax, get_value = coord_max
              end
              1: begin ; along z direction
                dirTitle = "z"
                widget_control, info.zmin, get_value = coord_min
                widget_control, info.zmax, get_value = coord_max
              end
            endcase
          end
        endcase

        coord_min = coord_min[0]
        coord_max = coord_max[0]

        if Keyword_Set(debug) then print, 'about to extract_line', plane, dir, x_query, y_query
                
        if (params_ptr.geometry EQ "CARTESIAN") then begin
          slice1d = extract3d_line(data_ptr, $
                                   POINT = [x_query, y_query], $
                                   XMERGE = xmerge_ptr, $
                                   YMERGE = ymerge_ptr, $
                                   ZMERGE = zmerge_ptr, $
                                   SLICE_DIR = plane, $
                                   DIRECTION = dir, $
                                   COORDS = coords1d)
        endif else begin
          print, "Can't do 1d slice: in 3d, can only slice if geometry is Cartesian"
        endelse
                
        if (coord_min EQ -1.d0) then coord_min = min(coords1d)
        if (coord_max EQ -1.d0) then coord_max = max(coords1d)

        if ( coord_min LT min(coords1d) AND coord_min NE -1.d0) then coord_min = min(coords1d)
        if ( coord_max GT max(coords1d) ) then coord_max = max(coords1d)
                
        for i = 0, (size(coords1d))[1]-1, 10 do begin
          print, coords1d[i], slice1d[i]
        endfor

        data1dMin = +1e99
        data1dMax = -1e99
                
        for i = 0, (size(coords1d))[1]-1 do begin
          if (coords1d[i] GE coord_min AND coords1d[i] LE coord_max ) then begin
            if ( slice1d[i] LT data1dMin ) then data1dMin = slice1d[i]
            if ( slice1d[i] GT data1dMax ) then data1dMax = slice1d[i]
          endif
        endfor
                
        ; get the plotting options -- log is the 0th one
        widget_control, info.opt_buttons, get_value = opts
                
        log = opts[0]

        if Keyword_Set(debug) then print, 'current window = ', !D.WINDOW
                
        ; save the setting for the master plot window (0)
        p0 = !P & x0 = !X & y0 = !Y
        window, 1
        if Keyword_Set(debug) then print, 'new current window = ', !D.WINDOW

        ; make some plot limits based on this slice
        delta = data1dMax - data1dMin
        data1dMin = data1dMin - 0.2*delta
        data1dMax = data1dMax + 0.2*delta
                
        if Keyword_Set(debug) then print, 'slice dataMin, dataMax: ',data1dMin, data1dMax
                
        case plane of
          0: begin
            x = x_query
            y = y_query
            z = zmerge_ptr[0]
          end 
          1: begin
            x = x_query 
            y = ymerge_ptr[0]
            z = y_query
          end 
          2: begin
            x = xmerge_ptr[0]
            y = x_query
            z = y_query
          end
        end

        title = dirTitle + " slice through x = " + $
                strcompress(string(x, FORMAT='(g12.5)')) + $
                " y = " + $
                strcompress(string(y, FORMAT='(g12.5)')) + $
                " and z = " + $
                strcompress(string(z, FORMAT='(g12.5)'))

        plot, coords1d, slice1d, $
              XRANGE = [coord_min, coord_max], $
              YRANGE = [data1dMin,data1dMax], $
              YSTYLE = 1, XSTYLE = 1, $
              XTITLE = dirTitle + ' [cm]', $
              YTITLE = varnames(ivar), $
              TITLE = title


        ; switch back to the master plot window
        wset, 0
        !P = p0 & !X = x0 & !Y = y0

        widget_control, info.status, set_value = 'status: awaiting orders'

      endif else begin ; 'ndim' was not equal to 3
        if dir EQ 0 then begin
          dirTitle = "x"

          widget_control, info.xmin, get_value = coord_min
          widget_control, info.xmax, get_value = coord_max

        endif else if dir EQ 1 then begin
          dirTitle = "y"

          widget_control, info.ymin, get_value = coord_min
          widget_control, info.ymax, get_value = coord_max

        endif

        coord_min = coord_min[0]
        coord_max = coord_max[0]

        ivar = widget_info(info.vars, /droplist_select)

        if Keyword_Set(debug) then print, 'plot a line through ', x_query, y_query

        dummy_var = create_variable(data_ptr, varnames(ivar),debug)

        if (params_ptr.geometry EQ "CARTESIAN" OR $
            params_ptr.geometry EQ "CYLINDRICAL") then begin

          if Keyword_Set(debug) then print, 'about to extract_line', dir, x_query, y_query
          slice1d = extract_line(dummy_var, $
                                 TREE = tree_ptr, $
                                 PARAMETERS = params_ptr, $
                                 POINT = [x_query,y_query], $
                                 DIRECTION = dir, $
                                 COORDS = coords1d)

        endif else if (params_ptr.geometry EQ "SPHERICAL") then begin

          x_sph = sqrt(x_query^2 + y_query^2)
          y_sph = atan(x_query, y_query)

          slice1d = extract_line_polar(dummy_var, $
                                       TREE = tree_ptr, $
                                       PARAMETERS = params_ptr, $
                                       POINT = [x_sph, y_sph], $
                                       DIRECTION = dir, $
                                       COORDS = coords1d)
        endif

        if (coord_min EQ -1.d0) then coord_min = min(coords1d)
        if (coord_max EQ -1.d0) then coord_max = max(coords1d)

        if ( coord_min LT min(coords1d) AND coord_min NE -1.d0) then coord_min = min(coords1d)
        if ( coord_max GT max(coords1d) ) then coord_max = max(coords1d)

        for i = 0, (size(coords1d))[1]-1, 10 do begin
          print, coords1d[i], slice1d[i]
        endfor

        data1dMin = +1e99
        data1dMax = -1e99

        for i = 0, (size(coords1d))[1]-1 do begin
          if (coords1d[i] GE coord_min AND coords1d[i] LE coord_max ) then begin
            if ( slice1d[i] LT data1dMin ) then data1dMin = slice1d[i]
            if ( slice1d[i] GT data1dMax ) then data1dMax = slice1d[i]
          endif
        endfor

        ; get the plotting options -- log is the 0th one
        widget_control, info.opt_buttons, get_value = opts

        log = opts[0]

        if Keyword_Set(debug) then print, 'current window = ', !D.WINDOW

        ; save the setting for the master plot window (0)
        p0 = !P & x0 = !X & y0 = !Y

        window, 1
        if Keyword_Set(debug) then  print, 'new current window = ', !D.WINDOW

        ; make some plot limits based on this slice
        if (log EQ 0) then begin
          delta = data1dMax - data1dMin
          data1dMin = data1dMin - 0.2*delta
          data1dMax = data1dMax + 0.2*delta
        endif else begin
          data1dMin = alog10(0.8*data1dMin)
          data1dMax = alog10(1.2*data1dMax)
        endelse

        if Keyword_Set(debug) then print, 'slice dataMin, dataMax: ',data1dMin, data1dMax

        title = dirTitle + " slice through x = " + $
                strcompress(string(x_query, FORMAT='(g12.5)')) + $
                " and y = " + $
                strcompress(string(y_query, FORMAT='(g12.5)'))

        if (log EQ 0) then begin
          plot, coords1d, slice1d, $
          XRANGE = [coord_min, coord_max], $
          YRANGE = [data1dMin,data1dMax], $
          YSTYLE = 1, XSTYLE = 1, $
          XTITLE = dirTitle + ' [cm]', $
          YTITLE = varnames(ivar), $
          TITLE = title
        endif else begin
          plot, coords1d, alog10(abs(slice1d)), $
          XRANGE = [coord_min, coord_max], $
          YRANGE = [data1dMin,data1dMax], $
          YSTYLE = 1, XSTYLE = 1, $
          XTITLE = dirTitle + ' [cm]', $
          YTITLE = varnames(ivar), $
          TITLE = title
        endelse

        undefine, dummy_var

        ; switch back to the master plot window
        wset, 0
        !P = p0 & !X = x0 & !Y = y0

        widget_control, info.status, set_value = 'status: awaiting orders'

      endelse
    endif ; 'iout' was not equal to 0
  end ; '1dslice'


  ;-----------------------------------------------------------------------------
  ; exit pressed
  ;-----------------------------------------------------------------------------
  'fexit': widget_control, ev.top, /destroy

  else: ; catch-all case

endcase ; 'uval'

end ; of xflash_event


;==============================================================================
; create the widget that handles 1, 2, and 3-d FLASH data.
;
; Once a prototype file is selected, the main information is read from
; this file to determine the variables stored, existence of particles,
; dimensionality, ...  This information will be used to make different
; parts of the main widget sensitive
;
;==============================================================================

pro xflash3, DEBUG=debug

; common all the controls into the widget
common variables, varnames, native_num

common save_parts, particle_widget
common save_color_info, current_color_id, current_color_name
common save_default_info, current_default_id, current_default_name
common save_number_info, current_number_id, current_number_name

common save_path, current_path
common save_filename, filename_old, filename_last_read, filename_last_read_id, $
                      old_start_sfx, old_end_sfx

print, 'initializing xflash3'

version = !VERSION.RELEASE
print, '... IDL version = ', version

if (version LE 5.3) then begin
    print, '... GIF graphics will be used'
endif else begin
    print, '... PNG graphics will be used'
endelse

; get the xflash directory from the xflash_dir environmental variable
xflash_dir = get_xflash_path()
; find the current path
spawn, 'pwd', current_path
current_path = current_path[0]
print, '... current path is ', current_path
print, '... XFLASH_DIR is ', xflash_dir

; initialize variables
filename_old = ' '
filename_last_read = ' '
filename_last_read_id = ' '

old_start_sfx = -99
old_end_sfx = -99

; if debug is not set, then make it zero
if (NOT Keyword_Set(debug)) then debug = 0

; create the main widget base
mainBase = widget_base(/column, title = 'xflash 3.0', MBAR=bar)

; divide this base into vertical sections
fileBase    = widget_base(mainBase, /column, /frame)
outputBase  = widget_base(mainBase, /column, /frame)
blockDistBase = widget_base(mainBase, /row, /frame)
varBase     = widget_base(mainBase, /row,    /frame)
optBase    = widget_base(mainBase, /column, /frame)
rangeBase   = widget_base(mainBase, /row,    /frame)
zoomBase   = widget_base(mainBase, /column, /frame)
advancedBase    = widget_base(mainBase, /row,    /frame)
endBase    = widget_base(mainBase, /row)
statusBase    = widget_base(mainBase, /row, /frame)


;------------------------------------------------------------------------------
; load the problem specific information and get the problem names
;------------------------------------------------------------------------------
xflash_defaults, 'Generic', $
NUM_DEFAULTS=num_defaults, DEFAULT_NAMES=default_names, $
CONTOURS=contours, NUM_CONTOURS=num_contours, $
VECTOR=vector, $
PARTICLE=particle_widget, $
CTR_LIM=ctr_lim

;------------------------------------------------------------------------------
; create the menus
;------------------------------------------------------------------------------

; ---- file menu --------------------------------------------------------------
fileMenu = widget_button(bar, VALUE='File', /MENU)
fopen = widget_button(fileMenu, VALUE='Open prototype...', UVALUE='fopen')
finfo = widget_button(fileMenu, VALUE='Information', UVALUE='finfo', $
                     SENSITIVE=0)
fexit = widget_button(fileMenu, VALUE='Exit', UVALUE='fexit')


; ---- default menu -----------------------------------------------------------
defaults_menu = widget_button(bar, VALUE='Defaults', /MENU)

def_item = lonarr(num_defaults)

; put a 'x' next to the currently selected problem name.  Pad all
; problems with 2 characters

for i = 0, num_defaults-1 do begin

  if (i EQ 0) then begin
    def_item[i] = widget_button(defaults_menu, $
                                VALUE='x '+ default_names[i], $
                                uvalue='def_item')
  endif else begin
    def_item[i] = widget_button(defaults_menu, $
                                VALUE='  '+ default_names[i], $
                                uvalue='def_item')
  endelse

endfor

; store the value of the widget that is the current problem 
current_default_id = def_item[0]
current_default_name = default_names[0]


; ---- color menu -------------------------------------------------------------
colorMenu = widget_button(bar, VALUE='Colormap', /MENU)

; get the available colors
dummy = color_index('Grayscale', GET_NAMES=colors)

numcolors = (size(colors))[1]

colorItem = lonarr(numcolors)

; the names on the menu will have two extra characters to indicate the 
; current selection

for i = 0, numcolors-1 do begin

  if (i EQ 0) then begin 
    colorItem[i] = widget_button(colorMenu, $
                                 VALUE='x ' + colors[i], $
                                 uvalue='colorItem')
  endif else begin 
    colorItem[i] = widget_button(colorMenu, $
                                 VALUE='  ' + colors[i], $
                                 uvalue='colorItem')
  endelse

endfor

; store the value of the widget that is the current selected color
current_color_id = colorItem[0]
current_color_name = colors[0]


; ---- number menu ------------------------------------------------------------
numberMenu = widget_button(bar, VALUE='X/Y plot count')

xycounts = ['1', '1x2', '2x1', '2x2', '2x3', '3x2', '3x3']

numcounts = (size(xycounts))[1]

xycountItem = lonarr(numcounts)

for i = 0, numcounts-1 do begin

  if (i EQ 0) then begin
    xycountItem[i] = widget_button(numberMenu, $
                                   VALUE = 'x ' + xycounts[i], $
                                   uvalue='xycountItem')
  endif else begin
    xycountItem[i] = widget_button(numberMenu, $
                                   VALUE = '  ' + xycounts[i], $
                                   uvalue='xycountItem')
  endelse
                                   
endfor

current_number_id = xycountItem[0]
current_number_name = xycounts[0]


;------------------------------------------------------------------------------
; file selection options
;------------------------------------------------------------------------------

fprototype = widget_label(fileBase, /ALIGN_LEFT, $
                          value='Prototype File: not yet defined', $
                          /DYNAMIC_RESIZE)

suf_base = widget_base(fileBase, /row)

fstartSuffix = cw_field(suf_base, title = 'suffix: ', value = '1', $
                        xsize = 4, uvalue = 'stsfx')
fendSuffix   = cw_field(suf_base, title = 'to', value = '1', $
                        xsize = 4, uvalue = 'endsfx')
fstep        = cw_field(suf_base, title = '   step', value = '1', $
                        xsize = 4, uvalue = 'step')

widget_control, fstartSuffix, sensitive=0
widget_control, fendSuffix, sensitive=0
widget_control, fstep, sensitive=0


;------------------------------------------------------------------------------
; output options
;------------------------------------------------------------------------------

if (!VERSION.RELEASE LE 5.3) then begin
  output = ['screen', 'postscript', 'gif']
endif else begin
  output = ['screen', 'postscript', 'png']
endelse

out_buttons = cw_bgroup(outputBase, output, column = 3, $
                        label_left = 'Output: ', $
                        /exclusive, set_value = 0, uvalue = 'output')

; create a sub-base for the basename and suffix range to plot
size_base = widget_base(outputBase, /row)

hsize = cw_field(size_base, title = 'Plot size, horizontal: ', $
                 value = 800, xsize = 5, uvalue = 'hsize')
vsize = cw_field(size_base, title = ' vertical: ', value = 600, $
                 xsize = 5, uvalue = 'vsize')

widget_control, out_buttons, sensitive=0
widget_control, hsize, sensitive=0
widget_control, vsize, sensitive=0


;------------------------------------------------------------------------------
; variables
;------------------------------------------------------------------------------

varnames = ['----']
vars = widget_droplist(varBase, title = 'Variables: ', uvalue = 'var', $
                       value = varnames)

widget_control, vars, sensitive = 0


;------------------------------------------------------------------------------
; data range
;------------------------------------------------------------------------------

varMin = cw_field(rangeBase, title = 'Data range: ', $
                  xsize = 12, uvalue = 'varMin', value = ctr_lim[0,0])

varMax = cw_field(rangeBase, title = 'to: ', xsize = 12, uvalue = 'varMax', $
                  value = (ctr_lim[0,1]))

widget_control, varMin, sensitive=0
widget_control, varMax, sensitive=0

dataOpt = ['auto']

autoData = cw_bgroup(rangeBase, dataOpt, $
                     /nonexclusive, $
                     set_value = [1], uvalue = 'auto')

widget_control, autoData, sensitive = 0

;------------------------------------------------------------------------------
; parallel block distribution
;------------------------------------------------------------------------------

dataOpt = ['Enable']
enableBlockDist = cw_bgroup(blockDistBase, dataOpt, $
                            label_left = 'Parallel Block Distribution: ', $
                            /nonexclusive, $
                            set_value=[0], $
                            uvalue = 'enableBlockDist')
widget_control, enableBlockDist, sensitive=0

levels = ['----']
levelDroplist = widget_droplist(blockDistBase, title = 'Level: ', uvalue = 'blockDistLevel', $
                                value = levels)

showProcNums = cw_bgroup(blockDistBase, ['Show Proc Numbers'], $
                         /nonexclusive, $
                         set_value=[0], $
                         uvalue = 'showProcNums')

widget_control, levelDroplist, sensitive = 0
widget_control, showProcNums, sensitive = 0

;------------------------------------------------------------------------------
; general plot options
;------------------------------------------------------------------------------

options = ['log ', $
           'abs. value', $
           'max', $
           'show blocks', $
           'annotate', $
           'colorbar', $
           'show ticks']

opt_buttons = cw_bgroup(optBase, options, column = 4, $
                        label_left = 'Options: ', $
                        /nonexclusive, set_value = [0, 0, 0, 0, 1, 1, 1], $
                        uvalue = 'opt')

widget_control, opt_buttons, sensitive = 0


;------------------------------------------------------------------------------
; zoom
;------------------------------------------------------------------------------

slice_options = ['x-y', 'x-z', 'y-z']

slicePlane = cw_bgroup(zoomBase, slice_options, column = 3, $
                       label_left = 'Slice Plane: ', $
                       /exclusive, set_value = 0, uvalue = 'slice')

zoomLabel = widget_label(zoomBase, /align_left, value = 'Zoom:' + $
                         ' (set = -1 for default)')


xbase = widget_base(zoomBase, /row)
ybase = widget_base(zoomBase, /row)
zbase = widget_base(zoomBase, /row)
zbuttonbase = widget_base(zoomBase, /row)

xmin = cw_field(xbase, title = '  xrange: ', xsize = 9, $
                uvalue = 'xmin', value=-1)

xmax = cw_field(xbase, title = 'to ',        xsize = 9, $
                uvalue = 'xmax', value=-1)


ymin = cw_field(ybase, title = '  yrange: ', xsize = 9, $
                uvalue = 'ymin', value=-1)

ymax = cw_field(ybase, title = 'to ',        xsize = 9, $
                uvalue = 'ymax', value=-1)


zmin = cw_field(zbase, title = '  zrange: ', xsize = 9, $
                uvalue = 'zmin', value=-1)

zmax = cw_field(zbase, title = 'to ',        xsize = 9, $
                uvalue = 'zmax', value=-1)


zoombox = widget_button(zbuttonbase, value = 'Zoom Box', uvalue = 'zoombox')
reset = widget_button(zbuttonbase, value = 'Reset', uvalue = 'reset')

widget_control, xmin, sensitive = 0
widget_control, xmax, sensitive = 0

widget_control, ymin, sensitive = 0
widget_control, ymax, sensitive = 0

widget_control, zmin, sensitive = 0
widget_control, zmax, sensitive = 0

widget_control, slicePlane, sensitive = 0
widget_control, zoomLabel, sensitive = 0

widget_control, zoombox, sensitive = 0
widget_control, reset, sensitive = 0


;------------------------------------------------------------------------------
; advanced options
;------------------------------------------------------------------------------

contourOpts = widget_button(advancedBase, value = 'Contour Options', $
                            uvalue = 'ctropt')
widget_control, contourOpts, sensitive = 0


vectorOpts = widget_button(advancedBase, value = 'Vector Options', $
                           uvalue = 'vecopt')
widget_control, vectorOpts, sensitive = 0

particleOpts = widget_button(advancedBase, value = 'Particle Options', $
                             uvalue = 'partopt')
widget_control, particleOpts, sensitive = 0


histogramOpts = widget_button(advancedBase, value = 'Histogram Options', $
                              uvalue = 'histopt')
widget_control, histogramOpts, sensitive = 0

labelOpts = widget_button(advancedBase, value = 'Floating Label', $
                          uvalue = 'labelopt')

widget_control, labelOpts, sensitive = 0

;------------------------------------------------------------------------------
; plot buttons and status
;------------------------------------------------------------------------------

plot = widget_button(endBase, value = 'Plot', uvalue = 'plot')
histogram = widget_button(endBase, value = 'Histogram', uvalue = 'histogram')
query = widget_button(endBase, value = 'Query', uvalue = 'query')
slice1d = widget_button(endBase, value = '1-d Slice', uvalue = '1dslice')

widget_control, query, sensitive = 0
widget_control, slice1d, sensitive = 0
widget_control, plot, sensitive = 0
widget_control, histogram, sensitive = 0

status = widget_label(statusBase, value = 'status: ', $
                      uvalue = 'status', /align_left, /dynamic_resize)


;------------------------------------------------------------------------------
; draw the widget
;------------------------------------------------------------------------------
widget_control, mainBase, /realize


; setup a structure to hold the widget information
info = {mainBase:mainBase, $        ; main base
        fileBase:fileBase, $        ; file base widget
        fileMenu:fileMenu, $
        fopen:fopen, $
        finfo:finfo, $
        fexit:fexit, $
        fprototype:fprototype, $
        fstartSuffix:fstartSuffix, $  ; field widget for starting suffix
        fendSuffix:fendSuffix, $      ; field widget for ending suffix
        fstep:fstep, $                ; field widget for step
        colorMenu:colorMenu, $
        colorItem:colorItem, $
        numberMenu:numberMenu, $
        xycountItem:xycountItem, $
        outputBase:outputBase, $    ; output base widget
        out_buttons:out_buttons, $    ; button group for output
        hsize:hsize, $                ; field widget for horizontal image size
        vsize:vsize, $                ; field widget for vertical image size
        varBase:varBase, $          ; variable base widget
        vars:vars, $                  ; droplist widget for variable selection
        blockDistBase:blockDistBase, $  ; processor distribution base widget
        enableBlockDist:enableBlockDist, $ ; the button to enable block distribution
        levelDroplist:levelDroplist, $ ; droplist widget for proc dist level selection
        showProcNums:showProcNums, $ ; the button to enable proc numbers 
        optBase:optBase, $          ; option base widget
        opt_buttons:opt_buttons, $    ; button group for options
        rangeBase:rangeBase, $          ; contour base widget
        varMin:varMin, $              ; field widget for min contour level
        varMax:varMax, $              ; field widget for max contour level
        autoData:autoData, $
        advancedBase:advancedBase, $  
        contourOpts:contourOpts, $          ; contour options button
        vectorOpts:vectorOpts, $          ; vector options button
        particleOpts:particleOpts, $          ; particle options button
        histogramOpts:histogramOpts, $
	labelOpts:labelOpts, $
        zoomBase:zoomBase, $        ; zoom base widget
        slicePlane:slicePlane, $
        zoomLabel:zoomLabel, $
        xmin:xmin, $                ; field widget for min x coord
        ymin:ymin, $                ; field widget for min y coord
        xmax:xmax, $                ; field widget for max x coord
        ymax:ymax, $                ; field widget for max y coord
        zmin:zmin, $
        zmax:zmax, $
        zoombox:zoombox, $
        reset:reset, $
        endBase:endBase, $          ; end base widget
        plot:plot, $                  ; button for plotting
        histogram:histogram, $
        query:query, $                ; button for querying
        slice1d:slice1d, $                ; button for slicing
        statusBase:statusBase, $
        status:status, $                ; status bar
        debug:debug }                ; debug information

; register the info structure with the widget base
widget_control, mainBase, set_uvalue=info, /no_copy

xmanager, 'xflash', /no_block, mainBase

end
