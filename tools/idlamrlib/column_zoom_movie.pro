pro column_zoom_movie, amr, plane, filebase=filebase, zooms=zooms, $
                       outputmethod=outputmethod, $ 
                       sranges=sranges, target=target,$
                       framesonly=framesonly, $
                       anisotropic=anisotropic, $
                       xrange=xrange, yrange=yrange, nocolorbar=nocolorbar,	$
                       verbose=verbose, scalerange=scalerange,			$
                       maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
                       pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
                       smeardist=smeardist, colort=colort,_extra = extra

; This program makes a movie of zooming in on a column-density plot
; zooms is a N-element fltarr of zoom factors relative to full-zoomout
; sranges is a Nx2 fltarr of min-max values of the colormap
; The output is N .ppm files with the base name of filebase
; target is an optional coordinate of where to zoom towards
; All options of column_amr are supported, although scalerange is
; obsolete, replaced by sranges
; and xrange,yrange are replaced by target and zooms

;pro 
if keyword_set(colort) then begin
  loadct,colort
  ctr = intarr(256)
  ctg = intarr(256)
  ctb = intarr(256)
  tvlct,ctr,ctg,ctb,/get
  if colort eq 5.1 then begin
      loadct,5
      tvlct,ctr,ctg,ctb,/get
      ctr[150:255] = 255
  endif
endif

size_zooms = size(zooms)
size_sranges = size(sranges)
error=0
if size_zooms[0] ne 1 then begin
    print, 'Error: zooms must be a 1-D array.'
    error=1
endif

if size_sranges[0] ne 2 then begin
    print, 'Error: sranges must be a 2-D array.'
    error=1
endif
if size_sranges[2] ne 2 then begin
    print, 'Error: sranges must have a dimension of X by 2.'
    error=1
endif
if size_zooms[1] ne size_sranges[1] then begin
    print, 'Error: zooms and sranges must be the same size.'
    error=1
endif
if error ne 0 then begin
    print, 'Exiting because of error.'
    print, 'size_zooms = ', + size_zooms
    print, 'size_sranges = ' + size_sranges
    print, 'plane = ' + plane
    print, 'filebase = ' + filebase
    print, 'zooms = ' + zooms
    print, 'sranges = ' + sranges
    print, 'target = ' + target
    return
endif

outM = 'x' ; output method
outS = ''  ; file suffix
outD = 'x' ; plotting method for set_plot
if keyword_set(outputmethod) then begin
    if outputmethod eq 'ps' then begin
        outM = 'ps'
        outS = '.ps'
        outD = 'ps'
    endif else if outputmethod eq 'ppm' then begin
        outM = 'ppm'
        outS = '.ppm'
        outD = 'x'
    endif else if outputmethod eq 'x' then begin
        outM = 'x'
        outS = ''
        outD = 'x'
    endif else if outputmethod eq 'cgm' then begin
        outM = 'cgm'
        outS = '.cgm'
        outD = 'cgm'
    endif else if outputmethod eq 'z' then begin
        outM = 'z'
        outS = '.ppm'
        outD = 'z'
    endif else if outputmethod eq 'jpeg' then begin
        outM = 'jpeg'
        outS = '.jpeg'
        outD = 'z'
    endif
endif

size = size_zooms[1]

if plane eq 0 then begin
    dim1=1
    dim2=2
    dim3=0
endif else if plane eq 1 then begin
    dim1=0
    dim2=2
    dim3=1
endif else if plane eq 2 then begin
    dim1=0
    dim2=1
    dim3=2
endif else begin
    print, 'Error: plane must be 0, 1, or 2'
    return
endelse

if keyword_set(target) then begin
    targ = target
endif else begin
    targ = [ $ 
             (amr.boxmin[dim1] + amr.boxmax[dim1])/2.0, $
             (amr.boxmin[dim2] + amr.boxmax[dim2])/2.0  ]
endelse
boxsize = [ $ 
             (-amr.boxmin[dim1] + amr.boxmax[dim1]), $
             (-amr.boxmin[dim2] + amr.boxmax[dim2])  ]

if keyword_set(framesonly) then begin
    istart = framesonly[0]
    iend = framesonly[1]
endif else begin
  istart = 0
  iend = size-1
endelse
for i=istart, iend do begin
    if keyword_set(filebase) then begin
        outfile = filebase        
    endif else begin
        outfile = 'Column_Zoom_Movie_Frame'
    endelse

    if i lt 1000 then outfile = outfile + '0'
    if i lt 100 then outfile = outfile + '0'
    if i lt 10 then outfile = outfile + '0'
    outfile = strcompress(outfile + string(i)+outS,/remove_all)

    print,'Outfile = '+outfile
    xr = boxsize[0] / zooms[i]
    yr = boxsize[1] / zooms[i]
    xbound = [targ[0] - xr/2.0, targ[0] + xr/2.0]
    ybound = [targ[1] - yr/2.0, targ[1] + yr/2.0]
    print, 'xr,'+string(xr)+':yr,'+string(yr)

    ; move clipping box to better fit domain
    if xbound[0] lt amr.boxmin[dim1] then begin
        ; shift right
        print, 'shift right ' + string(error)
        error = amr.boxmin[dim1] - xbound[0]
        xbound = xbound + error
    endif else if xbound[1] gt amr.boxmax[dim1] then begin
        ; shift left
        print, 'shift left ' + string(error)
        error = amr.boxmax[dim1] - xbound[1]
        xbound = xbound + error
    endif
    if ybound[0] lt amr.boxmin[dim2] then begin
        ; shift up
        print, 'shift up ' + string(error)
        error = amr.boxmin[dim2] - ybound[0]
        ybound = ybound + error
    endif else if ybound[1] gt amr.boxmax[dim2] then begin
        ; shift down
        print, 'shift down ' + string(error)
        error = amr.boxmax[dim2] - ybound[1]
        ybound = ybound + error
    endif

    ; reduce clipping box to not exceed domain
    if xbound[0] lt amr.boxmin[dim1] then begin
        xbound[0] = amr.boxmin[dim1]
        print, 'clip left ' + string(amr.boxmin[dim1] - xbound[0])
    endif 
    if xbound[1] gt amr.boxmax[dim1] then begin
        xbound[1] = amr.boxmax[dim1]
        print, 'clip right ' + string(xbound[1] - amr.boxmax[dim1])
    endif
    if ybound[0] lt amr.boxmin[dim2] then begin
        ybound[0] = amr.boxmin[dim2]
        print, 'clip bottom ' + string(amr.boxmin[dim2] - ybound[0])
    endif
    if ybound[1] gt amr.boxmax[dim2] then begin
        ybound[1] = amr.boxmax[dim2]
        print, 'clip top ' + string(ybound[1] - amr.boxmax[dim2])
    endif
    print, 'Frame='+string(i)+', xbound=' + string(xbound[0]) + ', ybound=' + string(ybound[0])
    
    ; Plot to the appropriate device
    ; If plotting to offscreen, plot to X occasionally to show progress
    set_plot,outD,/copy
    if (outM eq 'ps' or outM eq 'cgm') then begin
        device,file=outfile,/color
    endif

    column_amr, amr, plane, anisotropic=anisotropic,		$
      xrange=xbound, yrange=ybound, nocolorbar=nocolorbar,	$
      verbose=verbose, scalerange=sranges[i,*],			$
      maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
      pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
      smeardist=smeardist, _extra = extra

    if (outM eq 'ps' or outM eq 'cgm') then begin
        device,/close
    endif
    if (outM eq 'ppm' or outM eq 'z' or outM eq 'jpeg') then begin
        ;thisframe = tvrd(true=1,/order)
        thisframe = tvrd(/order)
        if keyword_set(colort) then begin
            monoframe = thisframe
            sizimg1 = size(monoframe)
            thisframe = intarr(sizimg1[1],sizimg1[2],3,/nozero)
            thisframe[*,*,0] = ctr[monoframe]
            thisframe[*,*,1] = ctg[monoframe]
            thisframe[*,*,2] = ctb[monoframe]
            tv,thisframe,true=3,/order            
        endif
    endif
    if (outM eq 'ppm' or outM eq 'z') then begin
        write_ppm,outfile,thisframe
    endif
    if outM eq 'jpeg' then begin
        if keyword_set(colort) then begin
            write_jpeg,outfile,thisframe,quality=100,/order,true=3
        endif else begin
            write_jpeg,outfile,thisframe,quality=100,/order
        endelse
    endif

    ; convert ps files into ppms (with white background)
    ;spawnstr = 'gs -sDEVICE=ppm -sOutputFile=' $
    ;  + outfile+'.ppm' $
    ;  + ' -dBATCH -dNOPAUSE ' $
    ;  + outfile+'.ps' $
    ;  + ' quit'
    ;
    ;print,spawnstr
    ;spawn,spawnstr
    if outS eq '.ppm' then begin
        spawnstr = 'gzip ' + outfile
        print,spawnstr
        spawn,spawnstr
    endif

    if (outM ne 'x' and outM ne 'ppm') then begin
    ; periodically refresh the screen as a progress indicator
        if (i mod 10) eq 0 then begin
            set_plot,'x'
            column_amr, amr, plane, anisotropic=anisotropic,		$
              xrange=xbound, yrange=ybound, nocolorbar=nocolorbar,	$
              verbose=verbose, scalerange=sranges[i,*],			$
              maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
              pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
              smeardist=smeardist, _extra = extra
            set_plot,outD
        endif
    endif

    
    
    ;if keyword_set(filebase) then begin
    ;    frame = tvrd(true=1,/order)
    ;    write_ppm,outfile,frame
    ;endif
    ;device,/close
    ;set_plot,'x'

    ; occasionally display to screen also
    


endfor
end

; END of COLUMN_ZOOM_MOVIE


