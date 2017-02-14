pro raster_zoom_movie, amr, plane, slice, filebase=filebase, zooms=zooms, $ 
                       sranges=sranges, target=target,$
                       anisotropic=anisotropic, $
                       xrange=xrange, yrange=yrange, nocolorbar=nocolorbar, $
                       verbose=verbose, zrange=zrange, $
                       maxlevel=maxlevel, localrange=localrange, $
                       snaptogrid=snaptogrid, _extra = extra

; This program makes a movie of zooming in on a column-density plot
; zooms is a N-element fltarr of zoom factors relative to full-zoomout
; sranges is a Nx2 fltarr of min-max values of the colormap
; The output is N .ps files with the base name of filbase
; target is an optional coordinate of where to zoom towards
; All options of column_amr are supported, although scalerange is
; obsolete, replaced by sranges
; and xrange,yrange are replaced by target and zooms

;pro 

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

for i=0, size-1 do begin
    outfile = filebase + string(i) + '.ps'
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
    
    set_plot,'ps'
    device,file=outfile
    column_amr, amr, plane, anisotropic=anisotropic,		$
	xrange=xbound, yrange=ybound, nocolorbar=nocolorbar,	$
	verbose=verbose, scalerange=sranges[i,*],			$
	maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
	pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
	smeardist=smeardist, _extra = extra
    device,/close
    set_plot,'x'

    ; occasionally display to screen also
    if (i mod 10) eq 0 then begin
        column_amr, amr, plane, anisotropic=anisotropic,		$
          xrange=xbound, yrange=ybound, nocolorbar=nocolorbar,	$
          verbose=verbose, scalerange=sranges[i,*],			$
          maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
          pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
          smeardist=smeardist, _extra = extra
    endif


endfor
end

; END of COLUMN_ZOOM_MOVIE
