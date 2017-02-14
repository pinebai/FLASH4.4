; col_any_frame,tree=tree,angles=[0,0,0],xrange=xr,yrange=yr,zrange=zr,outM='X',target=targ,maxlevel=ml,log=1

pro col_any_frame,$
                  tree=tree, $ ; amr tree to pass to column_any
                  angles=angles, $ ; set of angles to pass to column_any
                  axisoff=axisoff, $ ; true if you want the simulation axes drawn
                  xrange=xrange, $ ; plot window ranges
                  yrange=yrange, $
                  zrange=zrange, $
                  colorbar=colorbar, $ ; true if you want a colorbar
                  colortable=colortable, $ ; set to override the current colortable with another
                  scalerrange=scalerrange, $ ; force image values to be clamped to [min,max]
                  log=log, $ ; take the log of the image
                  outM=outM, $ ; string to identify which type of output
                  filename=filename, $ ; eg 'file0321.jpeg'
                  aftersmooth=aftersmooth, $ ;option to perform a smoother on the rendered image
                  _extra=extra

tstart = systime(1)
dw = !d.window
if not keyword_set(tree) then begin
    print,'Error in COL_ANY_FRAME!  TREE not specified!'
    return
endif

if not keyword_set(angles) then begin
    angles = double([0.0,0.0,0.0])
endif

if not keyword_set(outM) then outM = 'X'

if not keyword_set(xrange) then begin
    xrange=[ tree.boxmin[0],tree.boxmax[0] ]
endif else begin
    
endelse
if not keyword_set(yrange) then begin
    yrange=[ tree.boxmin[1],tree.boxmax[1] ]
endif else begin
    
endelse
if not keyword_set(zrange) then begin
    zrange=[ tree.boxmin[2],tree.boxmax[2] ]
endif else begin
    
endelse

imgsize_raster = [400,400,400]
imgsize_vector = 32768.0
if outM eq 'X' then begin
    imgsize=imgsize_raster
    if !d.window eq -1 then begin
        window,retain=2,xsize=640,ysize=480
    endif
    scalepix = 0
endif else if outM eq 'PS' then begin
    imgsize=imgsize_vector
    scalepix = 1
    if outM eq 'PS' and keyword_set(filename) then begin
        set_plot,'ps',/copy
        device,file=filename,/color
    endif
endif else if outM eq 'PPM' then begin
    imgsize=imgsize_raster
    if !d.window eq -1 then window, 0,retain=2,xsize=640,ysize=480
    scalepix = 0
    if !d.window eq -1 then window, 0,retain=2,xsize=640,ysize=480
endif else if outM eq 'JPEG' then begin
    imgsize=imgsize_raster
    scalepix = 0
    if !d.window eq -1 then window, 0,retain=2,xsize=640,ysize=480
endif

; find plotting device aspect_ratio
; so we can use /normal coords and get squares
AR = double(!d.x_size) / double(!d.y_size)
if AR ge 1.3 then begin
    boxsize = 0.85*double(!d.y_size)
endif else begin
    boxsize = 0.6*AR*double(!d.y_size)
endelse
if boxsize gt 2048 then begin
    print,'Warning, image size of '+string(boxsize)+' is being constrained to 2048'
    boxsize=2048
endif
print,'Boxsize = '+string(boxsize)
image = column_any(tree,angles,xrange=xrange,yrange=yrange,zrange=zrange,imgdim=[boxsize,boxsize,1],_extra=extra)

erase

; post-processing of image
if keyword_set(aftersmooth) then image = smooth(image,aftersmooth,/edge_truncate)

; clamp values to a range
if keyword_set( scalerrange) then begin
    wh = where(image lt scalerrange[0])
    if wh[0] ne -1 then image[wh] = scalerrange[0]
    wh = where(image gt scalerrange[1])
    if wh[0] ne -1 then image[wh] = scalerrange[1]
endif else begin
    ; if no scalerrange specified, treat zero as NaN
    ;print,bork
    wh = where(image gt double(0.0))
    mins = min(image[wh])
    wh = where(image lt mins)
    if wh[0] ne -1 then image[wh] = mins
    scalerrange=[mins,max(image)]
endelse
if scalerrange[0] gt scalerrange[1] then begin
    temp = scalerrange[0]
    scalerrange[0] = scalerrange[1]
    scalerrange[1] = temp
endif

if scalerrange[0] eq scalerrange[1] then begin
    scalerrange[0] = 0.95*scanerrange[0]
    scalerrange[1] = 1.05*scalerrange[1]
endif



if keyword_set(log) then begin
    image = alog10(image)
    scalerrange = alog10(scalerrange)
endif

image = bytscl(image) 
c = intarr(256,3) 
tvlct,c,/get	 
imc = [[[image]],[[image]],[[image]]] 
for ci=0,2 do imc[*,*,ci] = c[image[*,*],ci] 




pbang=!p.multi
!p.multi = [0,1,1,2] ; get to refresh screen, and call plot twice
; place image where it belongs
;print,bork
tv,imc,0.08,0.1,true=3,/normal



; stack plots in the z-dim of the frame


; draw box for image
if AR ge 1.3 then begin
    plot, [0,0], [0,0], /nodata, xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,position=[0.08,0.1,0.08+0.85/AR,0.1+0.85],/normal,/noerase
    boxsize = 0.85*double(!d.y_size)
endif else begin
    plot, [0,0], [0,0], /nodata, xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,position=[0.08,0.1,0.08+0.60,0.1+0.60*AR],/normal,/noerase
    boxsize = 0.6*AR*double(!d.y_size)
endelse


;!p.multi=pbang

; draw colorbar
if not keyword_set(axisoff) then begin
    if (abs(scalerrange[1]) lt 1.0e2) and (abs(scalerrange[1]) gt 1.0e-1) then $
      format = '(f9.2)' $
    else $
      format = '(e10.2)'

    if AR ge 1.3 then begin
        colorbar,position=[0.83,0.4,0.95,0.95],range=scalerrange,/vertical,format=format
    endif else begin
        colorbar,position=[0.83,0.1+0.3*AR,0.95,0.1+0.6*AR],range=scalerrange,/vertical,format=format
    endelse
endif

; draw box for little axes
absize=0.12
abcenter = [0.83+absize/2.0,0.1+absize*AR/2.0,0.0]
plot, [0,0], [0,0], /nodata,xrange=[-1.0*(xrange[1]-xrange[0])/8.0,(xrange[1]-xrange[0])/8.0], yrange=[-1.0*(yrange[1]-yrange[0])/8.0,(yrange[1]-yrange[0])/8.0],position=[0.83,0.1,0.83+absize,0.1+absize*AR],/normal,xticks=2,yticks=2,xstyle=1,ystyle=1,/noerase
;!p.multi=pbang

azimuth = angles[0]
altitude = angles[1]
pitch = angles[2]
; the identity quaternion
QI = [[1,0,0],$
      [0,1,0],$
      [0,0,1]]
; change camera azimuth angle
QA = [[ cos(azimuth),sin(azimuth),0],$
      [-sin(azimuth),cos(azimuth),0],$
      [            0,           0,1]]
; change camera polar angle
QP = [[1,            0,             0],$
      [0, cos(altitude),sin(altitude)],$
      [0,-sin(altitude),cos(altitude)]]
; roll the camera angle
QR = [[ cos(pitch),sin(pitch),0],$
      [-sin(pitch),cos(pitch),0],$
      [          0,         0,1]]
QROT = QR##(QP##(QA##QI)) ; rotate sim space to camera space (no transforms)
QXtoY = [[0,0,1],$
         [1,0,0],$
         [0,1,0]]
QXtoZ = [[0,1,0],$
         [0,0,1],$
         [1,0,0]]
QfixAR = [[1,0,0],$
          [0,AR,0],$
          [0,0,1]]

; normalized coords to draw an x-axis arrow
xarrow = [ [0.0, 0.0, 0.0],$  ; 0
           [1.0, 0.0, 0.0], $ ; 1
           [0.8, 0.1, 0.0],$  ; 2
           [0.8,-0.1, 0.0],$  ; 3
           [1.0, 0.0, 0.0],$  ; 4
           [0.8, 0.0, 0.1],$  ; 5
           [0.8, 0.0,-0.1],$  ; 6
           [1.0, 0.0, 0.0]]   ; 7

yarrow = fltarr(3,8)
zarrow = fltarr(3,8)


axisboxsize = absize
axisboxcenter = abcenter
for arrowi=0,7 do begin
    yarrow[*,arrowi] = QXtoY##xarrow[*,arrowi]
    zarrow[*,arrowi] = QXtoZ##xarrow[*,arrowi]
endfor
for arrowi=0,7 do begin
    xarrow[*,arrowi] = QROT##xarrow[*,arrowi]
    yarrow[*,arrowi] = QROT##yarrow[*,arrowi]
    zarrow[*,arrowi] = QROT##zarrow[*,arrowi]
endfor
for arrowi=0,7 do begin
    xarrow[*,arrowi] = QfixAR##xarrow[*,arrowi]
    yarrow[*,arrowi] = QfixAR##yarrow[*,arrowi]
    zarrow[*,arrowi] = QfixAR##zarrow[*,arrowi]
endfor

xarrow = xarrow*axisboxsize/2.0 + axisboxcenter#(dblarr(8)+1)
yarrow = yarrow*axisboxsize/2.0 + axisboxcenter#(dblarr(8)+1)
zarrow = zarrow*axisboxsize/2.0 + axisboxcenter#(dblarr(8)+1)
xarrow = xarrow[0:1,*]
yarrow = yarrow[0:1,*]
zarrow = zarrow[0:1,*]
plots,xarrow[0,*],xarrow[1,*],/normal
plots,yarrow[0,*],yarrow[1,*],/normal
plots,zarrow[0,*],zarrow[1,*],/normal
xyouts,xarrow[0,1],xarrow[1,1],'X',/normal,charsize=1.5,color=255.0 + 255.0*256.0
xyouts,yarrow[0,1],yarrow[1,1],'Y',/normal,charsize=1.5,color=255.0*256.0 + 255.0*65536.0
xyouts,zarrow[0,1],zarrow[1,1],'Z',/normal,charsize=1.5,color=255.0*65536.0 + 255.0
;bork

if keyword_set(filename) then begin
    if outM eq 'PPM' or outM eq 'JPEG' then begin
        tvr = tvrd(true=3)
        if outM eq 'PPM' then begin
            ; write a ppm
            write_ppm,filename,tvr
        endif else if OutM eq 'JPEG' then begin
            ; write a jpeg
            write_jpeg,filename,tvr,quality=100,true=3
        endif
    endif else if outM eq 'PS' then begin
        device,/close
        set_plot,'x'
    endif
endif
!p.multi = pbang
!d.window=dw

print,'COL_ANY_FRAME took '+strtrim(string(systime(1)-tstart),2)+' seconds'

end
                  
