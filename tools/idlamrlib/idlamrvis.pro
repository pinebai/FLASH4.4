; This file houses many GUI tools for visualization of AMR data
; The actual visualizarion routines are in other files.  This just
; contains the GUIs and wrappers

; IDLAMRVIS is the root program for launching GUI-based vis tools
; It launches slice, column tools, etc.
; It can be passed a filename,variable pair, or an AMR object created
; with GET_AMRCOMPONENT of created through data manipulation
; If called with no arguments, it prompt the user to select a plotfile
; If an AMR Object is not provided and variable is not set, it
; defaults to the variable 'density'
; INPUTS:
; filename: string value that indictes the plotfile to open
; variable: string value says what variable of the plotfile to read
; amr: in stead of reading a plotfile, use an existing AMR Object in memory
pro idlamrvis,filename=filename,amr=amr,variable=variable,slice=slice,column=column
initialmessage = 'Ready'
if not keyword_set(variable) then begin
    ; no variable requested, default to density
    component = 'density'
endif else begin
    ; user-specified AMR component
    component = variable[0]
endelse
zoom = 1.0
data = {idlamrvis, title:'IDLAMRVis', filename:'', component:component, center:dblarr(3), zoom:zoom, treeptr:ptr_new(), controlptr:ptr_new(), messageptr:ptr_new()}
dataptr = ptr_new(data,/no_copy)
base = WIDGET_BASE(/column, uvalue=dataptr, title=(*dataptr).title)
baseptr = ptr_new(base,/no_copy)
control = WIDGET_BASE(*baseptr,/row, uvalue=dataptr)
controlptr = ptr_new(control,/no_copy)
(*dataptr).controlptr = controlptr
close_base = WIDGET_BUTTON(*controlptr, value='Close', uvalue='Close')
create_column = WIDGET_BUTTON(*controlptr, value='Column Density', uvalue='Column Density')
create_slice = WIDGET_BUTTON(*controlptr, value='Slice', uvalue='Slice')
message = WIDGET_TEXT(*baseptr, event_pro='ReRenderVoxel', value=initialmessage, xsize=80)
messageptr = ptr_new(message,/no_copy)
(*dataptr).messageptr = messageptr
WIDGET_CONTROL, *baseptr, /REALIZE
if not keyword_set(amr) then begin
    if not keyword_set(filename) then begin
        initialmessage = 'Please select a plotfile'
        WIDGET_CONTROL,*messageptr, SET_VALUE=initialmessage
        filename = dialog_pickfile()
    endif else begin
                                ;initialmessage = 'Opening file '+filename
    endelse
    initialmessage = 'Opening file '+filename
    WIDGET_CONTROL,*messageptr, SET_VALUE=initialmessage
    
    (*dataptr).filename = filename
    amr = get_amrcomponent(filename, component)
    initialmessage = 'Making tree structure from plotfile '+filename
endif else begin
    ; AMR object supplied, don't load one
    initialmessage = 'Making tree structure from AMR data'
endelse
center = (amr.boxmax+amr.boxmin)/2.0
(*dataptr).center = center

WIDGET_CONTROL,*messageptr, SET_VALUE=initialmessage

tree = amr_tree(amr,free=1,quiet=1)
treeptr = ptr_new(tree,/no_copy)
(*dataptr).treeptr = treeptr

initialmessage = 'Ready'
WIDGET_CONTROL,*messageptr, SET_VALUE=initialmessage

;if keyword_set(column) then begin
;    WIDGET_CONTROL,*messageptr, SET_VALUE='Building Voxel.  Please be patient...'
;    col_any_widget,(*dataptr).treeptr, group_leader=*baseptr
;    WIDGET_CONTROL,*messageptr, SET_VALUE='Ready'
;endif
;if keyword_set(slice) then begin
;    WIDGET_CONTROL,*messageptr, SET_VALUE='Creating slice...'
;    slice_widget,(*dataptr).treeptr, group_leader=*baseptr
;    WIDGET_CONTROL,*messageptr, SET_VALUE='Ready'
;endif
XMANAGER, 'IDLAMRVis', *baseptr

end

;Event Handler for IDLAMRVis
pro idlamrvis_event,ev
WIDGET_CONTROL, ev.id, GET_UVALUE=uvalue
WIDGET_CONTROL, ev.top, GET_UVALUE=ps
case uvalue of
    'Close' : begin
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Closing application'
        ;wait,5
        WIDGET_CONTROL, ev.top, /DESTROY
    end
    'Column Density' : begin
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Building Voxel.  Please be patient...'
        col_any_widget,(*ps).treeptr, group_leader=ev.top
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Ready'
    end
    'Slice' : begin
        help,ev,/structures
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Creating slice...'
        slice_widget,(*ps).treeptr, group_leader=ev.top
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Ready'
    end
    else : begin
        print,'IDLAMRVis_Event could not handle unknown event '+string(uvalue)
    end
endcase
end

; GUI program/wrapper for visualization using slices through the data
; Spawns a main window with a plot area, and 2 child windows for
; colormap and geometry control
; INPUTS:
; treeptr: a pointer to an AMR Tree structure
; group_leader: Optional widget ID of the window that spawned this
;   Slice_widget
;   Setting group_leader ensures when the parent dies, so does slice_widget
pro slice_widget,treeptr, group_leader=group_leader
;help,*treeptr
x = treeptr
srangemin=double(1.0e-22)
srangemax=double(1.0)
data = {SliceWidget, $
        treeptr:treeptr, $
        drawx:fix(640), drawy:fix(480), ct:intarr(256,3), $
        cdepth:fix(8), $
        tv_offset:double([0.08,0.1]), tv_imgsize:fix([408,408]), $
        xr:dblarr(2), yr:dblarr(2), viewaxis:dblarr(3), target:dblarr(3), $
        angles:dblarr(3), $
        contours:fix(0), $
        dolinear:fix(0), $
        srange:dblarr(2),logsrangemin:long(0), logsrangeminslider:long(0), logsrangemintext:long(0), $
        logsrangemax:long(0), logsrangemaxslider:long(0), logsrangemaxtext:long(0), $
        xcenter:long(0),xcenterslider:long(0),xcentertext:long(0), $
        ycenter:long(0),ycenterslider:long(0),ycentertext:long(0), $
        zcenter:long(0),zcenterslider:long(0),zcentertext:long(0), $
        angform:'(f8.2)', $
        azimuthslider:long(0), azimuthtext:long(0), $
        altitudeslider:long(0), altitudetext:long(0), $
        pitchslider:long(0), pitchtext:long(0), $
        mydw:fix(0) ,$
        zoomrange:dblarr(2), zoom:double(0.0), logzoom:long(0), logzoomtext:long(0), logzoomslider:long(0), $
        rawimageptr:ptr_new(), imageptr:ptr_new(), $
        but_press:intarr(3,2)-1, but_release:intarr(3,2)-1}
dataptr = ptr_new(data,/no_copy)
titlestr = 'Slice '+  (*treeptr).name +' '+ (*treeptr).componentName
if keyword_set(group_leader) then begin
    base = WIDGET_BASE(/column, uvalue=dataptr, title=titlestr, group_leader=group_leader)
endif else begin
    base = WIDGET_BASE(/column, uvalue=dataptr, title=titlestr)
endelse
baseptr = ptr_new(base,/no_copy)
control = WIDGET_BASE(*baseptr,/row, uvalue=dataptr)
controlptr = ptr_new(control,/no_copy)
close = WIDGET_BUTTON(*controlptr, VALUE='Close', UVALUE='Close')
write_jpeg = WIDGET_BUTTON(*controlptr, VALUE='Write JPEG', UVALUE='Write JPEG')
write_ppm = WIDGET_BUTTON(*controlptr, VALUE='Write PPM', UVALUE='Write PPM')
view = WIDGET_BASE(*baseptr,/column, uvalue=dataptr)
viewptr = ptr_new(view,/no_copy)

gsr = (*(*dataptr).treeptr).srange
gsr = alog10(gsr)
srangecon = WIDGET_BASE(/column, uvalue=dataptr, title='Slice Scalar Range', group_leader=*baseptr)
logsrangemin = WIDGET_BASE(srangecon,/row, uvalue='logsrangemin')
logsrangeminslider = CW_FSLIDER(logsrangemin, value=gsr[0], minimum=gsr[0], maximum=gsr[1], title='Min', drag=1, xsize=500, uvalue='logsrangeminslider',/suppress_value)
logsrangemintext = WIDGET_TEXT(logsrangemin, /editable, uvalue='logsrangemintext', value=string(gsr[0],format='(e10.2)'),xsize=10, ysize=1)
(*dataptr).srange[0] = double(gsr[0])
(*dataptr).logsrangemin = logsrangemin
(*dataptr).logsrangeminslider = logsrangeminslider
(*dataptr).logsrangemintext = logsrangemintext

logsrangemax = WIDGET_BASE(srangecon,/row, uvalue='logsrangemax')
logsrangemaxslider = CW_FSLIDER(logsrangemax, value=gsr[1], minimum=gsr[0], maximum=gsr[1], title='Max', drag=1, xsize=500, uvalue='logsrangemaxslider',/suppress_value)
logsrangemaxtext = WIDGET_TEXT(logsrangemax, /editable, uvalue='logsrangemaxtext', value=string(gsr[1],format='(e10.2)'),xsize=10, ysize=1)
(*dataptr).srange[1] = gsr[1]
(*dataptr).logsrangemax = logsrangemax
(*dataptr).logsrangemaxslider = logsrangemaxslider
(*dataptr).logsrangemaxtext = logsrangemaxtext

(*dataptr).contours = 4
contoursbase = WIDGET_BASE(srangecon,/row,uvalue='contoursbase')
contourslabel = WIDGET_LABEL(contoursbase, value='Number of Contours')
contourstext = WIDGET_TEXT(contoursbase,/editable, uvalue='contourstext', value=string((*dataptr).contours), xsize=10, ysize=1)

(*dataptr).dolinear = 0
;dolinear = WIDGET_BUTTON(contoursbase, value='View Linear',uvalue='dolinear')
dolinear = WIDGET_BUTTON(contoursbase, value='images/loglinear.bmp',/bitmap, uvalue='dolinear')

boxmin = (*(*dataptr).treeptr).boxmin
boxmax = (*(*dataptr).treeptr).boxmax
(*dataptr).angles = [3*!PI/4.0,1*!PI/4.0,0*!PI/4.0]
(*dataptr).angles = [0,0,0]
target = (boxmax+boxmin)/2.0
(*dataptr).target = target
geocontrol = widget_base(/column,title='Slice Geometry',uvalue=dataptr, group_leader=*baseptr)
gsw=320
lrange = (*(*dataptr).treeptr).gridspacing[0,0] / (*(*dataptr).treeptr).gridspacing[0,(*(*dataptr).treeptr).maxlevel]
;print,lrange
lrange=alog10(lrange*4.0)
;print,lrange
zoomrange = [-1.0,lrange]
;print,zoomrange
(*dataptr).zoomrange = zoomrange
logzoom = WIDGET_BASE(geocontrol,/row, uvalue='logzoom')
logzoomslider = CW_FSLIDER(logzoom, value=0.0, uvalue='logzoomslider', minimum=zoomrange[0], maximum = zoomrange[1], title='Log of Zoom', drag=1, xsize=gsw,/suppress_value)
;logzoomsliderptr = ptr_new(logzoomslider,/no_copy)
logzoomtext = WIDGET_TEXT(logzoom, /editable, uvalue='logzoomtext', value=string(1.0),xsize=10, ysize=1)
;logzoomtextptr = ptr_new(logzoomtext,/no_copy)
(*dataptr).logzoom = logzoom
(*dataptr).logzoomslider = logzoomslider
(*dataptr).logzoomtext = logzoomtext

xcenter = WIDGET_BASE(geocontrol,/row,uvalue='xcenter')
xcenterslider = CW_FSLIDER(xcenter, value=target[0], minimum=boxmin[0], maximum=boxmax[0], title='X-Center', drag=1, xsize=gsw, uvalue='xcenterslider',/suppress_value)
xcentertext = WIDGET_TEXT(xcenter,/editable,uvalue='xcentertext', value=string(target[0], format='(e10.2)') ,xsize=10,ysize=1)
(*dataptr).xcenter = xcenter
(*dataptr).xcenterslider = xcenterslider
(*dataptr).xcentertext = xcentertext

ycenter = WIDGET_BASE(geocontrol,/row,uvalue='ycenter')
ycenterslider = CW_FSLIDER(ycenter, value=target[1], minimum=boxmin[1], maximum=boxmax[1], title='Y-Center', drag=1, xsize=gsw, uvalue='ycenterslider',/suppress_value)
ycentertext = WIDGET_TEXT(ycenter,/editable,uvalue='ycentertext', value=string(target[1], format='(e10.2)') ,xsize=10,ysize=1)
(*dataptr).ycenter = ycenter
(*dataptr).ycenterslider = ycenterslider
(*dataptr).ycentertext = ycentertext

zcenter = WIDGET_BASE(geocontrol,/row,uvalue='zcenter')
zcenterslider = CW_FSLIDER(zcenter, value=target[2], minimum=boxmin[2], maximum=boxmax[2], title='Z-Center', drag=1, xsize=gsw, uvalue='zcenterslider',/suppress_value)
zcentertext = WIDGET_TEXT(zcenter,/editable,uvalue='zcentertext', value=string(target[2], format='(e10.2)') ,xsize=10,ysize=1)
(*dataptr).zcenter = zcenter
(*dataptr).zcenterslider = zcenterslider
(*dataptr).zcentertext = zcentertext

angles = (*dataptr).angles
af = (*dataptr).angform
anglesliders = WIDGET_BASE(geocontrol,/column, uvalue='anglesliders')
azimuth = WIDGET_BASE(anglesliders,/row, uvalue='azimuth')
azimuthslider = CW_FSLIDER(azimuth,  value=angles[0], uvalue='azimuthslider',  minimum=-1.0*!PI, maximum=!PI,     title='Azimuth',  drag=1, xsize=gsw,/suppress_value)
azimuthtext = WIDGET_TEXT(azimuth, /editable, uvalue='azimuthtext', value=string(angles[0], format=af), xsize=10,ysize=1)
altitude = WIDGET_BASE(anglesliders,/row, uvalue='altitude')
altitudeslider = CW_FSLIDER(altitude, value=angles[1], uvalue='altitudeslider', minimum=0.0,      maximum=2.0*!PI, title='Altitude', drag=1, xsize=gsw,/suppress_value)
altitudetext = WIDGET_TEXT(altitude, /editable, uvalue='altitudetext', value=string(angles[1], format=af), xsize=10,ysize=1)
pitch = WIDGET_BASE(anglesliders,/row,uvalue='pitch')
pitchslider = CW_FSLIDER(pitch,    value=angles[2], uvalue='pitchslider',    minimum=-1.0*!PI, maximum=!PI,     title='Pitch',    drag=1, xsize=gsw,/suppress_value)
pitchtext = WIDGET_TEXT(pitch, /editable, uvalue='pitchtext', value=string(angles[2], format=af), xsize=10,ysize=1)
(*dataptr).azimuthslider = azimuthslider
(*dataptr).azimuthtext = azimuthtext
(*dataptr).altitudeslider = altitudeslider
(*dataptr).altitudetext = altitudetext
(*dataptr).pitchslider = pitchslider
(*dataptr).pitchtext = pitchtext

print,'!d.window:',!d.window
drawx = (*dataptr).drawx
drawy = (*dataptr).drawy
draw_slice = WIDGET_DRAW(*viewptr, retain=2, xsize=drawx,ysize=drawy, uvalue='plot_area',/button_events,/motion_events,/viewport_events)
print,'!d.window:',!d.window
olddw = !d.window
WIDGET_CONTROL, *baseptr,  /REALIZE
WIDGET_CONTROL, srangecon, /realize
WIDGET_CONTROL, geocontrol,/realize
colors = !d.n_colors
if colors ge 1.6e+7 then begin
    (*dataptr).cdepth=24
endif else begin
    (*dataptr).cdepth=8
endelse
cdepth = (*dataptr).cdepth
mydw = !d.window
wset,olddw
(*dataptr).mydw = mydw
print,'!d.window,mydw,olddw',!d.window,mydw,olddw
;plot,indgen(10),indgen(10)^2
read_ppm,'images/testpattern.pnm',testpat
help,testpat
;testpat = reverse(congrid(testpat,3,drawx,drawy),3)
wset,mydw
if cdepth ge 24 then begin
    tv,testpat,0,0,0,true=1
endif else begin
    tv,testpat,0,0,0
endelse
;xyouts,0,0,'Rendering Slice',charsize=3.0

wset,mydw
print,'!d.window:',!d.window

boxmin = (*treeptr).boxmin
boxmax = (*treeptr).boxmax
;zoom = (*dataptr).zoom
xr = [ boxmin[0],boxmax[0] ] / 1
yr = [ boxmin[1],boxmax[1] ] / 1
(*dataptr).xr = xr
(*dataptr).yr = yr
;image = slice_any(*treeptr,(*dataptr).angles,xrange=xr,yrange=yr)
;(*dataptr).rawimageptr = ptr_new(image)
;wh = where(image gt 0.0)
;min = min(image[wh])
;max = max(image)
;(*dataptr).srange = [min,max]
;wh = where(image eq 0.0)
;if wh[0] ne -1 then image[wh] = min
;;tv_colorized,bytscl(image)
;;tv_colorized,bytscl(alog10(image), min=alog10(min), max=alog10(max))

;image = bytscl(alog10(image), min=alog10(min), max=alog10(max))
c = intarr(256,3)
backup=c
ctnum=13
tvlct,backup,/get
loadct,ctnum,/silent
tvlct,c,/get
;orgc = c
tvlct,backup
;l = 47
;c[110:110+l, 1] = c[62:62+l, 1]
;c[62:62+l, 1] = 0
;c = c/2+orgc/2
;imc = [[[image]],[[image]],[[image]]]
(*dataptr).ct = c
;for i=0,2 do imc[*,*,i] = c[image[*,*],i]
;tv,imc,0,0,0,true=3
;(*dataptr).imageptr = ptr_new(imc)
ReDrawSlice,dataptr
;WIDGET_CONTROL,(*dataptr).logsrangeminslider, SET_VALUE=alog10(min)
;WIDGET_CONTROL, (*dataptr).logsrangemintext, SET_VALUE=string(min,format='(e10.2)')
;WIDGET_CONTROL,(*dataptr).logsrangemaxslider, SET_VALUE=alog10(max)
;WIDGET_CONTROL, (*dataptr).logsrangemaxtext, SET_VALUE=string(max,format='(e10.2)')
wset,olddw

XMANAGER, 'Slice_Widget', *baseptr
XMANAGER, 'Slice_Widget', srangecon
XMANAGER, 'Slice_Widget', geocontrol
print,'!d.window:',!d.window

;controlbase = widget_base(/column,uvalue=dataptr, title='Slice control')
;control_hello = widget_button(controlbase, value='Hello', uvalue='Hello')
;widget_control,controlbase,/realize
;xmanager, 'Slice_Widget', controlbase

end

; Event handler for slice_widget and the control windows it spwans
pro slice_widget_event,ev
WIDGET_CONTROL, ev.id, GET_UVALUE=uvalue
WIDGET_CONTROL, ev.top, GET_UVALUE=ps
;help,ps,/structures
;help,(*ps),/structures
;help,ev
case uvalue of
    'Hello' : begin
        ;help,ev,/structures
        ;help,ps, /structures
        ;help,*ps,/structures
    end
    'Close' : begin
        ;wait,5
        ;help,ev,/structures
        WIDGET_CONTROL, ev.top, /DESTROY
    end
    'Write PPM' : begin
        ;print,'Not able to write PPM'
        filename = dialog_pickfile(/write,filter=['*.ppm'])
        i = *((*ps).imageptr)
        ; for some reason, PPMs were coming out upside-down
        tvlct,r,g,b,/get
        si = size(i)
        help,i
        ;image = intarr(3,si[1],si[2])
        image = transpose(i,[2,0,1])
        help,image
        ;image[0,*,*] = r[i]
        ;image[1,*,*] = g[i]
        ;image[2,*,*] = b[i]
        image = reverse(image,3)
        help,image
        write_ppm,filename, image
        image = tvrd(true=1)
        help,image
        tv,image
        write_ppm,filename,image
    end
    'Write JPEG' : begin
        ;print,'Not able to write JPEG'
        i = *((*ps).imageptr)
        tvlct,r,g,b,/get
        help,i
        image = i
        ;image = [[[i]],[[i]],[[i]]]
        ;image[*,*,0] = r[i]
        ;image[*,*,1] = g[i]
        ;image[*,*,2] = b[i]
        help,image
        filename = dialog_pickfile(/write,filter=['*.jpeg'])
        write_jpeg,filename, image,quality=100,true=3
    end
    'contourstext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=textvalue
        str = double(textvalue)
        str = str[0]
        str = fix(str)
        if str ge 0 then begin
            (*ps).contours = str
            ReDrawSlice,ps,norerender=1,keepsrange=1
        endif else begin
            WIDGET_CONTROL, ev.id, SET_VALUE='0'
        endelse
    end
    'dolinear' : begin
        if (*ps).dolinear eq 0 then begin
            (*ps).dolinear = 1
            ;WIDGET_CONTROL, ev.id, SET_VALUE='View Log   '
            WIDGET_CONTROL, ev.id, SET_VALUE='images/linearlog.bmp',/bitmap
            ReDrawSlice,ps,keepsrange=1,norerender=1
        endif else begin
            (*ps).dolinear = 0
            ;WIDGET_CONTROL, ev.id, SET_VALUE='View Linear'
            WIDGET_CONTROL, ev.id, SET_VALUE='images/loglinear.bmp',/bitmap
            ReDrawSlice,ps,keepsrange=1,norerender=1
        endelse
    end
    'logsrangeminslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=logmin
        min = double(10.0)^logmin
        (*ps).srange[0] = min        
        ReDrawSlice,ps,norerender=1,keepsrange=1
        WIDGET_CONTROL, (*ps).logsrangemintext, SET_VALUE=string(min,format='(e10.2)')
    end
    'logsrangemintext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=textvalue
        str = double(textvalue)
        ;help,textvalue
        ;help,str
        str = str[0]
        textvalue = textvalue[0]
        if (str ne 0.0) or (textvalue eq '0.0') or (textvalue eq '0') then begin
            ; if str is valid (conversion worked)
            min = str
            ;max = (*ps).srange[1]
            (*ps).srange[0] = min
            ReDrawSlice,ps,keepsrange=1,norerender=1
            logmin = alog10(min)            
            WIDGET_CONTROL, (*ps).logsrangeminslider, SET_VALUE=logmin
        endif
    end


    'logsrangemaxslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=logmax
        max = double(10.0)^logmax
        (*ps).srange[1] = max
        ReDrawSlice,ps,keepsrange=1,norerender=1
        WIDGET_CONTROL, (*ps).logsrangemaxtext, SET_VALUE=string(max,format='(e10.2)')
    end
    'logsrangemaxtext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=textvalue
        str = double(textvalue)
        ;help,textvalue
        ;help,str
        str = str[0]
        textvalue = textvalue[0]
        if (str ne 0.0) or (textvalue eq '0.0') or (textvalue eq '0') then begin
            ; if str is valid (conversion worked)
            max = str
            (*ps).srange[1] = max
            ReDrawSlice,ps,keepsrange=1,norerender=1
            logmax = alog10(max)            
            WIDGET_CONTROL, (*ps).logsrangemaxslider, SET_VALUE=logmax
        endif
    end

    'logzoomslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=sliderpos
        min = (*ps).zoomrange[0]
        max = (*ps).zoomrange[1]
        zoom = 10.0^(sliderpos)
        
        WIDGET_CONTROL, (*ps).logzoomtext, SET_VALUE=string(zoom)
        if ev.drag eq 999 then begin
            ; quick and ugly rendering
            if zoom ge (*ps).zoom then begin
                oldimage = *((*ps).rawimageptr)
                rz = zoom/(*ps).zoom
                size = (*ps).tv_imgsize
                zoomi = congrid( (*((*ps).rawimageptr))[size[0]/2-size[0]/2/rz:size[0]/2+size[0]/2/rz, size[1]/2-size[1]/2/rz:size[1]/2+size[1]/2/rz], size[0],size[1])
                (*ps).rawimageptr = ptr_new(zoomi,/no_copy)
                xrb = (*ps).xr
                yrb = (*ps).yr
                zoomb = (*ps).zoom
                cx = ((*ps).xr[0]+(*ps).xr[1])/2.0
                cy = ((*ps).yr[0]+(*ps).yr[1])/2.0
                center = [cx,cy]
                treeptr = (*ps).treeptr
                boxmin = (*treeptr).boxmin
                boxmax = (*treeptr).boxmax
                span = [ boxmax[0]-boxmin[0] , boxmax[1]-boxmin[1] ] / zoom
                xr = [center[0]-span[0]/2.0, center[0]+span[0]/2.0]
                yr = [center[1]-span[1]/2.0, center[1]+span[1]/2.0]
                (*ps).xr = xr
                (*ps).yr = yr
                (*ps).zoom = zoom
                ReDrawSlice,ps,norerender=1
                (*ps).rawimageptr = ptr_new(oldimage,/no_copy)
                (*ps).xr = xrb
                (*ps).yr = yrb
                (*ps).zoom = zoomb
            endif
        endif
        if ev.drag ne 1 then begin ; let go of slider... rerender
            if 0 then begin
                
            endif else begin
                cx = ((*ps).xr[0]+(*ps).xr[1])/2.0
                cy = ((*ps).yr[0]+(*ps).yr[1])/2.0
                center = [cx,cy]
                treeptr = (*ps).treeptr
                boxmin = (*treeptr).boxmin
                boxmax = (*treeptr).boxmax
                span = [ boxmax[0]-boxmin[0] , boxmax[1]-boxmin[1] ] / zoom
                xr = [center[0]-span[0]/2.0, center[0]+span[0]/2.0]
                yr = [center[1]-span[1]/2.0, center[1]+span[1]/2.0]
                (*ps).xr = xr
                (*ps).yr = yr
                (*ps).zoom = zoom
                ReDrawSlice,ps
            endelse
        endif
    end
    'logzoomtext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=textvalue
        str = double(textvalue)
        help,textvalue
        help,str
        str = str[0]
        if (str[0] gt 0.0) then begin
            if 0 then begin
                min = (*ps).zoomrange[0]
                max = (*ps).zoomrange[1]
                zoom = str
                (*ps).zoom = str
                sliderpos = alog10(zoom)
                                ;help,ps,*ps,(*ps).logzoomslider
                
                WIDGET_CONTROL, long((*ps).logzoomslider), SET_VALUE=sliderpos
                treeptr = (*ps).treeptr
                boxmin = (*treeptr).boxmin
                boxmax = (*treeptr).boxmax
                zoom = (*ps).zoom
                xr = [ boxmin[0],boxmax[0] ] / zoom
                yr = [ boxmin[1],boxmax[1] ] / zoom
                image = slice_any( *treeptr,(*ps).angles, xrange=xr, yrange=yr)
                ptr_free,(*ps).rawimageptr
                (*ps).rawimageptr = ptr_new(image)
                wh = where(image gt 0.0)
                min = min(image[wh])
                max = max(image)
                (*ps).srange = [min,max]
                WIDGET_CONTROL,(*ps).logsrangeminslider, SET_VALUE=alog10(min)
                WIDGET_CONTROL, (*ps).logsrangemintext, SET_VALUE=string(min,format='(e10.2)')
                WIDGET_CONTROL,(*ps).logsrangemaxslider, SET_VALUE=alog10(max)
                WIDGET_CONTROL, (*ps).logsrangemaxtext, SET_VALUE=string(max,format='(e10.2)')
                wh = where(image eq 0.0)
                if wh[0] ne -1 then image[wh] = min
                                ;tv_colorized,bytscl(image)
                image = bytscl(alog10(image))
                wset,(*ps).mydw
                tv_colorized,image
                ptr_free,(*ps).imageptr
                (*ps).imageptr = ptr_new(image)
            endif else begin
                zoom = str
                sliderpos = alog10(zoom)
                                ;help,ps,*ps,(*ps).logzoomslider            
                WIDGET_CONTROL, long((*ps).logzoomslider), SET_VALUE=sliderpos
                cx = ((*ps).xr[0]+(*ps).xr[1])/2.0
                cy = ((*ps).yr[0]+(*ps).yr[1])/2.0
                center = [cx,cy]
                treeptr = (*ps).treeptr
                boxmin = (*treeptr).boxmin
                boxmax = (*treeptr).boxmax
                span = [ boxmax[0]-boxmin[0] , boxmax[1]-boxmin[1] ] / zoom
                xr = [center[0]-span[0]/2.0, center[0]+span[0]/2.0]
                yr = [center[1]-span[1]/2.0, center[1]+span[1]/2.0]
                (*ps).xr = xr
                (*ps).yr = yr
                (*ps).zoom = zoom
                ReDrawSlice,ps,oversample=2
            endelse
        endif
    end
    'xcenterslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=xcenter
        WIDGET_CONTROL, (*ps).xcentertext, SET_VALUE=string(xcenter,format='(e10.2)')
        if ev.drag ne 1 then begin
            (*ps).target[0] = xcenter
            ReDrawSlice,ps
        endif
    end
    'xcentertext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=xcenterstring
        str = double(xcenterstring)
        str = str[0]
        WIDGET_CONTROL, (*ps).xcenterslider, SET_VALUE=str
        (*ps).target[0] = str
        ReDrawSlice,ps
    end
    'ycenterslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=ycenter
        WIDGET_CONTROL, (*ps).ycentertext, SET_VALUE=string(ycenter,format='(e10.2)')
        if ev.drag ne 1 then begin
            (*ps).target[1] = ycenter
            ReDrawSlice,ps
        endif
    end
    'ycentertext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=ycenterstring
        str = double(ycenterstring)
        str = str[0]
        WIDGET_CONTROL, (*ps).ycenterslider, SET_VALUE=str
        (*ps).target[1] = str
        ReDrawSlice,ps
    end
    'zcenterslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=zcenter
        WIDGET_CONTROL, (*ps).zcentertext, SET_VALUE=string(zcenter,format='(e10.2)')
        if ev.drag ne 1 then begin
            (*ps).target[2] = zcenter
            ReDrawSlice,ps
        endif
    end
    'zcentertext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=zcenterstring
        str = double(zcenterstring)
        str = str[0]
        WIDGET_CONTROL, (*ps).zcenterslider, SET_VALUE=str
        (*ps).target[2] = str
        ReDrawSlice,ps
    end
    'azimuthslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=azimuth
        astr = string(azimuth/!DTOR, format=(*ps).angform)
        WIDGET_CONTROL, (*ps).azimuthtext, set_value=astr
        if ev.drag eq 0 then begin
            (*ps).angles[0] = azimuth
            ReDrawSlice,ps
        endif
    end
    'azimuthtext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=str
        str = double(str)
        str = str[0]
        str = str*!DTOR
        if str gt 1.0*!PI then str = 1.0*!PI
        if str lt -1.0*!PI then str = -1.0*!PI
        WIDGET_CONTROL, (*ps).azimuthslider, set_value=str
        (*ps).angles[0] = str
        astr = string(str/!DTOR, format=(*ps).angform)
        WIDGET_CONTROL, (*ps).azimuthtext, set_value=astr
        ReDrawSlice,ps
    end
    'altitudeslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=altitude
        astr = string(altitude/!DTOR, format=(*ps).angform)
        WIDGET_CONTROL, (*ps).altitudetext, set_value=astr
        if ev.drag eq 0 then begin            
            (*ps).angles[1] = altitude            
            ReDrawSlice,ps
        endif
    end
    'altitudetext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=str
        str = double(str)
        str = str[0]
        str = str*!DTOR
        if str gt 2.0*!PI then str = 2.0*!PI
        if str lt -0.0*!PI then str = -0.0*!PI
        WIDGET_CONTROL, (*ps).altitudeslider, set_value=str
        (*ps).angles[1] = str
        astr = string(str/!DTOR, format=(*ps).angform)
        WIDGET_CONTROL, (*ps).altitudetext, set_value=astr
        ReDrawSlice,ps
    end
    'pitchslider' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=pitch
        astr = string(pitch/!DTOR, format=(*ps).angform)
        WIDGET_CONTROL, (*ps).pitchtext, set_value=astr
        if ev.drag eq 0 then begin
            (*ps).angles[2] = pitch
            ReDrawSlice,ps
        endif
    end
    'pitchtext' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=str
        str = double(str)
        str = str[0]
        str = str*!DTOR
        if str gt !PI then str = !PI
        if str lt -1.0*!PI then str = -1.0*!PI
        WIDGET_CONTROL, (*ps).pitchslider, set_value=str
        (*ps).angles[2] = str
        astr = string(str/!DTOR, format=(*ps).angform)
        WIDGET_CONTROL, (*ps).pitchtext, set_value=astr
        ReDrawSlice,ps
    end
    'plot_area' : begin
        if ev.clicks eq 2 then begin
            help,ev,/structures
        endif
        if ev.type eq 2 then begin
            ; the mouse was moved
            if total((*ps).but_press[0,*]) ne -2 then begin
                ; button 0 pressed
                print,(*ps).but_press, ev.x,ev.y
                drag = [ev.x,ev.y]-(*ps).but_press[0,*]
                ;tv_colorized,(*(*ps).imageptr),drag[0],drag[1]
                help,(*(*ps).imageptr)
                wset,(*ps).mydw
                tv,(*(*ps).imageptr),(*ps).tv_offset[0]*!d.x_size+drag[0],(*ps).tv_offset[1]*!d.y_size+drag[1],0,true=3
            endif
        endif
        ; is the 2nd bit set?
        if (ishft(ev.press,-2) mod 2) or (ishft(ev.release,-2) mod 2) then begin
            print,'button2'
            help,ev,/structures
            if (ishft(ev.press,-2) mod 2) and ev.clicks then begin
                (*ps).but_press[2,*] = [ev.x,ev.y]
            endif else if (ishft(ev.release,-2) mod 2) then begin
                (*ps).but_release[2,*] = [ev.x,ev.y]
                print,'Drag2XY:', (*ps).but_release[2,*]-(*ps).but_press[2,*]
                (*ps).but_press[2,*] = [-1,-1]
            endif
        endif
        if (ishft(ev.press,-1) mod 2) or (ishft(ev.release,-1) mod 2) then begin
            print,'button1'
            help,ev,/structures
            if (ishft(ev.press,-1) mod 2) and ev.clicks then begin
                (*ps).but_press[1,*] = [ev.x,ev.y]
            endif else if (ishft(ev.release,-1) mod 2) then begin
                (*ps).but_release[1,*] = [ev.x,ev.y]
                drag1xy = (*ps).but_release[1,*]-(*ps).but_press[1,*]
                print,'Drag1XY:', drag1xy
                print,'But1 move:',(*ps).but_press[1,*], (*ps).but_release[1,*]
                len = sqrt(total(double(drag1xy)*double(drag1xy)))
                len = ceil(len)
                print,'len=',len
                elev = dblarr(len)
                pos = dblarr(len)
                idx = double(drag1xy[0])/double(len)
                idy = double(drag1xy[1])/double(len)
                simxy = [(*ps).xr[1]-(*ps).xr[0],(*ps).yr[1]-(*ps).yr[0] ]*double(drag1xy) / double((*ps).tv_imgsize)
                simlen = sqrt(total(double(simxy)*double(simxy)))
                xt = 'cm'
                reduce = 1.0
                AU = 1.49598e+13
                Pc = 3.08568025e+18
                if simlen ge AU*1.0E-03 then begin
                    ; 1AU=1.49598e+13cm, but handle mAU too
                    xt = 'AU'
                    reduce = AU
                endif
                ;if simlen gt 9.4605284e+17 then begin
                ;    xt = 'Light-Years'
                ;    reduce = 9.4605284e+17
                ;endif
                if simlen ge AU*1.0E+04 then begin
                    ; 1Pc=3.08568025e+18cm, but handle mPc too
                    xt='Parsecs'
                    reduce = Pc
                endif                
                x = 0
                y = 0
                xy0 = (*ps).but_press[1,*]-(*ps).tv_offset*[!d.x_size,!d.y_size]
                idp = simlen/(double(len-1))
                for i=0,len-1 do begin
                    x = xy0[0]+i*idx
                    y = xy0[1]+i*idy                    
                    elev[i] = (*(*ps).rawimageptr)[x,y]
                    pos[i] = idp*i
                endfor
                print,'XY=',x,y
                window,retain=2
                plot,pos/reduce,elev,psym=10,xtitle=xt,xstyle=1
                (*ps).but_press[1,*] = [-1,-1]
                
            endif
        endif
        ; is the 0th bit (button) set?
        if (ishft(ev.press,-0) mod 2) or (ishft(ev.release,-0) mod 2) then begin
            ; primary mouse button pressed or released
            ;help,*ps,/structures
            if ev.clicks eq 2 then begin
                ; double click.  recenter on position
                ; convert x,y into sim space
                pixxyz = double([ev.x,ev.y,0])
                ;print,'Pxyz:', pixxyz
                ;pixxyz = pixxyz - [(*ps).drawx,(*ps).drawy, 0]/2
                ;rxyz = pixxyz / double([(*ps).drawx,(*ps).drawy, 1])
                ;print,'Rxyz:',rxyz
                ;rxyz =rxyz*[(*ps).xr[0]-(*ps).xr[1],(*ps).yr[0]-(*ps).yr[1]]
                pixxy = pixxyz[0:1] - (*ps).tv_offset*[!d.x_size,!d.y_size] - (*ps).tv_imgsize/2.0
                print,'Pixxy:',pixxy
                rxy = pixxy / (*ps).tv_imgsize
                print,'Rxy:', rxy
                rxy = rxy * [(*ps).xr[0]-(*ps).xr[1],(*ps).yr[0]-(*ps).yr[1]]
                (*ps).xr = (*ps).xr - rxy[0]
                (*ps).yr = (*ps).yr - rxy[1]
                ReDrawSlice,ps,/keepsrange
                (*ps).but_press[0,*] = [-1,-1]
            endif
            if ev.press eq 1 then (*ps).but_press[0,*] = [ev.x,ev.y]
            if ev.release eq 1 then begin
                (*ps).but_release[0,*] = [ev.x,ev.y]
                deltaxy = [(*ps).but_release[0,*]-(*ps).but_press[0,*]]
                if deltaxy[0] ne 0 and deltaxy[1] ne 0 then begin
                    print,'DeltaXY: ',deltaxy
                                ; convert dxdy into simulation space
                    pdxy = [(*ps).but_release[0,*]-(*ps).but_press[0,*]]
                    rdxy = double([(*ps).but_release[0,*]-(*ps).but_press[0,*]]) / double([(*ps).tv_imgsize[0],(*ps).tv_imgsize[1]])
                    sdxy = rdxy * [(*ps).xr[0]-(*ps).xr[1],(*ps).yr[0]-(*ps).yr[1]]
                    xr = (*ps).xr
                    yr = (*ps).yr
                                ; nr0 is the part of the old image that is reused
                                ; nr1 is the vertical strip outside
                                ;   the old image, plus the corner
                                ; nr2 is the horizontal strip outside the old image
                    ;if 0 then begin
                    ;    imgmax = (size( (*(*ps).rawimageptr)))[1:2] - [1,1]
                    ;    if sdxy[0] gt 0.0 then begin
                    ;        nr1x = [ xr[1], xr[1]+sdxy[0]]
                    ;        nr2x = [ xr[0]+sdxy[0], xr[1]+sdxy[0]]
                    ;        sourcex = [pdxy[0],imgmax[0]]
                    ;        destx = [0,imgmax[0]-pdxy[0]]
                    ;    endif else begin
                    ;        nr1x = [xr[0]+sdxy[0], xr[1]+sdxy[0]]
                    ;        nr2x = [xr[0], xr[1]+sdxy[0]]
                    ;        nr0x = [xr[0], xr[1]+sdxy[0]]
                    ;        sourcex = [0, imgmax[0]-pdxy[0]]
                    ;        destx = [pdxy[0], imgmax[0]]
                    ;    endelse
                    ;    if sdxy[1] gt 0.0 then begin
                    ;        nr1y = [yr[0]+sdxy[1]. yr[1]+dsxy[1]]
                    ;        nr2y = [yr[1], yr[1]+sdxy[1]]
                    ;        sourcey = [pdxy[1], imgmax[1]]
                    ;        desty = [0, imgmax[1]-pdxy[0]]
                    ;    endif else begin
                    ;        nr1y = [yr[0]+sdxy[1], yr[1]+sdxy[1]]
                    ;        nr2y = [yr[0]+sdxy[1], yr[0]]
                    ;        sourcey = [0, imgmax[1]-pdxy[1]]
                    ;        desty = [pdxy[1], imgmax[1]]
                    ;    endelse
                    ;    treeptr = (*ps).treeptr
                    ;    nr1 = slice_any( *treeptr,(*ps).angles, xrange=nr1x, yrange=nr1y, imgdim = [abs(pdxy[0]), imgmax[1]+1])
                    ;    nr2 = slice_any( *treeptr,(*ps).angles, xrange=nr2x, yrange=nr2y, imgdim = [imgmax[0]+1-abs(pdxy[0]), abs(pdxy[1])] )
                    ;    help,nr1
                    ;    help,nr2
                    ;    composite = *((*ps).rawimageptr)
                    ;    composite[destx[0]:destx[1], desty[0]:desty[1]] = composite[sourcex[0]:sourcex[1],sourcey[0]:sourcey[1]]
                    ;    if (destx[0] eq 0) and (desty[0] eq 0) then begin
                    ;        composite[destx[1]+1:imgmax[0], 0:imgmax[1]] = nr1
                    ;        composite[0:destx[1], desty[1]+1:imgmax[1]] = nr2
                    ;        plot,composite
                    ;        wait,3
                    ;    endif
                    ;    if (destx[0] eq 0) and (desty[0] ne 0) then begin
                    ;    endif
                    ;    if (destx[0] ne 0) and (desty[0] eq 0) then begin
                    ;    endif
                    ;    if (destx[0] ne 0) and (desty[0] ne 0) then begin
                    ;    endif
                    ;endif
                    
                    (*ps).xr = (*ps).xr + sdxy[0]
                    (*ps).yr = (*ps).yr + sdxy[1]
                    ;(*ps).but_press[0,*] = [-1,-1]
                    ReDrawSlice,ps,/keepsrange
                endif
                (*ps).but_press[0,*] = [-1,-1]
            endif
        endif
    end
    else : begin
        print,'Slice_Widget_Event could not handle unknown event '+string(uvalue)
    end
endcase
end

; Renders a slice of a AMR Tree at arbitrary angle, using the settings
; in the data from an instance of SLICE_WIDGET
; INPUTS:
; p: pointer to a data structure of type SliceWidget.  The fields in
; this data determine how the slice is to be rendered
; oversample: Optional int that specifies that the slice is to be
; rendered into a higher resolution than *p calls for, and then
; interpolated down to the resolution requested by *p
pro ReRenderSlice,p,oversample=oversample
; p = pointer to WIDGET_SLICE uvalue structure
treeptr = (*p).treeptr
boxmin = (*treeptr).boxmin
boxmax = (*treeptr).boxmax
zoom = (*p).zoom
xr = (*p).xr
yr = (*p).yr
target = (*p).target
print,xr,yr
angles = (*p).angles
;wset,(*p).mydw
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
os=1
if keyword_set(oversample) then begin
    os = fix(oversample)
    os = max([os,1])
    os = min([os,10])
endif
(*p).tv_imgsize = [boxsize,boxsize]
help,boxsize
image = slice_any( *treeptr,angles, xrange=xr, yrange=yr,target=target, imgdim=[boxsize*os,boxsize*os])
help,image
image = rebin(image,boxsize,boxsize)
help,image
;print,angles,xr,yr,target,boxsize
rip = (*p).rawimageptr
;help,rip
ptr_free,rip
rip = ptr_new(image,/no_copy)
;help,rip
(*p).rawimageptr = rip
;help,(*p).rawimageptr
end

; Draws a slice of an AMR Tree for the SLICE_WIDGET program.  Actual
; rendering is done by calling RERENDERSLICE, which calls SLICE_ANY
; INPUTS:
; p: Pointer to the SliceWidget data structure of the instance of
; SLICE_WIDGET that called this function
; keepsrange: Set to keep the same scalar range set in *p  Default is
; to reset the scalar range based on the min/max of the new slice
; norerender: Set to prevent rerendering the slice.  The last
; rendering of the slice is used (much faster, good for tweaking srange)
; Oversample: Render into higher-resolution and then interpolate back
; to desired resolution.  Pretty but slow.
pro ReDrawSlice,p,keepsrange=keepsrange,norerender=norerender,oversample=oversample
if not keyword_set(norerender) then begin
    if keyword_set(oversample) then begin
        ReRenderSlice,p,oversample=oversample
    endif else begin
        ReRenderSlice,p
    endelse
endif
;treeptr = (*p).treeptr
;boxmin = (*treeptr).boxmin
;boxmax = (*treeptr).boxmax
zoom = (*p).zoom
xr = (*p).xr
yr = (*p).yr
target = (*p).target
print,xr,yr
angles = (*p).angles
;help,p,/structures
;help,*p,/structures
;help,(*p).rawimageptr,/structures
imp = (*p).rawimageptr
image = *imp

if not keyword_set(keepsrange) then begin
    wh = where(image gt 0.0)
    min = 0.0
    max = 0.0
    if wh[0] ne -1 then min = min(image[wh])
    max = max(image)
    (*p).srange = [min,max]
    WIDGET_CONTROL,(*p).logsrangeminslider, SET_VALUE=alog10(min)
    WIDGET_CONTROL, (*p).logsrangemintext, SET_VALUE=string(min,format='(e10.2)')
    WIDGET_CONTROL,(*p).logsrangemaxslider, SET_VALUE=alog10(max)
    WIDGET_CONTROL, (*p).logsrangemaxtext, SET_VALUE=string(max,format='(e10.2)')
endif else begin
    min = (*p).srange[0]
    max = (*p).srange[1]
endelse
print,'MinMax',[min,max]

wh = where(image lt min)
if wh[0] ne -1 then image[wh] = min
wh = where(image gt max)
if wh[0] ne -1 then image[wh] = max
print,'(*p).dolinear',(*p).dolinear
if (*p).dolinear eq 1 then begin
    image = bytscl(image, min=min, max=max)
endif else begin
    image = bytscl(alog10(image),min=alog10(min), max=alog10(max))
endelse

; apply a contour line at the 1/2 way point
;imhi = image gt 127
;cont = where(imhi - smooth(imhi,3))
;help,cont
;image[cont] = image[cont] / 2
cbar_div = 6
help,(*p).contours
if (*p).contours gt 0 then begin
    ;image = smooth(image,3)
    interval = 256.0 / ((*p).contours+1)
    steps = fix(image / interval)
    ;tvscl,steps
    ;wait,3
    cont = where(steps - smooth(steps,3))
    help,cont
    cbar_div = (*p).contours+1
endif

;image = bytscl(alog10(image), min=alog10(min), max=alog10(max))
c = (*p).ct
cdepth=(*p).cdepth
true=3
if cdepth ge 24 then begin
    true=3
endif else begin 
    true=0
endelse
if ((*p).contours eq 0) then begin
    if cdepth ge 24 then begin
        imc = [[[image]],[[image]],[[image]]]
        for i=0,2 do imc[*,*,i] = c[image[*,*],i]
    endif else begin
        imc = image
    endelse
endif else begin
    if cdepth ge 24 then begin
        imgr = image
        imgg = image
        imgb = image
        imgr[*,*] = c[image,0]
        imgg[*,*] = c[image,1]
        imgb[*,*] = c[image,2]
        help,imgr
        imgr[cont] = imgr[cont]/2
        imgg[cont] = imgg[cont]/2
        imgb[cont] = imgb[cont]/2
                                ;imgr = byte(imgr)
                                ;imgg = byte(imgg)
                                ;imgb = byte(imgb)
        imc = [[[imgr]],[[imgg]],[[imgb]]]
        help,imc
    endif else begin
        imc = image
        imc[cont] = byte(image)+byte(128)
    endelse
endelse
;for i=0,2 do imc[cont,i] = imc[cont,i] / 2
wset,(*p).mydw
erase
tv,imc,0.08,0.1,0,true=true,/normal
(*p).tv_offset = [0.08,0.1]
AR = double(!d.x_size) / double(!d.y_size)
if AR ge 1.3 then begin
    plot, [0,0], [0,0], /nodata, xrange=xr,yrange=yr,xstyle=1,ystyle=1,position=[0.08,0.1,0.08+0.85/AR,0.1+0.85],/normal,/noerase
    ;boxsize = 0.85*double(!d.y_size)
endif else begin
    plot, [0,0], [0,0], /nodata, xrange=xr,yrange=yr,xstyle=1,ystyle=1,position=[0.08,0.1,0.08+0.60,0.1+0.60*AR],/normal,/noerase
    ;boxsize = 0.6*AR*double(!d.y_size)
endelse
cbmin = min
cbmax = max
if (*p).dolinear eq 0 then begin
    cbmin = alog10(min)
    cbmax = alog10(max)
endif
if (abs(cbmax) lt 1.0e2) and (abs(cbmax) gt 1.0e-1) then $
  format = '(f9.2)' $
else $
  format = '(e10.2)'
if AR ge 1.3 then begin
    colorbar,position=[0.83,0.4,0.95,0.95],range=[cbmin,cbmax],/vertical,format=format,divisions=cbar_div
endif else begin
    colorbar,position=[0.83,0.1+0.3*AR,0.95,0.1+0.6*AR],range=[cbmin,cbmax],/vertical,format=format,divisions=cbar_div
endelse

; draw box for little axes
absize=0.125
abcenter = [0.83+absize/2.0,0.1+absize*AR/2.0,0.0]
;plot, [0,0], [0,0], /nodata,xrange=[-1.0*(xr[1]-xr[0])/8.0,(xr[1]-xr[0])/8.0], yrange=[-1.0*(yr[1]-yr[0])/8.0,(yr[1]-yr[0])/8.0],position=[0.83,0.1,0.83+absize,0.1+absize*AR],/normal,xticks=2,yticks=2,xstyle=1,ystyle=1,/noerase
;!p.multi=pbang

azimuth = (*p).angles[0]
altitude = (*p).angles[1]
pitch = (*p).angles[2]
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
aiv=0.75 ;arrow intensity variation, 1.0:(black-white), 0.0:(white-white)
xin = 255.0*((1.0-aiv)+aiv*(xarrow[2,7]+1.0)/2.0)
yin = 255.0*((1.0-aiv)+aiv*(yarrow[2,7]+1.0)/2.0)
zin = 255.0*((1.0-aiv)+aiv*(zarrow[2,7]+1.0)/2.0)
print,xin,yin,zin
xin = long(xin)
yin = long(yin)
zin = long(zin)
if cdepth ge 24 then begin
    xcolor=xin+xin*256.0+xin*65536.0
    ycolor=yin+yin*256.0+yin*65536.0
    zcolor=zin+zin*256.0+zin*65536.0
    color1 = 255.0 + 255.0*256.0
    color2 = 255.0*256.0 + 255.0*65536.0
    color3 = 255.0*65536.0 + 255.0
endif else begin
    xcolor = !p.color
    ycolor = !p.color-!p.color/3
    zcolor = !p.color-2*!p.color/3
    color1 = xcolor
    color2 = ycolor
    color3 = zcolor
endelse

for arrowi=0,7 do begin
    xarrow[*,arrowi] = QfixAR##xarrow[*,arrowi]
    yarrow[*,arrowi] = QfixAR##yarrow[*,arrowi]
    zarrow[*,arrowi] = QfixAR##zarrow[*,arrowi]
endfor


xarrow = xarrow*axisboxsize/2.0 + axisboxcenter#(dblarr(8)+1)
yarrow = yarrow*axisboxsize/2.0 + axisboxcenter#(dblarr(8)+1)
zarrow = zarrow*axisboxsize/2.0 + axisboxcenter#(dblarr(8)+1)
;xarrow = xarrow[0:1,*]
;yarrow = yarrow[0:1,*]
;zarrow = zarrow[0:1,*]



plots,xarrow[0,*],xarrow[1,*],/normal,color=xcolor
plots,yarrow[0,*],yarrow[1,*],/normal,color=ycolor
plots,zarrow[0,*],zarrow[1,*],/normal,color=zcolor
xyouts,xarrow[0,1],xarrow[1,1],'X',/normal,charsize=1.5,color=color1
xyouts,yarrow[0,1],yarrow[1,1],'Y',/normal,charsize=1.5,color=color2
xyouts,zarrow[0,1],zarrow[1,1],'Z',/normal,charsize=1.5,color=color3

                                ;tv_colorized,image
ptr_free,(*p).imageptr
help,imc
(*p).imageptr = ptr_new(imc)

end

; GUI program/wrapper for column-density tool.  Also includes a
; quick-annd dirty voxel-based mode.  Should be remade in the stype of
; SLICE_WIDGET
; INPUTS:
; treeptr: pointer to an AMR Tree structure
; group_leader: optional widget id of a parent window
; Use so that when the parent dies, this instance does too.
pro col_any_widget,treeptr,group_leader=group_leader
;help,*treeptr
x = treeptr
data = {ColAnyWidget, treeptr:treeptr, voxptr:ptr_new(), angles:dblarr(3), target:dblarr(3), zoom:double(1.0), log:1, mydw:fix(0), pixels:[512,512]}
;dataptr = ptr_new(data,/no_copy)

if keyword_set(group_leader) then begin
    base = WIDGET_BASE(/COLUMN, group_leader=group_leader)
endif else begin
    base = WIDGET_BASE(/COLUMN)
endelse
control = WIDGET_BASE(base,/row, uvalue='control')
angles = WIDGET_BASE(base,/column, uvalue='angles')
;otherbase = WIDGET_BASE(/COLUMN, uvalue=fltarr(50,50,50))
draw_voxel = WIDGET_DRAW(base, retain=2, xsize=data.pixels[0],ysize=data.pixels[1], uvalue='draw_voxel',/button_events,/motion_events,/viewport_events)
mydw = !d.window
data.mydw = mydw
print,'mydw',mydw
close = WIDGET_BUTTON(control, VALUE='Close', UVALUE='Close')
recalc_voxel = WIDGET_BUTTON(control, VALUE='ReCalc Voxel', UVALUE='ReCalc Voxel')
render_screen = WIDGET_BUTTON(control, VALUE='Render', UVALUE='Render')

sl_azimuth = CW_FSLIDER(angles, value=0.0, uvalue='azimuth', minimum=-1.0*!PI, maximum = !PI, title='Azimuth', drag=1, xsize=500)
sl_altiture = CW_FSLIDER(angles, value=0.0, uvalue='altitude', minimum=0.0, maximum=2.0*!PI, title='Altitude', drag=1, xsize=500)
sl_pitch = CW_FSLIDER(angles, value=0.0, uvalue='pitch', minimum=-1.0*!PI, maximum=!PI, title='Pitch', drag=1, xsize=500)

;tx_azimuth = WIDGET_TEXT(angles, event_pro='ReRenderVoxel', value=string(azimuth), xsize=10)
;tx_altitude = WIDGET_TEXT(angles, event_pro='ReRenderVoxel', value=string(altitude), xsize=10)
;tx_pitch = WIDGET_TEXT(angles, event_pro='ReRenderVoxel',
;value=string(pitch), xsize=10)



;s = {CAWStruct, vox:vox, size:128, altitude:double(0.0),azimuth:double(0.0),pitch:double(0.0)}
;s.vox[25:34, 55:74, 45:54] = 1.0
;s.vox[65:74, 55:74, 45:54] = 1.0
;s.vox[45:54, 25:74, 45:54] = 1.0
;s.vox[25:74, 45:54, 45:54] = 1.0

;ps = ptr_new(s,/no_copy)

olddw = !d.window
WIDGET_CONTROL, base, /REALIZE
mydw = !d.window
print,'mydw,olddw',mydw,olddw
data.mydw = mydw
;widget_control,otherbase,/realize
WIDGET_CONTROL,draw_voxel,GET_VALUE=windownumber
wset,windownumber
read_ppm,'images/testpattern.pnm',testpat
help,testpat
;testpat = reverse(congrid(testpat,3,drawx,drawy),3)
wset,mydw
tv,testpat,0,0,0,true=1

;help,base
;print,base
;print,'!D'
;print,!D
XMANAGER, 'Col_Any_Widget', base
;XMANAGER, 'Col_Any_Widget', otherbase
;ReRenderVoxel,*(data.voxptr),0,0,0,8,[!d.x_size, !d.y_size]
tp = data.treeptr
tps = *(data.treeptr)
ml = tps.maxlevel
xr = [tps.boxmin[0],tps.boxmax[0]]
yr = [tps.boxmin[1],tps.boxmax[1]]
zr = [tps.boxmin[2],tps.boxmax[2]]
vox = column_any(tps,[0,0,0],log=0,subsample=1.0,target=[0,0,0],maxlevel=ml,xrange=xr,yrange=yr,zrange=zr,imgdim=[128,128,128])
help,vox
data.voxptr = ptr_new(vox,/no_copy)

WIDGET_CONTROL, base, SET_UVALUE=ptr_new(data)
print,'*(data.voxptr) is:'
help,*(data.voxptr)
ReRenderVoxel,*(data.voxptr),0,0,0,4,[!d.x_size, !d.y_size],mydw
ReRenderVoxel,*(data.voxptr),0,0,0,2,[!d.x_size, !d.y_size],mydw
ReRenderVoxel,*(data.voxptr),0,0,0,1,[!d.x_size, !d.y_size],mydw
end


; Event handler for COL_ANY_WIDGET
pro col_any_widget_event,ev

;help,ev
;print,tag_names(ev)
;print,ev
;print,bork
WIDGET_CONTROL, ev.top, GET_UVALUE=dataptr
;help,ps
;help,*ps
help,dataptr,/structures
help,*dataptr,/structures
help,(*dataptr).voxptr
help,*((*dataptr).voxptr)
vox = *((*dataptr).voxptr)
WIDGET_CONTROL, ev.id, GET_UVALUE=uval
;print,uval
;widget_control, get_value=val
val = 0.0
imgsize = [!d.x_size, !d.y_size]
CASE uval OF
    'azimuth' : begin
        ;print,'did azimuth'
        WIDGET_CONTROL, ev.id, GET_VALUE=val
        ;print,val
        (*dataptr).angles[0] = val
        if ev.drag eq 1 then begin
            ; do a cheap-o rendering if mouse button not released
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],4,imgsize,(*dataptr).mydw
        endif else begin
                                ; do a slower rendering
                                ;ReRenderVoxel,vox,(*ps).azimuth,(*ps).altitude,(*ps).pitch,4,imgsize
            help,dataptr
            help,*dataptr
            help,(*dataptr).voxptr
            help,*((*dataptr).voxptr)
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],2,imgsize,(*dataptr).mydw
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],1,imgsize,(*dataptr).mydw
        endelse
    end
    'altitude' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=val
        (*dataptr).angles[1] = val
        if ev.drag eq 1 then begin
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],4,imgsize,(*dataptr).mydw
        endif else begin
            ;ReRenderVoxel,vox,(*ps).azimuth,(*ps).altitude,(*ps).pitch,4,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],2,imgsize,(*dataptr).mydw
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],1,imgsize,(*dataptr).mydw
        endelse
    end
    'pitch' : begin
        ;ReRenderVoxel
        WIDGET_CONTROL, ev.id, GET_VALUE=val
        (*dataptr).angles[2] = val
        if ev.drag eq 1 then begin
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],4,imgsize,(*dataptr).mydw
        endif else begin
            ;ReRenderVoxel,vox,(*ps).azimuth,(*ps).altitude,(*ps).pitch,4,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],2,imgsize,(*dataptr).mydw
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],1,imgsize,(*dataptr).mydw
        endelse
    end
    'draw_voxel' : begin
        print,'Event in window'
        ; button events are press or release
        ; movement events are mouse moving over the window
        ; vewport events are the window being dragged
        n = n_tags(ev)
        t = tag_names(ev)
        if (ishft(ev.press,-0) mod 2) or (ishft(ev.release,-0) mod 2) then begin
            ; primary mouse button
            if ev.clicks eq 2 then begin
                ; double-click
                ;pixxyz = double([ev.x,ev.y,0])
                pix = [ev.x,ev.y]
                rp = dblarr(3)
                ColWidget_pix_to_coord,dataptr,pix,rp
            endif
        endif
        ;if (ev.press ne 0) or (ev.release ne 0) then begin
            for i=0,n-1 do begin
                print,t[i]
            endfor
            print,ev
        ;endif
    end
    'Render' : begin
        ; make a window and re-render
        dw = !d.window
        print,!d.window
        ;window,/free,retain=2,xsize=640,ysize=480
        wset,(*dataptr).mydw
        col_any_frame,tree=(*(*dataptr).treeptr),angles=[(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2]],outM='X',log=1
        print,!d.window
        wset,dw
    end
  'Close': WIDGET_CONTROL, ev.top, /DESTROY
  else : begin
      print,'Col_Any_Widget_Event could not handle unknown event '+string(uval)
  end
ENDCASE
end

; Used so that a click in the COL_ANY_WIDGET window can be correlated
; to the coordinate system of the simulation
pro ColWidget_pix_to_coord,dp,pix,rp
; Convert pixel in a Col_Any_Widget draw window into a physical-unit
; coordinate
rp = pix
rp = rp - ((*dp).pixels)/2.0
rp = rp / (((*dp).pixels)/2.0)
help,dp
help,*dp
help,(*dp).treeptr
help,*(*dp).treeptr
rp = double([rp,0])
rp1 = double([rp,(*((*dp).treeptr)).boxmin[2]])
rp2 = double([rp,(*((*dp).treeptr)).boxmax[2]])
rp = rp * (  (*((*dp).treeptr)).boxmax - (*((*dp).treeptr)).boxmin)  / (*dp).zoom
rp1 = rp1 * (  (*((*dp).treeptr)).boxmax - (*((*dp).treeptr)).boxmin)  / (*dp).zoom
rp2 = rp2 * (  (*((*dp).treeptr)).boxmax - (*((*dp).treeptr)).boxmin)  / (*dp).zoom
print,rp
angles = (*dp).angles
azimuth=angles[0]
altitude = angles[1]
pitch = angles[2]
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
QROT = QR##(QP##(QA##QI))
QINV = invert(QROT,/double)
rp = QINV##rp
print,rp
print,'RPs:',rp1,rp2
;return,rp
end

; Used by COL_ANY_WIDGET to make a relatively fast but innacurate
; rendering of column density at an arbitrary angle
pro ReRenderVoxel,vox,azimuth,altitude,pitch,shrink,imgsize,mydw
print,'Called ReRenderVoxel'
sv = (size(vox))[1:3]
mss = sv[0]
if shrink ne 1 then begin
    voxt = congrid(vox,sv[0]/shrink, sv[1]/shrink, sv[2]/shrink)
endif else begin
    voxt = vox
endelse
miss = (machar()).xmin
wh = where(voxt gt 0.0)
if wh[0] ne -1 then miss = min(voxt[wh])
svt = (size(voxt))[1:3]
if azimuth ne 0.0 then begin
    for i=0,svt[2]-1 do voxt[*,*,i] = $
      rot(reform(voxt[*,*,i]),azimuth/!DTOR,/interp, cubic=cu, missing=miss)
endif
if altitude ne 0.0 then begin
    for i=0,svt[0]-1 do voxt[i,*,*] = $
      rot(reform(voxt[i,*,*]),altitude/!DTOR,/interp, cubic=cu, missing = miss)
endif
if pitch ne 0.0 then begin
    for i=0,svt[2]-1 do voxt[*,*,i] = $
      rot(reform(voxt[*,*,i]),pitch/!DTOR,/interp, cubic=cu, missing = miss)
endif
image = fltarr(svt[0],svt[1])
for i=0,svt[2]-1 do image = image + reform(voxt[*,*,i])
if (imgsize[0] ne svt[0]) or (imgsize[1] ne svt[1]) then begin
    wh = where(image gt 0.0)
    if wh[0] ne -1 then imin = min(image[wh])
    wh = where(image le 0.0)
    if wh[0] ne -1 then image[wh] = imin
    image = alog10(image)
    print,min(image),max(image)
    nm = (2.0*min(image) + max(image))/3.0
    image[where(image lt nm)] = nm
    image = congrid(image,imgsize[0], imgsize[1])
endif
wset,mydw
tv_colorized,bytscl(image)

end

;function cw_magicslider(parent, _extra=extra)
;data = {MagicSliderData, sliderpos:double(1.0), returnvalue:double(1.0), mode_ptr:ptr_new('exp10'), range:dblarr(2), cliprange:dblarr(2)}
;dataptr = ptr_new(data,/no_copy)
;base=widget_base(parent,/column,uvalue=dataptr,_extra=extra)
;slider = widget_slider(/drag,/suppress_value)
;text = widget_text(base)
;return,base
;end


; Program to show data along an arbitrary line.  Incomplete as of
; 9/30/2003
; INPUTS:
; amrptr: Pointer to an AMR object
; group_leader: Optional widget ID of a parent window.  Set so this
; dies when the parent dies.
pro line_widget,amrptr,group_leader=group_leader
;help,amrptr
;help,*amrptr
x = amrptr
data = {Line_Widget, $
        amrptr:amrptr, $
        ndim:(*amrptr).ndim, $
        drawx:fix(640), drawy:fix(480),mydw:fix(0), $
        srange:dblarr(2), $
        logsrangemin:long(0), logsrangeminslider:long(0), logsrangemintext:long(0), $
        logsrangemax:long(0), logsrangemaxslider:long(0), logsrangemaxtext:long(0), $
        contours:fix(0), $
        dolinear:fix(1), $
        p1:dblarr(3), p2:dblarr(3), $
        pf:'(e10.2)'}
dp = ptr_new(data,/no_copy)
ndim = (*dp).ndim
if ndim eq 1 then begin
(*dp).p1 = [(*amrptr).boxmin[0],0,0]
(*dp).p2 = [(*amrptr).boxmax[0],0,0]
endif
if ndim eq 2 then begin
(*dp).p1 = [(*amrptr).boxmin[0],(*amrptr).boxmin[1],0]
(*dp).p2 = [(*amrptr).boxmax[0],(*amrptr).boxmax[1],0]
endif
if ndim eq 3 then begin
(*dp).p1 = (*amrptr).boxmin
(*dp).p2 = (*amrptr).boxmax
endif
titlestr = 'Line Data on '+(*amrptr).name+':'+(*amrptr).componentName
if keyword_set(group_leader) then begin
    base = WIDGET_BASE(/column, uvalue=dp, title=titlestr, group_leader=group_leader)
endif else begin
    base = WIDGET_BASE(/column, uvalue=dp, title=titlestr)
endelse
draw = WIDGET_DRAW(base, retain=2, xsize=(*dp).drawx, ysize=(*dp).drawy, uvalue='plot_area',/button_events,/motion_events,/viewport_events)


scalarmin = min_amr(*amrptr)
scalarmax = max_amr(*amrptr)
gsr =[scalarmin,scalarmax]
gsr = alog10(gsr)
srangecon = WIDGET_BASE(/column, uvalue=dp, title='Line Scalar Range', group_leader=base)
logsrangemin = WIDGET_BASE(srangecon,/row, uvalue='logsrangemin')
logsrangeminslider = CW_FSLIDER(logsrangemin, value=gsr[0], minimum=gsr[0], maximum=gsr[1], title='Min', drag=1, xsize=500, uvalue='logsrangeminslider',/suppress_value)
logsrangemintext = WIDGET_TEXT(logsrangemin, /editable, uvalue='logsrangemintext', value=string(gsr[0],format='(e10.2)'),xsize=10, ysize=1)
(*dp).srange[0] = double(gsr[0])
(*dp).logsrangemin = logsrangemin
(*dp).logsrangeminslider = logsrangeminslider
(*dp).logsrangemintext = logsrangemintext

logsrangemax = WIDGET_BASE(srangecon,/row, uvalue='logsrangemax')
logsrangemaxslider = CW_FSLIDER(logsrangemax, value=gsr[1], minimum=gsr[0], maximum=gsr[1], title='Max', drag=1, xsize=500, uvalue='logsrangemaxslider',/suppress_value)
logsrangemaxtext = WIDGET_TEXT(logsrangemax, /editable, uvalue='logsrangemaxtext', value=string(gsr[1],format='(e10.2)'),xsize=10, ysize=1)
(*dp).srange[1] = gsr[1]
(*dp).logsrangemax = logsrangemax
(*dp).logsrangemaxslider = logsrangemaxslider
(*dp).logsrangemaxtext = logsrangemaxtext

(*dp).contours = 0
contoursbase = WIDGET_BASE(srangecon,/row,uvalue='contoursbase')
contourslabel = WIDGET_LABEL(contoursbase, value='Number of Contours')
contourstext = WIDGET_TEXT(contoursbase,/editable, uvalue='contourstext', value=string((*dp).contours), xsize=10, ysize=1)

(*dp).dolinear = 1
;dolinear = WIDGET_BUTTON(contoursbase, value='View Linear',uvalue='dolinear')
dolinear = WIDGET_BUTTON(contoursbase, value='images/loglinear.bmp',/bitmap, uvalue='dolinear')

geocontrol = WIDGET_BASE(/column, uvalue=dp, title='Line Geometry', group_leader=base)

point1 = WIDGET_BASE(geocontrol,/row, uvalue=dp)
point1_label = WIDGET_LABEL(point1, value='Point 1')
point1_set = WIDGET_BUTTON(point1, value='Set', uvalue='point1_set')
pf = (*dp).pf
if ndim ge 1 then begin
    x1_base = WIDGET_BASE(geocontrol,/row, uvalue=dp)
    x1_label = WIDGET_LABEL(x1_base, value='X')
    x1_input = WIDGET_TEXT(x1_base, /editable, uvalue='x1', value=string((*dp).p1[0], format=pf), xsize=10,ysize=1)
endif
if ndim ge 2 then begin
    y1_base = WIDGET_BASE(geocontrol,/row, uvalue=dp)
    y1_label = WIDGET_LABEL(y1_base, value='Y')
    y1_input = WIDGET_TEXT(y1_base, /editable, uvalue='y1', value=string((*dp).p1[1], format=pf), xsize=10,ysize=1)
endif
if ndim ge 3 then begin
    z1_base = WIDGET_BASE(geocontrol,/row, uvalue=dp)
    z1_label = WIDGET_LABEL(z1_base, value='Z')
    z1_input = WIDGET_TEXT(z1_base, /editable, uvalue='z1', value=string((*dp).p1[2], format=pf), xsize=10,ysize=1)
endif

point2 = WIDGET_BASE(geocontrol,/row, uvalue=dp)
point2_label = WIDGET_LABEL(point2, value='Point 2')
point2_set = WIDGET_BUTTON(point2, value='Set', uvalue='point2_set')
if ndim ge 1 then begin
    x2_base = WIDGET_BASE(geocontrol,/row, uvalue=dp)
    x2_label = WIDGET_LABEL(x2_base, value='X')
    x2_input = WIDGET_TEXT(x2_base, /editable, uvalue='x2', value=string((*dp).p2[0], format=pf), xsize=10,ysize=1)
endif
if ndim ge 2 then begin
    y2_base = WIDGET_BASE(geocontrol,/row, uvalue=dp)
    y2_label = WIDGET_LABEL(y2_base, value='Y')
    y2_input = WIDGET_TEXT(y2_base, /editable, uvalue='y2', value=string((*dp).p2[1], format=pf), xsize=10,ysize=1)
endif
if ndim ge 3 then begin
    z2_base = WIDGET_BASE(geocontrol,/row, uvalue=dp)
    z2_label = WIDGET_LABEL(z2_base, value='Z')
    z2_input = WIDGET_TEXT(z2_base, /editable, uvalue='z2', value=string((*dp).p2[2], format=pf), xsize=10,ysize=1)
endif





; Make X open the windows for this app

; Open the base window, and get its !d.window value so we can plot in
; the right place later
olddw = !d.window
WIDGET_CONTROL,base,/realize
mydw = !d.window
(*dp).mydw = mydw
wset,olddw
plot,[0,1],[2,3]

; Open the scalar control window
WIDGET_CONTROL,srangecon,/realize

WIDGET_CONTROL,geocontrol,/realize

; Bind each window to the event handler
help,base
help,srangecon
help,geocontrol
XMANAGER, 'line_widget',base,event_handler='line_widget_event'
XMANAGER, 'line_widget',srangecon,event_handler='line_widget_event'
XMANAGER, 'line_widget',geocontrol,event_handler='line_widget_event'

; End of line_widget
end


pro line_widget_event,ev
WIDGET_CONTROL, ev.id, GET_UVALUE=uvalue
WIDGET_CONTROL, ev.top, GET_UVALUE=ps
help,ev,/structures
case uvalue of
    'x1' : begin
        print,'Caught x1'
        help,ps
    end
    else: begin
        print,'Line_widget_event could not handle unknown event '+string(uvalue)
    end
endcase
; End of line_widget_event
end
