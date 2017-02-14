; widget to render a quick low quality column-density plot

pro idlamrvis,filename=filename
initialmessage = 'Ready'
component = 'density'
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
center = (amr.boxmax+amr.boxmin)/2.0
(*dataptr).center = center
initialmessage = 'Making tree structure from plotfile'+filename
WIDGET_CONTROL,*messageptr, SET_VALUE=initialmessage

tree = amr_tree(amr,free=1)
treeptr = ptr_new(tree,/no_copy)
(*dataptr).treeptr = treeptr

initialmessage = 'Ready'
WIDGET_CONTROL,*messageptr, SET_VALUE=initialmessage

XMANAGER, 'IDLAMRVis', *baseptr
end

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
        col_any_widget,(*ps).treeptr
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Ready'
    end
    'Slice' : begin
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Creating slice...'
        slice_widget,(*ps).treeptr
        WIDGET_CONTROL, *((*ps).messageptr), SET_VALUE='Ready'
    end
endcase
end

pro slice_widget,treeptr
help,tree
data = {SliceWidget, treeptr:treeptr, xrange:dblarr(2), yrange:dblarr(2), viewaxis:dblarr(3), centerpoint:dblarr(3), srange:dblarr(2)}
base = WIDGET_BASE(/column)
control = WIDGET_BASE(base,/row, uvalue='control')
end

pro col_any_widget,treeptr
help,*treeptr
data = {ColAnyWidget, treeptr:treeptr, voxptr:ptr_new(), angles:dblarr(3), target:dblarr(3), zoom:double(1.0), log:1}
;dataptr = ptr_new(data,/no_copy)

base = WIDGET_BASE(/COLUMN)
control = WIDGET_BASE(base,/row, uvalue='control')
angles = WIDGET_BASE(base,/column, uvalue='angles')
;otherbase = WIDGET_BASE(/COLUMN, uvalue=fltarr(50,50,50))
draw_voxel = WIDGET_DRAW(base, retain=2, xsize=512,ysize=512, uvalue='draw_voxel',/button_events,/motion_events,/viewport_events)

close = WIDGET_BUTTON(control, VALUE='Close', UVALUE='Close')
render_screen = WIDGET_BUTTON(control, VALUE='Render', UVALUE='Render')

sl_azimuth = CW_FSLIDER(angles, value=0.0, uvalue='azimuth', minimum=-1.0*!PI, maximum = !PI, title='Azimuth', drag=1, xsize=500)
sl_altiture = CW_FSLIDER(angles, value=0.0, uvalue='altitude', minimum=0.0, maximum=2.0*!PI, title='Altitude', drag=1, xsize=500)
sl_pitch = CW_FSLIDER(angles, value=0.0, uvalue='pitch', minimum=-1.0*!PI, maximum=!PI, title='Pitch', drag=1, xsize=500)

;tx_azimuth = WIDGET_TEXT(angles, event_pro='ReRenderVoxel', value=string(azimuth), xsize=10)
;tx_altitude = WIDGET_TEXT(angles, event_pro='ReRenderVoxel', value=string(altitude), xsize=10)
;tx_pitch = WIDGET_TEXT(angles, event_pro='ReRenderVoxel',
;value=string(pitch), xsize=10)
tp = data.treeptr
tps = *(data.treeptr)
ml = tps.maxlevel
xr = [tps.boxmin[0],tps.boxmax[0]]
yr = [tps.boxmin[1],tps.boxmax[1]]
zr = [tps.boxmin[2],tps.boxmax[2]]
vox = column_any(tps,[0,0,0],log=0,subsample=1.0,target=[0,0,0],maxlevel=ml,xrange=xr,yrange=yr,zrange=zr,imgdim=[128,128,128])
help,vox
data.voxptr = ptr_new(vox,/no_copy)


;s = {CAWStruct, vox:vox, size:128, altitude:double(0.0),azimuth:double(0.0),pitch:double(0.0)}
;s.vox[25:34, 55:74, 45:54] = 1.0
;s.vox[65:74, 55:74, 45:54] = 1.0
;s.vox[45:54, 25:74, 45:54] = 1.0
;s.vox[25:74, 45:54, 45:54] = 1.0

;ps = ptr_new(s,/no_copy)

WIDGET_CONTROL, base, SET_UVALUE=ptr_new(data)
print,'*(data.voxptr) is:'
help,*(data.voxptr)
WIDGET_CONTROL, base, /REALIZE
;widget_control,otherbase,/realize
WIDGET_CONTROL,draw_voxel,GET_VALUE=windownumber
wset,windownumber

help,base
print,base
print,'!D'
print,!D
XMANAGER, 'Col_Any_Widget', base
;XMANAGER, 'Col_Any_Widget', otherbase
end

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
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],4,imgsize
        endif else begin
                                ; do a slower rendering
                                ;ReRenderVoxel,vox,(*ps).azimuth,(*ps).altitude,(*ps).pitch,4,imgsize
            help,dataptr
            help,*dataptr
            help,(*dataptr).voxptr
            help,*((*dataptr).voxptr)
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],2,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],1,imgsize
        endelse
    end
    'altitude' : begin
        WIDGET_CONTROL, ev.id, GET_VALUE=val
        (*dataptr).angles[1] = val
        if ev.drag eq 1 then begin
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],4,imgsize
        endif else begin
            ;ReRenderVoxel,vox,(*ps).azimuth,(*ps).altitude,(*ps).pitch,4,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],2,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],1,imgsize
        endelse
    end
    'pitch' : begin
        ;ReRenderVoxel
        WIDGET_CONTROL, ev.id, GET_VALUE=val
        (*dataptr).angles[2] = val
        if ev.drag eq 1 then begin
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],4,imgsize
        endif else begin
            ;ReRenderVoxel,vox,(*ps).azimuth,(*ps).altitude,(*ps).pitch,4,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],2,imgsize
            ReRenderVoxel,vox,(*dataptr).angles[0],(*dataptr).angles[1],(*dataptr).angles[2],1,imgsize
        endelse
    end
    'draw_voxel' : begin
        print,'Event in window'
        ; button events are press or release
        ; movement events are mouse moving over the window
        ; vewport events are the window being dragged
        n = n_tags(ev)
        t = tag_names(ev)
        ;if (ev.press ne 0) or (ev.release ne 0) then begin
            for i=0,n-1 do begin
                print,t[i]
            endfor
            print,ev
        ;endif
    end
    'Render' : begin
        ; make a window and re-render
    end
  'Close': WIDGET_CONTROL, ev.top, /DESTROY
ENDCASE
end


pro ReRenderVoxel,vox,azimuth,altitude,pitch,shrink,imgsize
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

tv_colorized,bytscl(image)

end
