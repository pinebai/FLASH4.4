; this set of functions operates when the query button in xflash is
; pressed.  They allow you to select a point to find information about
; the data. 
;
; These are still under construction and not thouroughly tested with
; f77 formatted data
;

pro xquery, XVAL = x_zone, YVAL = y_zone, RVAL = r_zone, TVAL = theta_zone, $
            DXVAL = dx_zone, DYVAL = dy_zone, $
            REFINE = lrefine_zone, ZONE_INFO = zone_info, $
            VAR_NAMES = varnames, DIM = ndim, $
            SPHERICAL = spherical

if n_elements(spherical) EQ 0 then spherical = 0

; create the main widget base
base = widget_base(/column, title = 'xquery')

; divide this base vertically
base1 = widget_base(base, /column)
base2 = widget_base(base, /row)

; put the information onto the screen
nvar = (size(zone_info))[1]

label = indgen(nvar)

; print coordinate/resolution information

if (spherical) then begin

    r_label = widget_label(base1, /align_left, value = 'r: ' + $
                           strcompress(string(r_zone)) + ' cm')
   
    if (ndim EQ 2) then begin
        t_label  = widget_label(base1, /align_left, value = 'theta: ' + $
                               strcompress(string(theta_zone)) )
    endif

    blank1 = widget_label(base1, value = ' ')

    dr_label = widget_label(base1, /align_left, value = 'dr: ' + $
                           strcompress(string(dx_zone)) + ' cm')
    if (ndim EQ 2) then begin
        dt_label = widget_label(base1, /align_left, value = 'dt: ' + $
                               strcompress(string(dy_zone)) + ' deg')
    endif

endif else begin

    x_label  = widget_label(base1, /align_left, value = 'x: ' + $
                           strcompress(string(x_zone)) + ' cm')

    if (ndim EQ 2) then begin
        y_label  = widget_label(base1, /align_left, value = 'y: ' + $
                               strcompress(string(y_zone)) + ' cm')
    endif

    blank1 = widget_label(base1, value = ' ')

    dx_label = widget_label(base1, /align_left, value = 'dx: ' + $
                           strcompress(string(dx_zone)) + ' cm')

    if (ndim EQ 2) then begin
        dy_label = widget_label(base1, /align_left, value = 'dy: ' + $
                               strcompress(string(dy_zone)) + ' cm')
    endif

endelse

blank1 = widget_label(base1, value = ' ')

lrefine_label = widget_label(base1, /align_left, value = "refine level: " + $
                             strcompress(string(lrefine_zone)) )

blank1 = widget_label(base1, value = ' ')

; loop over the variables and print their value
for i = 0, nvar-1 do begin
    
    label[i] = widget_label(base1, /align_left, value = varnames[i] + ':' + $
                       strcompress(string(zone_info[i])) )

endfor

close  = widget_button(base2, value = 'close',  uvalue = 'close')

; draw the widget
widget_control, base, /realize
xmanager, 'xquery', base

end




;-----------------------------------------------------------------------------
; event handling routine for xquery -- just provide a mechanism for exiting
;-----------------------------------------------------------------------------
pro xquery_event, ev

widget_control, ev.top, get_uvalue=textwid
widget_control, ev.id, get_uvalue=uval

case uval of
    'close': widget_control, ev.top, /destroy
    else:      
endcase

end




pro query, X_CURS = x_curs, Y_CURS = y_curs, WIDGET_INFORMATION = info, $
           TREE_PTR = tree, PARAMS_PTR = params, DATA_PTR = data, $
           ORIENTATION = iswp, VAR_NAMES = varnames


print, 'in query, ', x_curs, y_curs

; get the limits of the current plot domain
plt_x_min = !x.crange[0]
plt_x_max = !x.crange[1]

plt_y_min = !y.crange[0]
plt_y_max = !y.crange[1]


iout = 0


case params.ndim of

 1: begin

; make sure the cursor was inside the plot domain
     if ((x_curs LT plt_x_min) OR (x_curs GT plt_x_max)) then begin

         widget_control, info.status, $
           set_value = 'status: outside of domain'

         iout = 1

     endif else begin
         
         block = (where(x_curs LE tree[*].bndBox[1,0] AND $
                        x_curs GT tree[*].bndBox[0,0] AND $
                        tree[*].nodetype EQ 1))[0]

; compute the minimum coord of each zone in the block
         x_block = (tree[block].bndBox[1,0] - tree[block].bndBox[0,0])* $
           findgen(params.nxb)/params.nxb + tree[block].bndBox[0,0]

; compute resolution
         dx_zone = (tree[block].bndBox[1,0] - tree[block].bndBox[0,0]) $
                  /params.nxb

; find the zone in the block which contains the cursor coords
         x_zone = max(where(x_curs GT x_block))
         y_zone = -1

         lrefine_zone = tree[block].lrefine

; in case we are in spherical geometry

         r_curs  = x_curs
         r_block = x_block
         dr_zone = dx_zone
         r_zone  = x_zone

; get the information in this zone
         zone_info = data[*,block,x_zone]
         
     endelse

 end

 2: begin

; make sure the cursor was inside the plot domain
     if ((iswp EQ 0 AND ((x_curs LT plt_x_min) OR (x_curs GT plt_x_max) OR $
                         (y_curs LT plt_y_min) OR (y_curs GT plt_y_max))) OR $
         (iswp EQ 1 AND ((y_curs LT plt_x_min) OR (y_curs GT plt_x_max) OR $
                         (x_curs LT plt_y_min) OR (x_curs GT plt_y_max))) OR  $
         (iswp EQ 2 AND ((x_curs LT plt_x_min) OR (x_curs GT plt_x_max) OR $
                         (y_curs GT plt_y_min) OR (y_curs LT plt_y_max)))) $
       then begin

         widget_control, info.status, $
           set_value = 'status: outside of domain'

         iout = 1

     endif else begin

; if we are in 2-d spherical coordinates, we must do a translation
; here
         if (params.geometry EQ "CARTESIAN" OR $
             params.geometry EQ "CYLINDRICAL") then begin

             block = (where(x_curs LE tree[*].bndBox[1,0] AND $
                            x_curs GT tree[*].bndBox[0,0] AND $
                            y_curs LE tree[*].bndBox[1,1] AND $
                            y_curs GT tree[*].bndBox[0,1] AND $
                            tree[*].nodetype EQ 1))[0]

; compute the minimum coord of each zone in the block
             x_block = (tree[block].bndBox[1,0] - tree[block].bndBox[0,0])* $
               findgen(params.nxb)/params.nxb + tree[block].bndBox[0,0]

             y_block = (tree[block].bndBox[1,1] - tree[block].bndBox[0,1])* $
               findgen(params.nyb)/params.nyb + tree[block].bndBox[0,1]

; compute resolution
             dx_zone = (tree[block].bndBox[1,0] - tree[block].bndBox[0,0]) $
                      /params.nxb

             dy_zone = (tree[block].bndBox[1,1] - tree[block].bndBox[0,1]) $
                      /params.nyb

; find the zone in the block which contains the cursor coords
             x_zone = max(where(x_curs GT x_block))
             y_zone = max(where(y_curs GT y_block))

             lrefine_zone = tree[block].lrefine

; get the information in this zone
             zone_info = data[*,block,x_zone,y_zone]
         
         endif else if (params.geometry EQ "SPHERICAL") then begin

             r_curs = sqrt(x_curs^2 + y_curs^2)
             theta_curs = atan(x_curs,y_curs)

             block = (where(r_curs LE tree[*].bndBox[1,0] AND $
                            r_curs GT tree[*].bndBox[0,0] AND $
                            theta_curs LE tree[*].bndBox[1,1] AND $
                            theta_curs GT tree[*].bndBox[0,1] AND $
                            tree[*].nodetype EQ 1))[0]
             
; compute the minimum coord of each zone in the block
             r_block = (tree[block].bndBox[1,0] - tree[block].bndBox[0,0])* $
               findgen(params.nxb)/params.nxb + tree[block].bndBox[0,0]

             t_block = (tree[block].bndBox[1,1] - tree[block].bndBox[0,1])* $
               findgen(params.nyb)/params.nyb + tree[block].bndBox[0,1]

; compute resolution
             dr_zone = (tree[block].bndBox[1,0] - tree[block].bndBox[0,0]) $
                      /params.nxb

             dt_zone = (tree[block].bndBox[1,1] - tree[block].bndBox[0,1]) $
                      /params.nyb * 180. / !pi

; find the zone in the block which contains the cursor coords
             r_zone = max(where(r_curs GT r_block))
             t_zone = max(where(theta_curs GT t_block))

             lrefine_zone = tree[block].lrefine

; get the information in this zone
             zone_info = data[*,block,r_zone,t_zone]
             
         endif else begin

             print, "Error: geometry not supported"
           
         endelse

     endelse

 end

endcase


print, 'about to pop up the widget'
if (iout EQ 0) then begin

    print, params.geometry

    if (params.geometry EQ "CARTESIAN" OR $
        params.geometry EQ "CYLINDRICAL") then begin
        xquery, XVAL = x_curs, YVAL = y_curs, DXVAL = dx_zone, DYVAL = dy_zone, REFINE = lrefine_zone, ZONE_INFO = zone_info, $
          VAR_NAMES = varnames, DIM = params.ndim
    endif else if (params.geometry EQ "SPHERICAL") then begin
        xquery, RVAL = r_curs, TVAL = theta_curs, $
                XVAL = x_curs, YVAL = y_curs, DXVAL = dr_zone, DYVAL = dt_zone, REFINE = lrefine_zone, ZONE_INFO = zone_info, $
          VAR_NAMES = varnames, DIM = params.ndim, /SPHERICAL
    endif

endif

end
