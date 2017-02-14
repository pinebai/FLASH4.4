; contour widget event handler
pro xvector_event, ev

common save_vector_info, vector

common widget_pass, variables

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; take action pased on which widget was changed
case uval of

  ; enabled
  'vec_plt': begin
    widget_control, info.vec_plot, get_value = enabled
    enabled = enabled[0]

    case enabled of
      0: begin
        widget_control, info.varx,    sensitive = 0
        widget_control, info.vary,    sensitive = 0
        widget_control, info.typ_vel, sensitive = 0
        widget_control, info.min_vel, sensitive = 0
        widget_control, info.max_vel, sensitive = 0
        widget_control, info.xskip,   sensitive = 0
        widget_control, info.yskip,   sensitive = 0
        widget_control, info.outline, sensitive = 0
      end
      1: begin
        widget_control, info.varx,    sensitive = 1
        widget_control, info.vary,    sensitive = 1
        widget_control, info.typ_vel, sensitive = 1
        widget_control, info.min_vel, sensitive = 1
        widget_control, info.max_vel, sensitive = 1
        widget_control, info.xskip,   sensitive = 1
        widget_control, info.yskip,   sensitive = 1
        widget_control, info.outline, sensitive = 1
      end
    endcase
  end

  ; *** done pressed ***
  'done': begin
    widget_control, info.vec_plot, get_value = enabled
    vector.enabled = enabled

    ivar = widget_info(info.varx, /droplist_select)
    vector.xcomp = variables(ivar)

    ivar = widget_info(info.vary, /droplist_select)
    vector.ycomp = variables(ivar)

    widget_control, info.typ_vel, get_value = temp
    vector.typical_vector = temp

    widget_control, info.min_vel, get_value = temp
    vector.min_vector = temp

    widget_control, info.max_vel, get_value = temp
    vector.max_vector = temp

    widget_control, info.xskip, get_value = temp
    vector.xskip = temp

    widget_control, info.yskip, get_value = temp
    vector.yskip = temp

    widget_control, info.outline, get_value = temp
    vector.outline = temp

    widget_control, ev.top, /destroy

  end

  else:
endcase

end



pro xvector, VARIABLES=varnames, NATIVE_NUMBER = native_num

common save_vector_info, vector

common widget_pass, variables


; we need to know the available variables in the event routine in
; order to map the droplist index back to the name
variables =  varnames

; create the main widget base
main_base = widget_base(/column, title = 'vector options')

; divide this base into vertical sections
enable_base      = widget_base(main_base, /column, /frame)
variable_base    = widget_base(main_base, /column, /frame)
information_base = widget_base(main_base, /column)
skip_base        = widget_base(main_base, /row, /frame)
option_base      = widget_base(main_base, /row, /frame)
end_base         = widget_base(main_base, /column, /frame)


vec_plot = cw_bgroup(enable_base, 'plot vectors', $
                     /nonexclusive, set_value = [vector.enabled], $
                     uvalue = 'vec_plt')

; only the atomic variables are allowed to be components
varx = widget_droplist(variable_base, title = 'x-component: ', $
                       uvalue = 'varx', value = varnames(0:native_num-1))

widget_control, varx, set_droplist_select = var_index(vector.xcomp)


vary = widget_droplist(variable_base, title = 'y-component: ', $
                       uvalue = 'vary', value = varnames(0:native_num-1))

widget_control, vary, set_droplist_select = var_index(vector.ycomp)


typ_vel = cw_field(information_base, title = ' typical vector: ', $
                   xsize = 12, uvalue = 'typ_vel', $
                   value = vector.typical_vector)

min_vel = cw_field(information_base, title = ' minimum vector: ', $
                   xsize = 12, uvalue = 'min_vel', $
                   value = vector.min_vector)


max_vel = cw_field(information_base, title = ' maximum vector: ', $
                   xsize = 12, uvalue = 'max_vel', $
                   value = vector.max_vector)

xskip = cw_field(information_base, title = ' x skip: ', $
                   xsize = 8, uvalue = 'xskip', $
                   value = vector.xskip)

yskip = cw_field(information_base, title = ' y skip: ', $
                   xsize = 8, uvalue = 'yskip', $
                   value = vector.yskip)

outline = cw_bgroup(option_base, 'draw outlines', $
                     /nonexclusive, set_value = [vector.outline], $
                     uvalue = 'vec_outline')

if NOT vector.enabled then begin
  widget_control, varx,    sensitive = 0
  widget_control, vary,    sensitive = 0
  widget_control, typ_vel, sensitive = 0
  widget_control, min_vel, sensitive = 0
  widget_control, max_vel, sensitive = 0
  widget_control, xskip,   sensitive = 0
  widget_control, yskip,   sensitive = 0
  widget_control, outline, sensitive = 0
endif


;---------
; buttons
;---------
done = widget_button(end_base, value = 'Done', uvalue = 'done')


;-----------------
; draw the widget
;-----------------
widget_control, main_base, /realize


; setup a structure to hold the widget information
info = {main_base:main_base, $               ; main base
        information_base:information_base, $ ; information base
        skip_base:skip_base, $               ; skip base
        vec_plot:vec_plot, $
        varx:varx, $
        vary:vary, $
        typ_vel:typ_vel, $
        min_vel:min_vel, $
        max_vel:max_vel, $
        xskip:xskip, $
        yskip:yskip, $
        outline:outline, $
        end_base:end_base, $                 ; end base widget
        done:done}                           ; button for done


; register the info structure with the widget base
widget_control, main_base, set_uvalue=info, /no_copy

xmanager, 'xvector', main_base

end
