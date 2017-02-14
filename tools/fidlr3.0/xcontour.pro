; contour widget event handler
pro xcontour_event, ev

common save_contour_info, contours

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; take action pased on which widget was changed
case uval of

    'ctr1_plt': begin
        widget_control, info.ctr1_plot, get_value = select
        select = select[0]
        
        case select of 
            0: begin
                widget_control, info.ctr1_vars, sensitive = 0
                widget_control, info.ctr1_value, sensitive = 0
                widget_control, info.ctr1_color, sensitive = 0
            end
            1: begin
                widget_control, info.ctr1_vars, sensitive = 1
                widget_control, info.ctr1_value, sensitive = 1
                widget_control, info.ctr1_color, sensitive = 1
            end
        endcase

    end

    'ctr2_plt': begin
        widget_control, info.ctr2_plot, get_value = select
        select = select[0]

        case select of 
            0: begin
                widget_control, info.ctr2_vars, sensitive = 0
                widget_control, info.ctr2_value, sensitive = 0
                widget_control, info.ctr2_color, sensitive = 0
            end
            1: begin
                widget_control, info.ctr2_vars, sensitive = 1
                widget_control, info.ctr2_value, sensitive = 1
                widget_control, info.ctr2_color, sensitive = 1
            end
        endcase

    end

    'ctr3_plt': begin
        widget_control, info.ctr3_plot, get_value = select
        select = select[0]
        
        case select of 
            0: begin
                widget_control, info.ctr3_vars, sensitive = 0
                widget_control, info.ctr3_value, sensitive = 0
                widget_control, info.ctr3_color, sensitive = 0
            end
            1: begin
                widget_control, info.ctr3_vars, sensitive = 1
                widget_control, info.ctr3_value, sensitive = 1
                widget_control, info.ctr3_color, sensitive = 1
            end
        endcase

    end

    'ctr4_plt': begin
        widget_control, info.ctr4_plot, get_value = select
        select = select[0]
        
        case select of 
            0: begin
                widget_control, info.ctr4_vars, sensitive = 0
                widget_control, info.ctr4_value, sensitive = 0
                widget_control, info.ctr4_color, sensitive = 0
            end
            1: begin
                widget_control, info.ctr4_vars, sensitive = 1
                widget_control, info.ctr4_value, sensitive = 1
                widget_control, info.ctr4_color, sensitive = 1
            end
        endcase

    end


; *** done pressed ***
    'done': begin

        widget_control, info.ctr1_plot, get_value = status
        contours[0].enabled = status

        contours[0].var = widget_info(info.ctr1_vars, /droplist_select)

        widget_control, info.ctr1_value, get_value = value
        contours[0].value = value

        contours[0].color = widget_info(info.ctr1_color, /droplist_select)


        widget_control, info.ctr2_plot, get_value = status 
        contours[1].enabled = status

        contours[1].var = widget_info(info.ctr2_vars, /droplist_select)

        widget_control, info.ctr2_value, get_value = value
        contours[1].value = value

        contours[1].color = widget_info(info.ctr2_color, /droplist_select)


        widget_control, info.ctr3_plot, get_value = status
        contours[2].enabled = status

        contours[2].var = widget_info(info.ctr3_vars, /droplist_select)

        widget_control, info.ctr3_value, get_value = value
        contours[2].value = value

        contours[2].color = widget_info(info.ctr3_color, /droplist_select)


        widget_control, info.ctr4_plot, get_value = status
        contours[3].enabled = status

        contours[3].var = widget_info(info.ctr4_vars, /droplist_select)

        widget_control, info.ctr4_value, get_value = value
        contours[3].value = value

        contours[3].color = widget_info(info.ctr4_color, /droplist_select)

; store all the contour information in the contours structure

        widget_control, ev.top, /destroy
    end

    else:
endcase

end




pro xcontour, VARIABLES = variables

common save_contour_info, contours

; get the list of colors
idummy = color('black', get_list = color_list)

; create the main widget base
main_base = widget_base(/column, title = 'contour options')

; divide this base into vertical sections
information_base = widget_base(main_base, /row)
contour1_base    = widget_base(main_base, /row, /frame)
contour2_base    = widget_base(main_base, /row, /frame)
contour3_base    = widget_base(main_base, /row, /frame)
contour4_base    = widget_base(main_base, /row, /frame)
end_base         = widget_base(main_base, /column, /frame)

;-----------------------------------------------------------------------------
; contour # 1
;-----------------------------------------------------------------------------
ctr1_plot = cw_bgroup(contour1_base, ' ', label_left = '1: ', $
                      /nonexclusive, set_value = [contours[0].enabled], $
                      uvalue = 'ctr1_plt')

ctr1_vars = widget_droplist(contour1_base, title = 'Variable: ', $
                       uvalue = 'ctr1_var', value = variables)

; set the current variable
widget_control, ctr1_vars, set_droplist_select = contours[0].var

ctr1_value = cw_field(contour1_base, title = ' value: ', xsize = 8, $
                      uvalue = 'ctr1_val', value = contours[0].value)


ctr1_color = widget_droplist(contour1_base, title = 'Color: ', $
                             uvalue = 'ctr1_color', value = color_list)
widget_control, ctr1_color, set_droplist_select = contours[0].color



;-----------------------------------------------------------------------------
; contour # 2
;-----------------------------------------------------------------------------
ctr2_plot = cw_bgroup(contour2_base, ' ', label_left = '2: ', $
                      /nonexclusive, set_value = [contours[1].enabled], $
                      uvalue = 'ctr2_plt')

ctr2_vars = widget_droplist(contour2_base, title = 'Variable: ', $
                       uvalue = 'ctr2_var', value = variables)
widget_control, ctr2_vars, set_droplist_select = contours[1].var

ctr2_value = cw_field(contour2_base, title = ' value: ', xsize = 8, $
                      uvalue = 'ctr2_val', value = contours[1].value)

ctr2_color = widget_droplist(contour2_base, title = 'Color: ', $
                             uvalue = 'ctr2_color', value = color_list)
widget_control, ctr2_color, set_droplist_select = contours[1].color



;-----------------------------------------------------------------------------
; contour # 3
;-----------------------------------------------------------------------------
ctr3_plot = cw_bgroup(contour3_base, ' ', label_left = '3: ', $
                      /nonexclusive, set_value = [contours[2].enabled], $
                      uvalue = 'ctr3_plt')

ctr3_vars = widget_droplist(contour3_base, title = 'Variable: ', $
                       uvalue = 'ctr3_var', value = variables)
widget_control, ctr3_vars, set_droplist_select = contours[2].var

ctr3_value = cw_field(contour3_base, title = ' value: ', xsize = 8, $
                      uvalue = 'ctr3_val', value = contours[2].value)

ctr3_color = widget_droplist(contour3_base, title = 'Color: ', $
                             uvalue = 'ctr3_color', value = color_list)
widget_control, ctr3_color, set_droplist_select = contours[2].color



;-----------------------------------------------------------------------------
; contour # 4
;-----------------------------------------------------------------------------
ctr4_plot = cw_bgroup(contour4_base, ' ', label_left = '4: ', $
                      /nonexclusive, set_value = [contours[3].enabled], $
                      uvalue = 'ctr4_plt')

ctr4_vars = widget_droplist(contour4_base, title = 'Variable: ', $
                       uvalue = 'ctr4_var', value = variables)
widget_control, ctr4_vars, set_droplist_select = contours[3].var

ctr4_value = cw_field(contour4_base, title = ' value: ', xsize = 8, $
                      uvalue = 'ctr4_val', value = contours[3].value)

ctr4_color = widget_droplist(contour4_base, title = 'Color: ', $
                             uvalue = 'ctr4_color', value = color_list)
widget_control, ctr4_color, set_droplist_select = contours[3].color



if NOT contours[0].enabled then begin
    widget_control, ctr1_vars, sensitive = 0
    widget_control, ctr1_value, sensitive = 0
    widget_control, ctr1_color, sensitive = 0
endif

if NOT contours[1].enabled then begin
    widget_control, ctr2_vars, sensitive = 0
    widget_control, ctr2_value, sensitive = 0
    widget_control, ctr2_color, sensitive = 0
endif

if NOT contours[2].enabled then begin
    widget_control, ctr3_vars, sensitive = 0
    widget_control, ctr3_value, sensitive = 0
    widget_control, ctr3_color, sensitive = 0
endif

if NOT contours[3].enabled then begin
    widget_control, ctr4_vars, sensitive = 0
    widget_control, ctr4_value, sensitive = 0
    widget_control, ctr4_color, sensitive = 0
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
        contour1_base:contour1_base, $       ; contour # 1 base
        ctr1_plot:ctr1_plot, $               ; contour # 1 select
        ctr1_vars:ctr1_vars, $               ; contour # 1 variable
        ctr1_value:ctr1_value, $             ; contour # 1 value 
        ctr1_color:ctr1_color, $             ; contour # 1 color
        contour2_base:contour2_base, $       ; contour # 2 base
        ctr2_plot:ctr2_plot, $               ; contour # 2 select
        ctr2_vars:ctr2_vars, $               ; contour # 2 variable
        ctr2_value:ctr2_value, $             ; contour # 2 value
        ctr2_color:ctr2_color, $             ; contour # 2 color
        contour3_base:contour3_base, $       ; contour # 3 base
        ctr3_plot:ctr3_plot, $               ; contour # 3 select
        ctr3_vars:ctr3_vars, $               ; contour # 3 variable
        ctr3_value:ctr3_value, $             ; contour # 3 value
        ctr3_color:ctr3_color, $             ; contour # 3 color
        contour4_base:contour4_base, $       ; contour # 4 base
        ctr4_plot:ctr4_plot, $               ; contour # 4 select
        ctr4_vars:ctr4_vars, $               ; contour # 4 variable
        ctr4_value:ctr4_value, $             ; contour # 4 value
        ctr4_color:ctr4_color, $             ; contour # 4 color
        end_base:end_base, $                 ; end base widget
        done:done}                           ; button for done


; register the info structure with the widget base
widget_control, main_base, set_uvalue=info, /no_copy

xmanager, 'xcontour', main_base

end











