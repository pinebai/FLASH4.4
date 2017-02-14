; contour widget event handler
pro xhist_event, ev

common save_histogram_info, nbins, max_scale

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; take action pased on which widget was changed
case uval of

; *** done pressed ***
    'done': begin
        widget_control, info.nbins_field, get_value = temp
        nbins = temp[0]

        widget_control, info.scale_field, get_value = temp
        max_scale = temp[0]

        widget_control, ev.top, /destroy
    end

    else:
endcase

end




pro xhist

common save_histogram_info, nbins, max_scale

; create the main widget base
main_base = widget_base(/column, title = 'hist options')

; divide this base into vertical sections
nbins_base = widget_base(main_base, /column, /frame)
end_base   = widget_base(main_base, /column, /frame)

nbins_field = cw_field(nbins_base, title = ' number of bins: ', $
                       xsize = 12, uvalue = 'nbins', $
                       value = nbins)


scale_field = cw_field(nbins_base, title = ' vertical range: 0 to ', $
                       xsize = 12, uvalue = 'scale', $
                       value = max_scale)

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
        nbins_field:nbins_field, $
        scale_field:scale_field, $
        end_base:end_base, $                 ; end base widget
        done:done}                           ; button for done


; register the info structure with the widget base
widget_control, main_base, set_uvalue=info, /no_copy

xmanager, 'xhist', main_base

end











