pro xpick_geometry_event, ev
common geometry_common, geometryChoices, choiceId

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; take action pased on which widget was changed
case uval of

; *** done pressed ***
    'done': begin
        choiceId = widget_info(info.geometryChoiceDroplist, /droplist_select)
        widget_control, ev.top, /destroy
    end

    else:     
endcase
end 



function xpick_geometry, base
common geometry_common, geomteryChoices, choiceId

geometryChoices = [ "CARTESIAN", "CYLINDRICAL", "POLAR", "SPHERICAL" ] 

; create the main widget base
mainBase = widget_base(group_leader = base, /column, /modal, title = 'pick geometry')

; divide this base into vertical sections
choiceBase    = widget_base(mainBase, /column, /frame)
endBase         = widget_base(mainBase, /column, /frame)

geometryChoiceDroplist = $
  widget_droplist(choiceBase, title = 'Geometry: ', $
                  uvalue = 'geomChoice', value = geometryChoices)

widget_control, geometryChoiceDroplist, set_droplist_select = 0

;---------
; buttons
;---------
done = widget_button(endBase, value = 'Done', uvalue = 'done')
widget_control, mainBase, default_button = done

;-----------------
; draw the widget
;-----------------
widget_control, mainBase, /realize

; setup a structure to hold the widget information
info = {mainBase:mainBase, $               ; main base
        choiceBase:choiceBase, $
        geometryChoiceDroplist:geometryChoiceDroplist, $
        endBase:endBase, $                 ; end base widget
        done:done}                           ; button for done


; register the info structure with the widget base
widget_control, mainBase, set_uvalue=info, /no_copy

xmanager, 'xpick_geometry', mainBase, EVENT_HANDLER='xpick_geometry_event'

choice = geometryChoices[choiceId]
return, choice
end
  
