; contour widget event handler
pro xlabel_event, ev

common save_label, label_info

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; take action pased on which widget was changed
case uval of


; enabled
    'flt_lbl': begin
        widget_control, info.use_label, get_value = temp
        enabled = temp[0]

        case enabled of
            0: begin
                widget_control, info.text_label, sensitive = 0
                widget_control, info.text_size, sensitive = 0
                widget_control, info.text_thick, sensitive = 0
                widget_control, info.minx, sensitive = 0
                widget_control, info.posx_lbl, sensitive = 0
                widget_control, info.maxx, sensitive = 0
                widget_control, info.miny, sensitive = 0
                widget_control, info.posy_lbl, sensitive = 0
		widget_control, info.maxy, sensitive = 0
		widget_control, info.multi_plots, sensitive = 0
		widget_control, info.multi_menu, sensitive = 0
		widget_control, info.txt_color, sensitive = 0
		for i=0,8 do $
		  widget_control, info.multi_plot[i], sensitive = 0
            end
            1: begin
                widget_control, info.minx, sensitive = 1
                widget_control, info.posx_lbl, sensitive = 1
                widget_control, info.maxx, sensitive = 1

                widget_control, info.miny, sensitive = 1
                widget_control, info.posy_lbl, sensitive = 1
                widget_control, info.maxy, sensitive = 1

		widget_control, info.text_label, sensitive = 1
		widget_control, info.text_size, sensitive = 1
		widget_control, info.text_thick, sensitive = 1

		widget_control, info.multi_plots, sensitive = 1
		widget_control, info.multi_plots, get_value = enabled_multi
		widget_control, info.txt_color, sensitive = 1
		for i=0, 8 do $
		  widget_control, info.multi_plot[i], sensitive = 0

            end
        endcase
    end

    'multi': begin
	widget_control, info.multi_plots, get_value = temp
	enabled = temp[0]

	case enabled of
	    0: begin
		widget_control, info.text_label, sensitive = 1
		widget_control, info.text_size, sensitive = 1
		widget_control, info.text_thick, sensitive = 1
		widget_control, info.multi_menu, sensitive = 0
		widget_control, info.multi_label, sensitive = 0
		for i=0,8 do $
	 	  widget_control, info.multi_plot[i], sensitive = 0
	    end
	    1: begin
		widget_control, info.text_label, sensitive = 0
		widget_control, info.text_size, sensitive = 1
		widget_control, info.text_thick, sensitive = 1
		widget_control, info.multi_menu, sensitive = 1
		widget_control, info.multi_label, sensitive = 1
		for i=0,8 do $
		  widget_control, info.multi_plot[i], sensitive = 0
	    end
	endcase
     end

     'multi_menu':begin
	xy = widget_info(info.multi_menu, /droplist_select)

	case xy of
	    ; 1x2
	    0: begin 
		widget_control, info.text_label, sensitive = 0
		widget_control, info.multi_plot[0], sensitive = 1
		widget_control, info.multi_plot[1], sensitive = 0
		widget_control, info.multi_plot[2], sensitive = 0
		widget_control, info.multi_plot[3], sensitive = 1
		widget_control, info.multi_plot[4], sensitive = 0
		widget_control, info.multi_plot[5], sensitive = 0
		widget_control, info.multi_plot[6], sensitive = 0
		widget_control, info.multi_plot[7], sensitive = 0
		widget_control, info.multi_plot[8], sensitive = 0
	    end
	    ; 2x1
	    1: begin
		widget_control, info.text_label, sensitive = 0
		widget_control, info.multi_plot[0], sensitive = 1
		widget_control, info.multi_plot[1], sensitive = 1
		widget_control, info.multi_plot[2], sensitive = 0
		widget_control, info.multi_plot[3], sensitive = 0
		widget_control, info.multi_plot[4], sensitive = 0
		widget_control, info.multi_plot[5], sensitive = 0
		widget_control, info.multi_plot[6], sensitive = 0
		widget_control, info.multi_plot[7], sensitive = 0
		widget_control, info.multi_plot[8], sensitive = 0
	    end
	    ; 2x2
	    2: begin
		widget_control, info.text_label, sensitive = 0
		widget_control, info.multi_plot[0], sensitive = 1
		widget_control, info.multi_plot[1], sensitive = 1
		widget_control, info.multi_plot[2], sensitive = 0
		widget_control, info.multi_plot[3], sensitive = 1
		widget_control, info.multi_plot[4], sensitive = 1
		widget_control, info.multi_plot[5], sensitive = 0
		widget_control, info.multi_plot[6], sensitive = 0
		widget_control, info.multi_plot[7], sensitive = 0
		widget_control, info.multi_plot[8], sensitive = 0
	    end
	    ; 2x3
	    3: begin
		widget_control, info.text_label, sensitive = 0
		widget_control, info.multi_plot[0], sensitive = 1
		widget_control, info.multi_plot[1], sensitive = 1
		widget_control, info.multi_plot[2], sensitive = 0
		widget_control, info.multi_plot[3], sensitive = 1
		widget_control, info.multi_plot[4], sensitive = 1
		widget_control, info.multi_plot[5], sensitive = 0
		widget_control, info.multi_plot[6], sensitive = 1
		widget_control, info.multi_plot[7], sensitive = 1
		widget_control, info.multi_plot[8], sensitive = 0
	    end
	    ; 3x2
	    4: begin
		widget_control, info.text_label, sensitive = 0
		widget_control, info.multi_plot[0], sensitive = 1
		widget_control, info.multi_plot[1], sensitive = 1
		widget_control, info.multi_plot[2], sensitive = 1
		widget_control, info.multi_plot[3], sensitive = 1
		widget_control, info.multi_plot[4], sensitive = 1
		widget_control, info.multi_plot[5], sensitive = 1
		widget_control, info.multi_plot[6], sensitive = 0
		widget_control, info.multi_plot[7], sensitive = 0
		widget_control, info.multi_plot[8], sensitive = 0
	    end
	    ; 3x3
	    5: begin
		widget_control, info.text_label, sensitive = 0
		widget_control, info.multi_plot[0], sensitive = 1
		widget_control, info.multi_plot[1], sensitive = 1
		widget_control, info.multi_plot[2], sensitive = 1
		widget_control, info.multi_plot[3], sensitive = 1
		widget_control, info.multi_plot[4], sensitive = 1
		widget_control, info.multi_plot[5], sensitive = 1
		widget_control, info.multi_plot[6], sensitive = 1
		widget_control, info.multi_plot[7], sensitive = 1
		widget_control, info.multi_plot[8], sensitive = 1
	    end
	endcase

    end
		


; *** done pressed ***
    'done': begin
        widget_control, info.use_label, get_value = enabled
        label_info.enabled = enabled[0]

	widget_control, info.multi_plots, get_value = enabled
	label_info.multi = enabled[0]

	widget_control, info.text_size, get_value = tmp
	label_info.size = tmp[0]

	widget_control, info.text_thick, get_value = tmp
	label_info.thick = tmp[0]

	label_info.multi_opt = widget_info(info.multi_menu, /droplist_select)

	if (label_info.multi) then begin
	  for i=0, 8 do begin
	     widget_control, info.multi_plot[i], get_value = tmp
	     label_info.label[i] = tmp
	  endfor
	   
	endif else begin
           widget_control, info.text_label, get_value = temp
           label_info.label[0] = temp[0]
	endelse

        widget_control, info.posx_lbl, get_value = temp
        label_info.posx = temp[0]

        widget_control, info.posy_lbl, get_value = temp
        label_info.posy = temp[0]

	label_info.color = widget_info(info.txt_color, /droplist_select)

	widget_control, ev.top, /destroy

    end

    else:
endcase

end




pro xlabel

common save_label, label_info

; get the list of colors
idummy = color('black', get_list = color_list)

; create the main widget base
main_base = widget_base(/column, title = 'floating label options')

; divide this base into vertical sections
information_base = widget_base(main_base, /column)
text_base         = widget_base(main_base, /row, /frame)
color_base	 = widget_base(main_base, /row, /frame)
slider_base      = widget_base(main_base, /row, /frame)
multi_base         = widget_base(main_base, /column, /frame)
end_base         = widget_base(main_base, /column, /frame)

use_label = cw_bgroup(information_base, 'use floating label', $
                     /nonexclusive, set_value = label_info.enabled, $
                     uvalue = 'flt_lbl')

text_label = cw_field(text_base, title = ' label text: ', $
                   xsize = 12, uvalue = 'text_label', $
                   value = label_info.label[0])

text_size  = cw_field(color_base, title = 'size: ', $
		xsize = 12, uvalue = 'text_size', $
		value = label_info.size)

text_thick = cw_field(color_base, title='line thickness: ', $
			xsize = 12, uvalue = 'text_thick', $
			value = label_info.thick)

txt_color = widget_droplist(color_base, title = 'color: ', uvalue = 'txt_color', $
			value = color_list)
widget_control, txt_color, set_droplist_select = label_info.color

minx = widget_label(slider_base, value = '0. ', $
                      uvalue = 'minx', /align_left)

posx_lbl = cw_fslider(slider_base, title = 'position in x ', $
                          maximum = 1.0, minimum = 0.0, value=label_info.posx, $
                          uvalue = 'posx')

maxx = widget_label(slider_base, value = '1.', $
                      uvalue = 'maxx', /align_left)

miny = widget_label(slider_base, value = '0. ', $
                      uvalue = 'miny', /align_left)

posy_lbl = cw_fslider(slider_base, title = 'position in y ', $
                          maximum = 1., minimum = 0., value=label_info.posy, $
                          uvalue = 'posy')

maxy = widget_label(slider_base, value = '1.', $
                      uvalue = 'maxy', /align_left)

multi_plots = cw_bgroup(multi_base, 'multiple plots', $
			/nonexclusive, set_value = label_info.multi, $
			uvalue = 'multi')
xycount = ['1x2','2x1','2x2','2x3','3x2','3x3'] 

multi_menu  = widget_droplist(multi_base, title = 'X/Y plot count', $
			uvalue = 'multi_menu', value = xycount)
widget_control, multi_menu, set_droplist_select = label_info.multi_opt

multi_label = widget_label(multi_base, /align_left, value = 'Labels (-1 for blanks):') 
widget_control, multi_label, sensitive = 0

multi_setup         = widget_base(multi_base, /row, /frame)
multi_setup2         = widget_base(multi_base, /row, /frame)
multi_setup3         = widget_base(multi_base, /row, /frame)

multi_plot = lonarr(9)

multi_plot[0] = cw_field(multi_setup, title='0:', $
			xsize =12, value = label_info.label[0],$
			uvalue=label_info.label[0])
multi_plot[1] = cw_field(multi_setup, title='1:', $
			xsize =12, value = label_info.label[1],$
			uvalue=label_info.label[1] )
multi_plot[2] = cw_field(multi_setup, title='2:', $
			xsize =12, value = label_info.label[2],$
			uvalue= label_info.label[2])
multi_plot[3] = cw_field(multi_setup2, title='3:', $
			xsize =12, value = label_info.label[3],$
			uvalue=label_info.label[3])
multi_plot[4] = cw_field(multi_setup2, title='4:', $
			xsize =12, value = label_info.label[4],$
			uvalue=label_info.label[4] )
multi_plot[5] = cw_field(multi_setup2, title='5:', $
			xsize =12, value = label_info.label[5],$
			uvalue=label_info.label[5] )
multi_plot[6] = cw_field(multi_setup3, title='6:', $
			xsize =12, value = label_info.label[6],$
			uvalue=label_info.label[6])
multi_plot[7] = cw_field(multi_setup3, title='7:', $
			xsize =12, value = label_info.label[7],$
			uvalue=label_info.label[7] )
multi_plot[8] = cw_field(multi_setup3, title='8:', $
			xsize =12, value = label_info.label[8],$
			uvalue=label_info.label[8] )


if NOT label_info.enabled then begin
    widget_control, use_label, sensitive = 1
    widget_control, text_label, sensitive = 0
    widget_control, text_size, sensitive = 0
    widget_control, text_thick, sensitive = 0
    widget_control, txt_color, sensitive = 0
    widget_control, minx, sensitive = 0
    widget_control, posx_lbl, sensitive = 0
    widget_control, maxx, sensitive = 0
    widget_control, miny, sensitive = 0
    widget_control, posy_lbl, sensitive = 0
    widget_control, maxy, sensitive = 0
    widget_control, multi_plots, sensitive = 0
    widget_control, multi_menu, sensitive = 0
    for i=0,8 do $
       widget_control, multi_plot[i], sensitive = 0
endif

if NOT label_info.multi then begin
    widget_control, multi_menu, sensitive = 0
    for i=0,8 do $
       widget_control, multi_plot[i], sensitive = 0
endif

if label_info.multi then begin
    ;xy = widget_info(info.multi_menu, /droplist_select)
    case label_info.multi_opt of
	0: begin
           widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 0
	   widget_control, multi_plot[2], sensitive = 0
	   widget_control, multi_plot[3], sensitive = 1
	   widget_control, multi_plot[4], sensitive = 0
	   widget_control, multi_plot[5], sensitive = 0
	   widget_control, multi_plot[6], sensitive = 0
	   widget_control, multi_plot[7], sensitive = 0
	   widget_control, multi_plot[8], sensitive = 0
	end
	1: begin
           widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 1
	   widget_control, multi_plot[2], sensitive = 0
	   widget_control, multi_plot[3], sensitive = 0
	   widget_control, multi_plot[4], sensitive = 0
	   widget_control, multi_plot[5], sensitive = 0
	   widget_control, multi_plot[6], sensitive = 0
	   widget_control, multi_plot[7], sensitive = 0
	   widget_control, multi_plot[8], sensitive = 0
	end
	2: begin
           widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 1
	   widget_control, multi_plot[2], sensitive = 0
	   widget_control, multi_plot[3], sensitive = 1
	   widget_control, multi_plot[4], sensitive = 1
	   widget_control, multi_plot[5], sensitive = 0
	   widget_control, multi_plot[6], sensitive = 0
	   widget_control, multi_plot[7], sensitive = 0
	   widget_control, multi_plot[8], sensitive = 0
	end
	3: begin
           widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 1
	   widget_control, multi_plot[2], sensitive = 0
	   widget_control, multi_plot[3], sensitive = 1
	   widget_control, multi_plot[4], sensitive = 1
	   widget_control, multi_plot[5], sensitive = 0
	   widget_control, multi_plot[6], sensitive = 1
	   widget_control, multi_plot[7], sensitive = 1
	   widget_control, multi_plot[8], sensitive = 0
	end
        4: begin
	   widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 1
	   widget_control, multi_plot[2], sensitive = 1
	   widget_control, multi_plot[3], sensitive = 0
	   widget_control, multi_plot[4], sensitive = 0
	   widget_control, multi_plot[5], sensitive = 0
	   widget_control, multi_plot[6], sensitive = 0
	   widget_control, multi_plot[7], sensitive = 0
	   widget_control, multi_plot[8], sensitive = 0
	end
	5: begin
           widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 1
	   widget_control, multi_plot[2], sensitive = 1
	   widget_control, multi_plot[3], sensitive = 1
	   widget_control, multi_plot[4], sensitive = 1
	   widget_control, multi_plot[5], sensitive = 1
	   widget_control, multi_plot[6], sensitive = 0
	   widget_control, multi_plot[7], sensitive = 0
	   widget_control, multi_plot[8], sensitive = 0
	end
	6: begin
           widget_control, text_label, sensitive = 0
	   widget_control, multi_plot[0], sensitive = 1
	   widget_control, multi_plot[1], sensitive = 1
	   widget_control, multi_plot[2], sensitive = 1
	   widget_control, multi_plot[3], sensitive = 1
	   widget_control, multi_plot[4], sensitive = 1
	   widget_control, multi_plot[5], sensitive = 1
	   widget_control, multi_plot[6], sensitive = 1
	   widget_control, multi_plot[7], sensitive = 1
	   widget_control, multi_plot[8], sensitive = 1
	end
    endcase
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
        slider_base: slider_base, $
	multi_base:multi_base, $
        text_base:text_base, $
        use_label:use_label,$
        text_label:text_label, $
	text_size:text_size, $
	text_thick:text_thick, $
	txt_color:txt_color, $
        minx:minx, $
        posx_lbl:posx_lbl, $
        maxx:maxx, $
        miny:miny, $
	posy_lbl:posy_lbl, $
	maxy:maxy,$
	multi_plots:multi_plots, $
	multi_plot:multi_plot, $
	multi_menu:multi_menu, $
	multi_label:multi_label, $
        end_base:end_base, $                 ; end base widget
        done:done}                           ; button for done


; register the info structure with the widget base
widget_control, main_base, set_uvalue=info, /no_copy

xmanager, 'xlabel', main_base

end











