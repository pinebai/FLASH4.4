;  This file contains the particle widget event handler AND the
;   actual setup for the particle widget. 

;----------------------------------------------------
; particle widget event handler
pro xparticle_event, event

common save_parts, particle_widget

; get state information
;widget_control, event.top, get_uvalue=info
; lumley says this, with pointers
widget_control, event.top, get_uvalue=ptrInfo
info = *ptrInfo


; identify widget which caused the event
widget_control, event.id, get_uvalue=widgetWhich

; take action pased on which widget was changed
case widgetWhich of

  ; "plot particles" has been selected
  'part_plt': begin
    widget_control, info.part_plot, get_value = temp
    enabled = temp[0]

    ; enable/disable the various groups
    case enabled of
      0: begin   ; deselect
        widget_control, info.typ_vel, sensitive = 0
        widget_control, info.smaller, sensitive = 0
        widget_control, info.psym_size, sensitive = 0
        widget_control, info.larger, sensitive = 0
        widget_control, info.vel_plot, sensitive = 0
        widget_control, info.tag, sensitive = 0
        widget_control, info.data_plot, sensitive = 0
        widget_control, info.dummy_plot, sensitive = 0
        widget_control, info.plot_traj, sensitive = 0
        widget_control, info.traj_color, sensitive = 0
        widget_control, info.part_num, sensitive = 0
        widget_control, info.traj_min, sensitive = 0
        widget_control, info.traj_max, sensitive = 0
        widget_control, info.traj_log, sensitive = 0
      end
      1: begin  ; select
        widget_control, info.smaller, sensitive = 1
        widget_control, info.psym_size, sensitive = 1
        widget_control, info.larger, sensitive = 1

        widget_control, info.vel_plot, sensitive = 1

        widget_control, info.typ_vel, sensitive = 0

        widget_control, info.tag, sensitive = 1
        widget_control, info.data_plot, sensitive = 1
        ; widget_control, info.data_plot, get_value = enabled_data

        widget_control, info.dummy_plot, sensitive = 0
        widget_control, info.plot_traj, sensitive = 1
        widget_control, info.part_num, sensitive = 0
        widget_control, info.traj_color, sensitive = 0
        widget_control, info.traj_min, sensitive = 0
        widget_control, info.traj_max, sensitive = 0
        widget_control, info.traj_log, sensitive = 0
      end
    endcase   ; selected/deselected button for "plot particles" 
  end   ; of event handling for button "plot particles"

  ; event handling for selection of "plot velocity vectors"
  'vel_plt': begin
    widget_control, info.vel_plot, get_value = temp
    enabled = temp[0]
        
    ; case selection for enable/disable "plot velocity vectors"
    case enabled of
      0: begin
        widget_control, info.typ_vel, sensitive = 0
        widget_control, info.plot_traj, sensitive = 1
      end
      1: begin
        widget_control, info.typ_vel, sensitive = 1
        widget_control, info.plot_traj, sensitive = 0
      end
    endcase 
  end

  'data_plt': begin
    widget_control, info.data_plot, get_value = temp
    enabled = temp[0]

    case enabled of
      0: begin
        widget_control, info.dummy_plot, sensitive = 0
        widget_control, info.plot_traj, sensitive = 1
      end
      1: begin
        widget_control, info.dummy_plot, sensitive = 1
        widget_control, info.plot_traj, sensitive = 0
      end
    endcase
  end

  'traj':begin
    widget_control, info.plot_traj, get_value = temp
    enabled = temp[0]

    case enabled of
      0: begin
        widget_control, info.part_num, sensitive = 0
        widget_control, info.traj_color, sensitive = 0
        widget_control, info.vel_plot, sensitive = 1
        widget_control, info.data_plot, sensitive = 1
      end
      1: begin
        widget_control, info.part_num, sensitive = 1
        widget_control, info.traj_color, sensitive = 1
        widget_control, info.vel_plot, sensitive = 0
        widget_control, info.data_plot, sensitive = 0
      end
    endcase
  end

  'traj_color':begin
    widget_control, info.traj_color, get_value = temp
    enabled = temp[0]

    case enabled of
      0: begin
        widget_control, info.traj_min, sensitive = 0
        widget_control, info.traj_max, sensitive = 0
        widget_control, info.traj_log, sensitive = 0
      end
      1: begin
        widget_control, info.traj_min, sensitive = 1
        widget_control, info.traj_max, sensitive = 1
        widget_control, info.traj_log, sensitive = 1
      end
    endcase
   end

  ; *** done pressed ***
  'done': begin
    widget_control, info.part_plot, get_value = enabled
    particle_widget.enabled = enabled[0]

    widget_control, info.vel_plot, get_value = temp
    particle_widget.plot_vel = temp[0]

    widget_control, info.typ_vel, get_value = temp
    particle_widget.typical_velocity = temp[0]

    widget_control, info.psym_size, get_value = temp
    particle_widget.sym_size = temp[0]/4.

    widget_control, info.tag, get_value = temp
    particle_widget.show_tag = temp[0]

    widget_control, info.data_plot, get_value = enabled
    particle_widget.data_enabled = enabled[0]

    widget_control, info.dummy_plot, get_value= enabled
    particle_widget.dummy = enabled[0]

    widget_control, info.plot_traj, get_value = temp
    particle_widget.traj = temp[0]

    widget_control, info.traj_color, get_value = temp
    particle_widget.traj_color = temp[0]

    widget_control, info.part_num, get_value = temp
    particle_widget.partnum = temp[0]

    widget_control, info.traj_min, get_value = temp
    particle_widget.color_min = temp[0]

    widget_control, info.traj_max, get_value = temp
    particle_widget.color_max = temp[0]

    widget_control, info.traj_log, get_value = temp
    particle_widget.log = temp[0]

    ; destroy widget hierarchy
    widget_control, event.top, /destroy
  end

  else:
endcase  ; whichWidget was selected

; save state information
*ptrInfo = info

end  ; xparticle_event

;---------------------------------------------------------------------------------
;  the layout of the particle widge

pro xparticle 

common save_parts, particle_widget

; create the main widget base
main_base = widget_base(/column, title = 'particle options')

; divide this base into vertical sections
information_base = widget_base(main_base, /column)
vel_base         = widget_base(main_base, /column, /frame)
slider_base      = widget_base(main_base, /row, /frame)
tag_base         = widget_base(main_base, /column, /frame)
data_base	 = widget_base(main_base, /column, /frame)
traj_base        = widget_base(main_base, /column, /frame)
end_base         = widget_base(main_base, /column, /frame)

part_plot = cw_bgroup(information_base, 'plot particles', $
                     /nonexclusive, set_value = particle_widget.enabled, $
                     uvalue = 'part_plt')

vel_plot = cw_bgroup(vel_base, 'plot velocity vectors', $
;                    /nonexclusive, set_value = [1], $
                     /nonexclusive, set_value = particle_widget.plot_vel, $
                     uvalue = 'vel_plt')

typ_vel = cw_field(vel_base, title = ' typical velocity: ', $
                   xsize = 12, uvalue = 'typ_vel', $
                   value = particle_widget.typical_velocity)

smaller = widget_label(slider_base, value = 'smaller ', $
                       uvalue = 'smaller', /align_left)

temp = particle_widget.sym_size*4
psym_size = widget_slider(slider_base, title = 'symbol size ', $
                          maximum = 8, minimum = 1, value=temp, $
                          uvalue = 'psym')

larger = widget_label(slider_base, value = 'larger', $
                      uvalue = 'larger', /align_left)

data_plot = cw_bgroup(data_base, ' plot data', $
                      /nonexclusive, set_value = particle_widget.data_enabled,$
                      uvalue='data_plt')

dummy_plot = cw_bgroup(data_base, 'use dummy', $
                       /nonexclusive, set_value = particle_widget.dummy,$
                       uvalue='dummy')

tag = cw_bgroup(tag_base, 'show particle tags', $
                /nonexclusive, set_value = particle_widget.show_tag, $
                uvalue = 'tag_plt')

plot_traj = cw_bgroup(traj_base, 'Plot single particle trajectory', $
                      /nonexclusive, set_value = particle_widget.traj, $
                      uvalue = 'traj')

part_num = cw_field(traj_base, title ='particle number: ', $
                    xsize = 12, uvalue = 'part_num', /INT, $
                    value = particle_widget.partnum)

traj_color = cw_bgroup(traj_base, 'Grayscale by temperature', $
                       /nonexclusive, set_value = particle_widget.traj_color, $
                       uvalue = 'traj_color')

color_base = widget_base(traj_base, /row)

traj_min = cw_field(color_base, title='Grayscale Min: ',$
                    value = particle_widget.color_min, uvalue='clr_min')
traj_max = cw_field(color_base, title='Max: ', value = particle_widget.color_max, $
                    uvalue = 'clr_max')
traj_log = cw_bgroup(color_base, 'Log',$
                     /nonexclusive, set_value = particle_widget.log, uvalue='traj_log')


if not particle_widget.enabled then begin
  widget_control, typ_vel, sensitive = 0
  widget_control, smaller, sensitive = 0
  widget_control, psym_size, sensitive = 0
  widget_control, larger, sensitive = 0
  widget_control, vel_plot, sensitive = 0
  widget_control, tag, sensitive = 0
  widget_control, data_plot, sensitive = 0
  widget_control, dummy_plot, sensitive = 0
  widget_control, plot_traj, sensitive = 0
  widget_control, traj_color, sensitive = 0
  widget_control, part_num, sensitive = 0
  widget_control, traj_min, sensitive = 0
  widget_control, traj_max, sensitive = 0
  widget_control, traj_log, sensitive = 0
endif

if not particle_widget.plot_vel then begin
  widget_control, typ_vel, sensitive = 0
endif 

if not particle_widget.data_enabled then begin
  widget_control, dummy_plot, sensitive = 0
endif

if not particle_widget.traj then begin
  widget_control, traj_color, sensitive = 0
  widget_control, part_num, sensitive = 0
  widget_control, traj_min, sensitive = 0
  widget_control, traj_max, sensitive = 0
  widget_control, traj_log, sensitive = 0
endif

if not particle_widget.traj_color then begin
  widget_control, traj_min, sensitive = 0
  widget_control, traj_max, sensitive = 0
  widget_control, traj_log, sensitive = 0
endif

if particle_widget.plot_vel then begin
  widget_control, plot_traj, sensitive = 0
endif

if particle_widget.traj then begin
  widget_control, vel_plot, sensitive = 0
  widget_control, data_plot, sensitive = 0
endif

if particle_widget.data_enabled then begin
   widget_control, plot_traj, sensitive = 0
endif


;---------
; buttons
;---------
done = widget_button(end_base, value = 'Done', uvalue = 'done')


;-----------------
; draw the widget
;-----------------
widget_control, main_base, /realize


;--------------------------------------------------
; setup a structure to hold the widget information
;--------------------------------------------------
info = {main_base:main_base, $               ; main base
        information_base:information_base, $ ; information base
        slider_base: slider_base, $
        part_plot:part_plot, $
        vel_plot:vel_plot, $
        typ_vel:typ_vel, $
        smaller:smaller, $
        psym_size:psym_size, $
        larger:larger, $
        tag:tag, $
        data_plot:data_plot, $
        dummy_plot:dummy_plot,$
        plot_traj:plot_traj, $
        traj_color:traj_color, $
        traj_min:traj_min, $
        traj_max:traj_max, $
        traj_log:traj_log, $
        part_num:part_num, $
        end_base:end_base, $                 ; end base widget
        done:done}                           ; button for done


; register the info structure with the widget base
; use pointers as per lumley p 435
ptrInfo = ptr_new(info)
widget_control, main_base, set_uvalue=ptrInfo     ; was , /no_copy

;  set up the event handler as xparticle_event.pro
xmanager, 'xparticle', main_base

; cleanup 
result = *ptrInfo
ptr_free, ptrInfo

end   ; pro xparticle
