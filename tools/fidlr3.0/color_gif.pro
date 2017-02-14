; quick hack to get color gif off the screen 
; get the plot on screen, and then type 
; IDL> color_gif, 'filename.gif'

pro color_gif, filename

; save the various plot options -- so we don't mess anything up . . .
old_plot = !p

set_plot, 'Z'
tvlct, red, green, blue, /get

set_plot, 'X'
!p = old_plot

a = tvrd()

write_gif, filename, a, red, green, blue

end

