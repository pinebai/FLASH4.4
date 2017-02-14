; 24 bit color is necessary for concurrent multiple colormaps
device, retain=2, pseudo=8, decomposed = 0, true_color=24
window, /free, /pixmap, colors=-100         ; create window to allocate colors
plot, [0]                                  
wdelete,!d.window                          ; delete the window
device, set_character_size = [6,9]         ; set the vector font size

print, 'Number of colors allocated is ', !d.n_colors
