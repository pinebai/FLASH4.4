
pro open_eps, name

   set_plot, 'ps'
   device, encapsulated=1, /portrait, font_size=13, /color
   !p.font = -1
   xyouts, 0, 0, '!5'
   device, xsize=11.0, ysize=8.0
   device, filename=name+'.eps'

end ;----------------------------------------------------------------------------

pro close_ps
   device, /close
   set_plot, 'x'
end ;----------------------------------------------------------------------------

pro plot_mom_sinks

  cmd = 'python clean_sinks_evol.py'
  print, cmd & spawn, cmd

  cmd = 'python clean_flashdat.py'
  print, cmd & spawn, cmd

  mom_sinks


  ; read mom_sinks.dat
  tab_sinks = dblarr(8,100000)
  openr, 1, 'mom_sinks.dat'
  line = dblarr(8)
  line_number = 0l
  readf, 1, format='(A)'
  while not eof(1) do begin
    readf, 1, line
    tab_sinks[*,line_number] = double(line)
    line_number = line_number + 1
  endwhile
  close, 1
  tab_sinks = tab_sinks[*,0:line_number-1]

  ; read flash.dat
  tab_grid = dblarr(12,100000)
  openr, 1, 'SMT.dat_cleaned'
  line = dblarr(12)
  line_number = 0l
  readf, 1, format='(A)'
  while not eof(1) do begin
    readf, 1, line
    tab_grid[*,line_number] = double(line)
    line_number = line_number + 1
  endwhile
  close, 1
  tab_grid = tab_grid[*,0:line_number-1]

  ; plot mom x
  setcolors
  outfile = 'momentum_x'
  open_eps, outfile
  x1 = 1e3
  x2 = 1e6
  y1 = -8.
  y2 =  8.
  plot, [x1,x2], [y1,y2], xstyle=1, ystyle=1, background=255, /nodata, /xlog, $
    xtitle=textoidl('time [yr]'), ytitle=textoidl('!8p_x!X [10^{36}_{ }g_{ }cm_{ }s^{-1}]')
  oplot, [x1,x2], [0.,0.], linestyle=1
  color = 9
  oplot, tab_grid [0,*]/3.154d7, tab_grid [6,*]/1.d36, linestyle=0, color=color
  draw_label, 'gas', 0.52, 0.75, color=color, /right
  color = 5
  oplot, tab_sinks[0,*]/3.154d7, tab_sinks[3,*]/1.d36, linestyle=0, color=color
  draw_label, 'sinks', 0.52, 0.3, color=color, /right
  close_ps
  print, outfile+'.eps  written.'

  ; plot mom y
  setcolors
  outfile = 'momentum_y'
  open_eps, outfile
  x1 = 1e3
  x2 = 1e6
  y1 = -3.
  y2 =  7.
  plot, [x1,x2], [y1,y2], xstyle=1, ystyle=1, background=255, /nodata, /xlog, $
    xtitle=textoidl('time [yr]'), ytitle=textoidl('!8p_y!X [10^{36}_{ }g_{ }cm_{ }s^{-1}]')
  oplot, [x1,x2], [1.989,1.989], linestyle=1
  color = 9
  oplot, tab_grid [0,*]/3.154d7, tab_grid [7,*]/1.d36, linestyle=0, color=color
  draw_label, 'gas', 0.5, 0.3, color=color, /right
  color = 5
  oplot, tab_sinks[0,*]/3.154d7, tab_sinks[4,*]/1.d36, linestyle=0, color=color
  draw_label, 'sinks', 0.5, 0.75, color=color, /right
  close_ps
  print, outfile+'.eps  written.'


stop

  

  return

end
