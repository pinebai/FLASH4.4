; CF 2009-2011
; draws a label 'string' in an active plot
; posx in NORMAL coords [0,1] (unless /DATA is set)
; posy in NORMAL coords [0,1] (unless /DATA is set)
; dx, dy, and textdx are in units of normal coordinates
pro draw_label, string, posx, posy, LINESTYLE=linestyle, COLOR=color, THICK=thick, PSYM=psym, SYMSIZE=symsize, orientation=orientation, $
                RIGHT=right, CENTRE=centre, CHARSIZE=charsize, CHARTHICK=charthick, DX=dx, DY=dy, TEXTDX=textdx, DATA=data

  if not keyword_set(orientation) then orientation = 0. ; in degrees
  if not keyword_set(charsize) then charsize = 1.0
  if not keyword_set(charthick) then charthick = 1.0
  xch = double(!d.x_ch_size)/!d.x_size*charsize ; x-charsize in norm coords
  ych = double(!d.y_ch_size)/!d.y_size*charsize ; y-charsize in norm coords
  if not keyword_set(color) then color = 0
  if not keyword_set(thick) then thick = 1.0
  if not keyword_set(symsize) then symsize = 1.0
  if not keyword_set(dx) then dx = 4.*xch
  if not keyword_set(dy) then dy = 0.7*ych/2.
  if not keyword_set(textdx) then textdx = 1.1
  if keyword_set(right) then begin ; right-justified
     alignment=1.0
     shift=textdx*dx
     l=-dx
  endif else begin ; left-justified
     alignment=0.0
     shift=-textdx*dx
     l=+dx
  endelse
  if keyword_set(centre) then begin
     alignment=0.5
  endif
  pt = [posx, posy]
  if keyword_set(data) then pt = CONVERT_COORD(pt[0], pt[1], /DATA, /TO_NORMAL)
  p1 = CONVERT_COORD(pt[0]+shift, pt[1]+dy, /NORMAL, /TO_DATA)
  p2 = CONVERT_COORD(pt[0]+l+shift, pt[1]+dy, /NORMAL, /TO_DATA)
  pm = CONVERT_COORD(pt[0]+l/2.+shift, pt[1]+dy, /NORMAL, /TO_DATA)     
  if keyword_set(linestyle) then oplot, [p1[0],p2[0]], [p1[1],p2[1]], linestyle=linestyle, color=color, thick=thick
  if keyword_set(psym) then oplot, [pm[0], pm[0]], [pm[1], pm[1]], color=color, psym=psym, symsize=symsize, thick=thick
  xyouts, pt[0], pt[1], string, /NORMAL, alignment=alignment, color=color, charsize=charsize, charthick=charthick, orientation=orientation

return
end ;----------------------------------------------------------------------------
