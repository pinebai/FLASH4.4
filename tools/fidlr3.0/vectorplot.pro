; NAME:
;       vectorplot
;
; PURPOSE:
;       This procedure plots velocity vectors throughout
;       the computational domain
;
; CATEGORY:
;       Plotting, Two-dimensional.
;
; CALLING SEQUENCE:
;       VECTORPLOT, VELX, VELY, POSX, POSY [, X, Y]
;
; INPUTS:
;       VELX:  A 2D array, containing the x-components
;              of the velocity at each cooridnate pair.
;       VELY:  A 2D array, containing the y-components
;              of the velocity at each cooridnate pair.
;
; OPTIONAL INPUTS:
;       X:   Optional abcissae values. X must be a vector.
;       Y:   Optional ordinate values. Y must be a vector. If only X
;            is specified, then Y is taken equal to be equal to X.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       FRACTION:   The fraction of the vectors to plot. They are
;                   taken at random from the complete sample.    
;                   Default is FRACTION = 1.0, use all vectors
;
;       LENGTH:     The maximum vectorlength relative to the plot data
;                   window.   
;                   Default = 0.08
;
;       COLOR:      The color for the vectors, axes and titles.
;                   Default=!P.COLOR
;
;       OVER:       Plot over the previous plot
;
;       MINMAG:     The minimum magnitude vector to plot
;                   Default is plot all
;
;       MAXMAG:     The maximum magnitude vector to plot
;                   Default is 1.e30
;
;       TYPVEL:     The magnitude to scale all velocities to.  This is
;                   useful if you want the arrow length to mean the
;                   same thing across plots.  
;                   The default is the maximum velocity
;
;       XSKIP &     Plot every xskip or yskip vector in x or y respectively.
;       YSKIP:      This allows you to thin out the velocity field in
;                   the x and y directions independently. 
;      
;       LEGEND:     The x & y NORMAL coordinates to plot a velocity legend
;
;       OUTLINE:    Draw a white outline around the arrows --
;                   increases contrast
;
;
;       COUNTER:    Tracks which file is being read in.  Important for trajectories.
;
;
;       DIR:        The xflash directory.  Important for loading a seperate colormap.
;
;       Plot        All other keywords available to PLOT are also used
;       Keywords:   by this procedure.
;
; OUTPUTS:
;       This procedure plots the velocity vectors (VELX,VELY) at the
;       positions of the particles, (POSX,POSY). If X and Y are not
;       specified, then the size of the plot is such that all vectors
;       just fit within in the plot data window.
;
; SIDE EFFECTS:
;       Plotting on the current device is performed.
;
; EXAMPLE:
;       Generate some particle positions and velocities.
;
;         POSX=RANDOMU(seed,200)
;         POSY=RANDOMU(seed,200)
;         VELX=RANDOMU(seed,200)-0.5
;         VELY=RANDOMU(seed,200)-0.5
;
;       Plot the particle velocities.
;
;         VECTORPLOT, VELX, VELY, POSX, POSY
;
; MODIFICATION HISTORY:
;       Written by:  Joop Schaye (jschaye@astro.rug.nl), Sep 1996.
;
;       Modified:    Theo Brauers (th.brauers@fz-juelich.de) Oct. 1997
;                    use with maps, incl. debug
;             
;                    Michael Zingale (zingale@oddjob.uchicago.edu)
;                    Aug. 1998, added minimum and maximum cutoff for 
;                    velocity, legend plotting, clipping of velocity
;                    vectors outside the plot window, option of
;                    skipping in x and y directions, and scaling to a
;                    typical velocity instead of the max if desired.
;

PRO vectorplot, velx, vely, posx, posy, x, y, $
                OVER=over, LENGTH=length, COLOR=color, _EXTRA=extra, $
                MINMAG = minmag, MAXMAG = maxmag, XSKIP = xskip, $
                YSKIP = yskip, TYPVEL = typvel, LEGEND = legend, $
                LEGCLR = legclr, TAG=tag, SYM_SIZE=sym_size, $
                OUTLINE=outline, COUNTER = counter, DIR = dir


;debug, '1.10 T.B. 1997-OCT-20' 

forward_function sci_notat

cust={customize, $
      length: 0.08, $  ; Maximum vector length relative to plot region. (*)
      lengtharrow: 0.3, $  ; Length of arrowhead legs relative to vectorlength.
      angle: 22.5 }  ; 1/2 times the angle between the arrowhead legs.

; (*) Not used if keyword LENGTH is present

if (n_elements(minmag) EQ 0) then minmag = 0
if (n_elements(maxmag) EQ 0) then maxmag = 1.e30

if (n_elements(xskip) EQ 0) then xskip = 0
if (n_elements(yskip) EQ 0) then yskip = 0

if (n_elements(sym_size) EQ 0) then sym_size=0
if (n_elements(outline) EQ 0) then outline = 0


;---------------------
; Some error handling
;---------------------

on_error,2  ; Return to caller if an error occurs.

nparams=n_params()
if nparams NE 4 then begin
  if (nparams NE 5 AND nparams NE 6) then begin
    message,'Wrong number of parameters!',/continue
    message,'Syntax: VECTORPLOT, VELX, VELY, POSX, POSY [, X, Y]', /noname,/noprefix
  endif

  if nparams EQ 5 then y=x
  sizex=size(x)
  sizey=size(y)
  if (sizex(0) NE 1 OR sizey(0) NE 1) then $
    message,'X and Y must be vectors!'
endif

sizevelx=size(velx)
sizevely=size(vely)
sizeposx=size(posx)
sizeposy=size(posy)

ndims = sizevelx[0]

if (total(sizevelx(0:ndims)-sizevely(0:ndims)) NE 0 $
  OR total(sizevelx(0:ndims)-sizeposx(0:ndims)) NE 0 $
  OR total(sizevelx(0:ndims)-sizeposy(0:ndims)) NE 0) then $
  message,'All arguments must have the same dimension and size!'

;---------------------------------------------------
; thin out the velocities if xskip or yskip are set
;---------------------------------------------------

xskip = fix(xskip)
yskip = fix(yskip)
if xskip LT 1 then xskip = 1
if yskip LT 1 then yskip = 1

if n_elements(xskip) GT 0 then begin

  if xskip GT n_elements(posx[*,0]) then begin
    msg1 = strjoin(["Unable to plot vectors: value of xskip (", strtrim(xskip,1), ") is larger"])
    msg2 = strjoin(["than the number of cells in the x-direction, (", strtrim(n_elements(posx[*,0]),1), ") ."])
    msg3 = ""
    msg4 = "Use 'Vector Options' button to adjust."

    result = dialog_message([msg1, msg2, msg3, msg4])
    return
  endif

  xindices = lindgen(n_elements(posx[*,0])/xskip)*xskip

  posx = posx[xindices,*]
  posy = posy[xindices,*]

  velx = velx[xindices,*]
  vely = vely[xindices,*]
endif

if n_elements(yskip) GT 0 then begin

  ; DEV not sure at all if I'm checking the right values
  ; here or delivering the right warning message -nttaylor
  if yskip GT n_elements(posx[*,0]) then begin
    msg1 = strjoin(["Unable to plot vectors: value of yskip (", strtrim(yskip,1), ") is larger"])
    msg2 = strjoin(["than the number of cells in the y-direction, (", strtrim(n_elements(posx[*,0]),1), ") ."])
    msg3 = ""
    msg4 = "Use 'Vector Options' button to adjust."

    result = dialog_message([msg1, msg2, msg3, msg4])
    return
  endif

  yindices = lindgen(n_elements(posx[0,*])/yskip)*yskip

  posx = posx[*,yindices]
  posy = posy[*,yindices]

  velx = velx[*,yindices]
  vely = vely[*,yindices]
endif

;--------------
; Prepare plot
;--------------

nvecs=n_elements(velx)  ; Number of cells
vel=sqrt(velx^2+vely^2)  ; Total velocity.

maxvel=max(vel)  ; Maximum velocity.
if n_elements(typvel) EQ 0 then typvel = maxvel

; Compute maximum length of vectors.
if n_elements(length) LE 0 then length=cust.length
minposx=min(posx)
maxposx=max(posx)
minposy=min(posy)
maxposy=max(posy)
length=length*((maxposx-minposx) > (maxposy-minposy))

; Convert velocities.
vx=length*velx/typvel
vy=length*vely/typvel
vel=length*temporary(vel)/typvel

; modified -- MZ
; check to see if a minimum or maximum vector length was specified
if n_elements(minmag) EQ 0 then begin
  minmag = -1
endif else begin
  minmag = length*minmag/typvel
endelse

if n_elements(maxmag) EQ 0 then begin
  maxmag = 1.e30
endif else begin
  maxmag = length*maxmag/typvel
endelse

; Make sure no vectors extend beyond the plot data window.
x1=posx+vx  ; End of vector.
y1=posy+vy
if (nparams EQ 4 and n_elements(over) EQ 0) then begin
  minposx=min(x1)<minposx
  maxposx=max(x1)>maxposx
  minposy=min(y1)<minposy
  maxposy=max(y1)>maxposy
endif


angle=cust.angle*!dtor  ; Convert from degrees to radians.
sinangle=sin(angle)  ; Need these.
cosangle=cos(angle)


;-----------
; Plot axes
;-----------

if n_elements(color) EQ 0 then color=!p.color

if n_elements(over) EQ 0 then begin
  if nparams EQ 4 then $
    plot,[minposx,maxposx],[minposy,maxposy], $
    /nodata,/xstyle,/ystyle,COLOR=color,_EXTRA=extra $
  else plot,x,y,/nodata,/xstyle,/ystyle,COLOR=color,_EXTRA=extra
endif

;--------------
; Plot vectors
;--------------

if n_elements(fraction) GT 0 then begin
  if fraction EQ 1.0 then goto, plotall
  nrgood=long(fraction*nvecs)  ; # of vectors to plot.
  if nrgood EQ 0 then return
  ; Compute indices of vectors to plot. I use two lines to get more
  ; random "random numbers".
  good=long(randomu(seed,nrgood+1)*(nvecs-1.0))
  good=good(1:*)
  vx=temporary(vx(good))
  vy=temporary(vy(good))
  px=posx(good)  ; Can't use temporary if we wan't to keep the data.
  py=posy(good)
  x1=temporary(x1(good))
  y1=temporary(y1(good))
  nvecs=nrgood
endif else begin
  plotall:
    px=posx
    py=posy
endelse

for i=0l, nvecs-1l do begin

  ; Note that we cannot put the next three lines outside the loop,
  ; because we want the arrow size to be relative to the vector length.
  r=cust.lengtharrow*vel(i)  ; Length of arrow head.
  rsin=r*sinangle
  rcos=r*cosangle

  ; Draw basis, arrow leg, same arrow leg, other arrow leg.
  ; One arrow leg is drawn twice, because we need to return to the end
  ; of the vector to draw the other leg.

  ; modified -- MZ

;    print, 'now minmag is: ', minmag
;    print, 'now maxmag is: ', maxmag
;    print, 'and vel(i) is: ', vel(i)

  if ((vel(i) GE minmag) and (vel(i) LE maxmag)) then begin
    if (outline EQ 1) then begin
      plots, [px(i),x1(i),x1(i)-(vx(i)*rcos+vy(i)*rsin)/vel(i), $
              x1(i),x1(i)-(vx(i)*rcos-vy(i)*rsin)/vel(i)], $
             [py(i),y1(i),y1(i)-(vy(i)*rcos-vx(i)*rsin)/vel(i), $
              y1(i),y1(i)-(vy(i)*rcos+vx(i)*rsin)/vel(i)], $
             COLOR=color('white'), thick=3, noclip = 0
    endif

;    print, 'entering plots#1 where color is: ', color
    plots, [px(i),x1(i),x1(i)-(vx(i)*rcos+vy(i)*rsin)/vel(i), $
            x1(i),x1(i)-(vx(i)*rcos-vy(i)*rsin)/vel(i)], $
           [py(i),y1(i),y1(i)-(vy(i)*rcos-vx(i)*rsin)/vel(i), $
            y1(i),y1(i)-(vy(i)*rcos+vx(i)*rsin)/vel(i)], $
           COLOR=color, noclip = 0
  endif
endfor


if n_elements(legend) EQ 2 then begin

  if n_elements(legclr) EQ 0 then legclr = 0

  ptemp = convert_coord(legend[0], legend[1], /normal, /to_data)
  px = ptemp[0]
  py = ptemp[1]

  vx = length
  vy = 0
    
  vel = vx

  r = cust.lengtharrow*vel
  rcos = r*cosangle
  rsin = r*sinangle

  x1 = px + vx
  y1 = py + vy

  plots,[px,x1,x1-(vx*rcos+vy*rsin)/vel,x1,x1-(vx*rcos-vy*rsin)/vel], $
    [py,y1,y1-(vy*rcos-vx*rsin)/vel, y1,y1-(vy*rcos+vx*rsin)/vel], $
    COLOR=legclr


  vx = 0
  vy = length
    
  vel = vy

  x1 = px + vx
  y1 = py + vy

  plots,[px,x1,x1-(vx*rcos+vy*rsin)/vel,x1,x1-(vx*rcos-vy*rsin)/vel], $
    [py,y1,y1-(vy*rcos-vx*rsin)/vel, y1,y1-(vy*rcos+vx*rsin)/vel], $
    COLOR = legclr


  xyouts, x1 + length, py + vy/2, sci_notat(typvel) + ' cm/s'
endif

end  ; End of procedure VECTORPLOT.
