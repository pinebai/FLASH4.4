; NAME:
;       partvelvec
;
; PURPOSE:
;       This procedure plots particles and the velocity vectors
;       of particles (at the positions of the particles).
;
; CATEGORY:
;       Plotting, Two-dimensional.
;
; CALLING SEQUENCE:
;       PARTVELVEC, VELX, VELY, POSX, POSY [, X, Y]
;
; INPUTS:
;       VELX:  A 1D array, containing the x-components
;              of the particle velocities.
;       VELY:  A 1D array, containing the y-components
;              of the particle velocities.
;       POSX:  A 1D array, containing the x-components
;              of the particle positions.
;       POSY:  A 1D array, containing the y-components
;              of the particle positions.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       ABSCISSA:   Optional abscissae values. must be a vector.
;
;       ORDINATE:   Optional ordinate values. must be a vector.
;                   If only 'ABSCISSA' is specified, 'ORDINATE'
;                   is taken to be equal to 'ABSCISSA'
;
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
;       PARTICLE_WIDGET:  The structure containing information from the paricle
;                   menu.
;
;       PART_DATA:  The entire particle structure so that individual particle
;                   temperatures and eventually other values can be used.
;
;       PART_TEMP:  The particle temperatures color-scaled to the min and max
;                   from the xflash gui.  this is used over the dummy colormap
;                   when plotting particles colored by temperature.
;
;       COUNTER:    Tracks which file is being read in.  Important for trajectories.
;
;       OLD_PARTS:  An array that holds the position and sometimes temperature
;                   of the single particle at each timestep.
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
;         PARTVELVEC, VELX, VELY, POSX, POSY
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
;                    Kim Robinson (kim@flash.uchicago.edu)
;                    Aug. 2004, added coloring of particles by 
;                    temperature, tracing single particle trajectory,
;                    and coloring single particle trajectory by 
;                    temperature.
;-

PRO partvelvec, particles, PARTICLE_WIDGET=particle_widget, $
                ABSCISSA=abscissa, ORDINATE=ordinate, $
                OVER=over, FRACTION=fraction, LENGTH=length, COLOR=color, $
                MINMAG = minmag, MAXMAG = maxmag, $
                LEGEND = legend, LEGCLR = legclr, $
                PART_TEMP = part_temp, COUNTER = counter, $
                PARTICLE_TAG=particleTag, OLD_PARTS = old_parts, DIR = dir


;debug, '1.10 T.B. 1997-OCT-20' 

;--------------------------------------------
; Various settings, modify these to customize
;---------------------------------------------

forward_function sci_notat, number

cust={customize, $
      length: 0.08, $  ; Maximum vector length relative to plot region. (*)
      lengtharrow: 0.3, $  ; Length of arrowhead legs relative to vectorlength.
      angle: 22.5 }  ; 1/2 times the angle between the arrowhead legs.

; (*) Not used if keyword LENGTH is present

if (n_elements(minmag) EQ 0) then minmag = 0
if (n_elements(maxmag) EQ 0) then maxmag = 1.e30

if (n_elements(xskip) EQ 0) then xskip = 0
if (n_elements(yskip) EQ 0) then yskip = 0

if (n_elements(outline) EQ 0) then outline = 0

velx = particles.velx
vely = particles.vely
posx = particles.posx
posy = particles.posy
tags = particles.tag

; 'particle_widget' represents the state of the "particle options" widget
showtags     = particle_widget.show_tag
showvectors  = particle_widget.plot_vel
sym_size     = particle_widget.sym_size
traj         = particle_widget.traj
traj_color   = particle_widget.traj_color
partnum      = particle_widget.partnum
data_enabled = particle_widget.data_enabled
color_min    = particle_widget.color_min
color_max    = particle_widget.color_max
typvel       = particle_widget.typical_velocity

;---------------------
; Some error handling
;---------------------

on_error,2  ; Return to caller if an error occurs.

if (n_elements(abscissa) NE 0) then begin
  if (n_elements(ordinate) EQ 0) then ordinate = abscissa
  if (((size(abscissa))[0] NE 1) OR ((size(ordinate))[0] NE 1)) then $
  message,'ABSCISSA and ORDINATE must be vectors!'
endif

if n_elements(fraction) GT 0 then $
  if (fraction LT 0.0 OR fraction GT 1.0) then $
  message,'Fraction has to be between 0.0 and 1.0.'

;--------------
; Prepare plot
;--------------

numparticles = n_elements(velx)     ; Number of particles.
vel          = sqrt(velx^2+vely^2)  ; Total velocity.
maxvel       = max(vel)             ; Maximum velocity.

if n_elements(typvel) EQ 0 then typvel = maxvel

; Compute maximum length of vectors.
if n_elements(length) LE 0 then length=cust.length
minposx = min(posx)
maxposx = max(posx)
minposy = min(posy)
maxposy = max(posy)
; NB: in IDL, "a > b" means "max(a, b)"
length  = length * ((maxposx-minposx) > (maxposy-minposy))

; Convert velocities.
velx = length * velx / typvel
vely = length * vely / typvel
vel  = length * temporary(vel) / typvel

; modified -- MZ
; check to see if a minimum or maximum vector length was specified
if n_elements(minmag) EQ 0 then begin
  minmag = -1
endif else begin
  minmag = length * minmag / typvel
endelse

if n_elements(maxmag) EQ 0 then begin
  maxmag = 1.e30
endif else begin
  maxmag = length * maxmag / typvel
endelse

; Make sure no vectors extend beyond the plot data window.
x1=posx+velx  ; End of vector.
y1=posy+vely
if n_elements(over) EQ 0 then begin
  minposx = min(x1) < minposx
  maxposx = max(x1) > maxposx
  minposy = min(y1) < minposy
  maxposy = max(y1) > maxposy
endif

angle    = cust.angle * !dtor  ; Convert from degrees to radians.
sinangle = sin(angle)       ; Need these.
cosangle = cos(angle)


;-----------
; Plot axes
;-----------

if n_elements(color) EQ 0 then color=!p.color

if n_elements(over) EQ 0 then begin
  if n_elements(abscissa) EQ 0 then $
    plot,[minposx,maxposx],[minposy,maxposy], $
    /nodata,/xstyle,/ystyle,COLOR=color $
  else plot, abscissa, ordinate, /nodata, /xstyle, $
    /ystyle, COLOR=color
endif

;--------------
; Plot vectors
;--------------

if n_elements(fraction) GT 0 then begin
  if fraction EQ 1.0 then begin
    nrgood=long(fraction*numparticles)  ; # of vectors to plot.
    if nrgood EQ 0 then return
    ; Compute indices of vectors to plot. I use two lines to get more
    ; random "random numbers".
    good=long(randomu(seed,nrgood+1)*(numparticles-1.0))
    good=good(1:*)
    velx=temporary(velx(good))
    vely=temporary(vely(good))
    posx=posx(good)  ; Can't use temporary if we wan't to keep the data.
    posy=posy(good)
    x1=temporary(x1(good))
    y1=temporary(y1(good))
    numparticles=nrgood
  endif
endif

; define a user symbol that is a circle
NPTS = 24
tsym = findgen(NPTS)*2.*!pi/NPTS
xsym = cos(tsym)
ysym = sin(tsym)

usersym, xsym, ysym, /fill
;----------------------------------------------
; bypass loop if plotting single particle
;----------------------------------------------

if (traj) then begin
  ; find particle by tag
  particleIndex = (where(tags EQ particleTag))[0]

  if (where(tag_names(particles) EQ 'ptemp') GE 0) then begin
    ; This data-set has a field 'ptemp', which is supposed to
    ; be the temperature at the particle's position, so use it
    part_temp = particles.ptemp
  endif else begin
    ; This data-set doesn't have a 'ptemp', so we can't plot
    ; color according to temperature. Make an appropriately-
    ; sized array full of zeros to act as a stand-in for the
    ; missing 'ptemp'.
    ; This is just another hack made necessary by the fact that
    ; fidlr doesn't read information from the file at hand until
    ; plot-time. That is, the user shouldn't even have the option
    ; to select 'plot color by temperature' if the data-structure
    ; contains no 'ptemp', but there's no way to know whether it
    ; does or not when you're selecting options from the particle
    ; widget because 'read_amr', as I mentioned, won't be called
    ; until it's time to plot. Madness!
    part_temp = intarr((size(posx))[1])
  endelse

  oplot, [posx[particleIndex]], [posy[particleIndex]], psym=8, symsize=sym_size, color=color

  ;  tmp = ((counter+step)/step)-1-start
  if (traj_color) then begin
    if (counter EQ 0) then begin
      print, 'in here and [part_temp[particleIndex]] is: ', [part_temp[particleIndex]]
      old_parts[0,0]=[posx[particleIndex]]
      old_parts[1,0]=[posy[particleIndex]]
      old_parts[2,0]=[part_temp[particleIndex]]
    endif else begin
      if (counter GT 0) then begin
        old_parts[0,counter]=[posx[particleIndex]]
        old_parts[1,counter]=[posy[particleIndex]]
        old_parts[2,counter]=[part_temp[particleIndex]]

        end_traj = counter
        start_traj = counter-1

        iclrmap = color_index('Grayscale', MIN_VALUE=colorMin, MAX_VALUE=colorMax)
        loadct, iclrmap, FILE = dir + 'flash_colors.tbl', /SILENT

        color_traj = scale_color((old_parts[2,*]), $
          VARMAX = color_max, VARMIN = color_min, $
          COLORMAP_MIN=colorMin, COLORMAP_MAX=colorMax)

        while (end_traj NE 0) do begin
          oplot, [old_parts[0,end_traj],old_parts[0,start_traj]],$
            [old_parts[1,end_traj],old_parts[1,start_traj]],$
            COLOR=color_traj(start_traj), thick=5, noclip = 0

          end_traj = start_traj
          start_traj = start_traj -1 
        endwhile
      endif
    endelse

  endif else begin

    if (counter EQ 0) then begin
      old_parts[0,0]=[posx[particleIndex]]
      old_parts[1,0]=[posy[particleIndex]]
    endif else begin
      if (counter GT 0) then begin

        old_parts[0,counter]=[posx[particleIndex]]
        old_parts[1,counter]=[posy[particleIndex]]

        end_traj = counter
        start_traj = counter-1

        while (end_traj NE 0) do begin
          oplot, [old_parts[0,end_traj],old_parts[0,start_traj]],$
            [old_parts[1,end_traj],old_parts[1,start_traj]],$
            COLOR=color('gray'), thick=5, noclip = 0

          end_traj = start_traj
          start_traj = start_traj -1 
        endwhile
      endif
    endelse
  endelse
endif else begin

  for i=0l, numparticles-1l do begin

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
      if (showvectors EQ 1) then begin
        if (outline EQ 1) then begin
          plots, [posx(i),x1(i),x1(i)-(velx(i)*rcos+vely(i)*rsin)/vel(i), $
                  x1(i),x1(i)-(velx(i)*rcos-vely(i)*rsin)/vel(i)], $
                 [posy(i),y1(i),y1(i)-(vely(i)*rcos-velx(i)*rsin)/vel(i), $
                  y1(i),y1(i)-(vely(i)*rcos+velx(i)*rsin)/vel(i)], $
                 COLOR=color('white'), thick=3, noclip = 0
        endif

;        print, 'entering plots#1 where color is: ', color
        plots, [posx(i),x1(i),x1(i)-(velx(i)*rcos+vely(i)*rsin)/vel(i), $
                x1(i),x1(i)-(velx(i)*rcos-vely(i)*rsin)/vel(i)], $
               [posy(i),y1(i),y1(i)-(vely(i)*rcos-velx(i)*rsin)/vel(i), $
                y1(i),y1(i)-(vely(i)*rcos+velx(i)*rsin)/vel(i)], $
               COLOR=color, noclip = 0

      endif


      ;----------------------
      ; plot particles
      ;----------------------

      if (data_enabled EQ 1) then begin
        plots, posx(i), posy(i), noclip = 0, COLOR=part_temp(i), psym=8, symsize = sym_size
      endif else begin
        oplot, [posx(i)], [posy(i)], psym=8, symsize = sym_size, color=color, noclip=0
      endelse

      if (showtags EQ 1) then begin
        xyouts, posx(i), posy(i), ' '+ strtrim(string(i),2), noclip = 0
      endif
    endif
  endfor
endelse

if n_elements(legend) EQ 2 then begin

  if n_elements(legclr) EQ 0 then legclr = 0

  ptemp = convert_coord(legend[0], legend[1], /normal, /to_data)
  posx = ptemp[0]
  posy = ptemp[1]

  velx = length
  vely = 0
    
  vel = velx

  r = cust.lengtharrow*vel
  rcos = r*cosangle
  rsin = r*sinangle

  x1 = posx + velx
  y1 = posy + vely

  plots,[posx,x1,x1-(velx*rcos+vely*rsin)/vel,x1,x1-(velx*rcos-vely*rsin)/vel], $
    [posy,y1,y1-(vely*rcos-velx*rsin)/vel, y1,y1-(vely*rcos+velx*rsin)/vel], $
    COLOR=legclr


  velx = 0
  vely = length
    
  vel = vely

  x1 = posx + velx
  y1 = posy + vely

  plots,[posx,x1,x1-(velx*rcos+vely*rsin)/vel,x1,x1-(velx*rcos-vely*rsin)/vel], $
    [posy,y1,y1-(vely*rcos-velx*rsin)/vel, y1,y1-(vely*rcos+velx*rsin)/vel], $
    COLOR = legclr


  xyouts, x1 + length, posy + vely/2, sci_notat(typvel) + ' cm/s'
endif

end                             ; End of procedure PARTVELVEC.
