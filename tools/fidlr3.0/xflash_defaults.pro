;-----------------------------------------------------------------------------
; define some helper functions first
;-----------------------------------------------------------------------------

pro set_contour, contours, num_contours, variable, value, ctr_color

; find the first available contour
ictr = 0
iset = 0

for i = 0, num_contours-1 do begin
  if contours[i].enabled EQ 0 AND iset EQ 0 then begin
    iset = 1
    ictr = i
  endif
endfor


if var_index(variable) GE 0 then begin
  contours[ictr].enabled = 1
  contours[ictr].var     = var_index(variable)
  contours[ictr].value   = value
  contours[ictr].color   = color(ctr_color)
endif

end


pro set_limits, ctr_lim, variable, min, max

if var_index(variable) GE 0 then ctr_lim[var_index(variable),*] = [min, max]

end

;-----------------------------------------------------------------------------
; store the default contour ranges, etc. for the various problems
;
; 3-9-99  MZ
;-----------------------------------------------------------------------------

pro xflash_defaults, pblm_name, NUM_DEFAULTS=num_defaults, $
                     DEFAULT_NAMES=problem, $
                     CONTOURS=contours, NUM_CONTOURS=num_contours, $
                     VECTOR=vector, $
                     PARTICLE=particle, $
                     NBINS=nbins, $
                     HIST_SCALE=max_scale, $
                     CTR_LIM=ctr_lim, $
                     ISWP=iswp
                         



; allocate some storage
max_var = 100
ctr_lim = fltarr(max_var,2)

; set up the data structure for the contour options
contours = {enabled:0, $    
           var:0, $
           value:0.0, $
           color:0}

num_contours = 4
contours = replicate(contours, num_contours)

; setup the data structure for the particle options
particle = {enabled:0, $
            plot_vel:1, $
            sym_size: 1.0, $
            show_tag: 1, $
            typical_velocity:10.}

; setup the data structure for the vector options
vector = {enabled:0, $          ; set = 1 to plot vectors
          xcomp:'none', $
          ycomp:'none', $
          xskip:1, $            ; plot every xskip vectors in x dir
          yskip:1, $            ; plot every yskip vectors in y dir
          outline:0, $          ; draw a white outline around the vectors
          typical_vector:10., $ ; typical vector -- to scale 
          min_vector:1., $      ; minimum vector to plot
          max_vector:100.}      ; maximum vector to plot


; setup the default histogram info
nbins = 25
max_scale = 1.0

; set some defaults -- these are good for the abundances
ctr_lim[*,0] = .1
ctr_lim[*,1] = 1.


; list the problem names
problem = ['Generic',     $   
           'X-ray burst', $  
           'X-ray (cyl)', $  
           'Flame',       $  
           'Sedov',       $  
           'RT',          $  
           'Mach Step',   $  
           'Bondi Accretion', $  
           'Detonation']     

num_defaults = (size(problem))[1]

ipblm = (where(pblm_name EQ problem))[0]

if (ipblm EQ -1) then begin
  print, 'ERROR: default problem does not exist'
  return
endif


;-----------------------------------------------------------------------------
; Defining a new problem:
;  
;     To define a new problem, add the problem name to the problem
;     array above.  The problem number (ipblm) is then the array
;     index (starting from zero).  
;
;     The case statement below defines the problem specific information.
;     The following information must be set:
;
;         pblm_name   -- a short version of the problem name, used
;                           in xplot_amr to do problem specific
;                           plotting, if desired.
;
;         contour limits -- the contour limits are stored in the array,
;                           ctr_lim[ivar,0:1], where ivar is the variable
;                           number,
;
;                           The easiest way to get the variable
;                           number, since it can change with the data
;                           file, is to use the index function with the
;                           string abbreviation of the variable.  For
;                           example, to se the temperature contour
;                           limits, you would use:
;
;                             ctr_lim[index('temp',*)] = [1.e8, 1.e10]
;
;                           The second argument controls the
;                           minimum (0) and maximum (1) contour level
;                           
;
;         vskip          -- the default number of zones to skip in the
;                           x and y direction when drawing vector 
;                           vectors
; 
;         typvel         -- the typical vector in the problem, used
;                           to scale the vectors
;;
;         iswp           -- iswp = 1 puts the x coordinate on the
;                           vertical axis    
;                           iswp = 2 flips resulting images 180 degrees
;
; 
;      The easiest way to set these is to copy the case block of
;      another problem and set the number before the colon, ex:
; 
;           1: begin
;
;      to reflect the array index of the new problem.  To make the
;      changes take effect, this routine must be recompiled by idl, by
;      typing 
; 
;           .run xflash_defaults
;
;      at the idl prompt.
;  
;-----------------------------------------------------------------------------

case pblm_name of 

  ; *** X-ray burst ***
  'X-ray burst': begin

    ; set the contour limits
    set_limits, ctr_lim, 'temp',    1.e8,  6.e9
    set_limits, ctr_lim, 'dens',    1.e-5, 1.e8
    set_limits, ctr_lim, 'pres',    1.e16, 1.e26
    set_limits, ctr_lim, 'velx',    1.e5,  1.e10
    set_limits, ctr_lim, 'vely',    1.e5,  1.e10
    set_limits, ctr_lim, 'tot_vel', 1.e5,  1.e10
    set_limits, ctr_lim, 'ener',    1.e15, 1.e20
    set_limits, ctr_lim, 'snd_spd', 1.e5,  1.e10
    set_limits, ctr_lim, 'gamc',    4./3., 5./3.
    set_limits, ctr_lim, 'game',    4./3., 5./3.
    set_limits, ctr_lim, 'enuc',    1.e22, 1.e25

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 32
    vector.yskip          = 32
    vector.typical_vector = 5.e9
    vector.min_vector     = 1.e8
    vector.max_vector     = 3.e10

    ; plot 'x' in the file along the y axis
    iswp = 1

    ; define some reference lines
    set_contour, contours, num_contours, 'dens', 10, 'dkblue'

    set_contour, contours, num_contours, 'He4 ', .95, 'ltgreen'
    set_contour, contours, num_contours, 'he4 ', .95, 'ltgreen'

    set_contour, contours, num_contours, 'Ni56', .95, 'ltblue'
    set_contour, contours, num_contours, 'ni56', .95, 'ltblue'
  end

  ; *** X-ray burst ***
  'X-ray (cyl)': begin

    ; set the contour limits
    set_limits, ctr_lim, 'temp',    1.e8,  6.e9
    set_limits, ctr_lim, 'dens',    1.e-5, 1.e8
    set_limits, ctr_lim, 'pres',    1.e16, 1.e26
    set_limits, ctr_lim, 'velx',    1.e5,  1.e10
    set_limits, ctr_lim, 'vely',    1.e5,  1.e10
    set_limits, ctr_lim, 'tot_vel', 1.e5,  1.e10
    set_limits, ctr_lim, 'ener',    1.e15, 1.e20
    set_limits, ctr_lim, 'snd_spd', 1.e5,  1.e10
    set_limits, ctr_lim, 'gamc',    4./3., 5./3.
    set_limits, ctr_lim, 'game',    4./3., 5./3.
    set_limits, ctr_lim, 'enuc',    1.e22, 1.e25

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 32
    vector.yskip          = 32
    vector.typical_vector = 5.e9
    vector.min_vector     = 1.e8
    vector.max_vector     = 3.e10

    iswp = 0

    ; define some reference lines
    set_contour, contours, num_contours, 'dens', 10, 'dkblue'

    set_contour, contours, num_contours, 'He4 ', .95, 'ltgreen'
    set_contour, contours, num_contours, 'he4 ', .95, 'ltgreen'

    set_contour, contours, num_contours, 'Ni56', .95, 'ltblue'
    set_contour, contours, num_contours, 'ni56', .95, 'ltblue'
  end

  ; *** flame ***
  'Flame': begin

    ; set the contour limits
    set_limits, ctr_lim, 'temp',    1.e8,  1.e10
    set_limits, ctr_lim, 'dens',    2.e9,  8.e9
    set_limits, ctr_lim, 'pres',    1.e16, 1.e26
    set_limits, ctr_lim, 'ener',    1.e15, 1.e20
    set_limits, ctr_lim, 'velx',    1.e4,  1.e8
    set_limits, ctr_lim, 'vely',    1.e4,  1.e8
    set_limits, ctr_lim, 'tot_vel', 1.e4,  1.e8
    set_limits, ctr_lim, 'snd_spd', 1.e7,  1.e9
    set_limits, ctr_lim, 'gamc',    4./3., 5./3.
    set_limits, ctr_lim, 'game',    4./3., 5./3.
    set_limits, ctr_lim, 'enuc',    1.e22, 1.e25

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 32
    vector.yskip          = 32
    vector.typical_vector = 5.e6
    vector.min_vector     = 1.e4
    vector.max_vector     = 3.e10

    iswp = 0
  end

  ; *** blast wave ***
  'Sedov': begin

    set_limits, ctr_lim, 'temp',    1, 100
    set_limits, ctr_lim, 'dens',    1, 100
    set_limits, ctr_lim, 'pres',    1, 100
    set_limits, ctr_lim, 'ener',    1, 100
    set_limits, ctr_lim, 'velx',    1, 100
    set_limits, ctr_lim, 'vely',    1, 100
    set_limits, ctr_lim, 'tot_vel', 1, 100
    set_limits, ctr_lim, 'snd_spd', 1, 100
    set_limits, ctr_lim, 'gamc',    1, 100
    set_limits, ctr_lim, 'game',    1, 100
    set_limits, ctr_lim, 'enuc',    1, 100

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 4
    vector.yskip          = 4
    vector.typical_vector = 10.
    vector.min_vector     = 1.
    vector.max_vector     = 1000.

    iswp = 0
  end

  ; *** Rayleigh-Taylor ***
  'RT': begin

    set_limits, ctr_lim, 'temp',    1, 100
    set_limits, ctr_lim, 'dens',    1, 100
    set_limits, ctr_lim, 'pres',    1, 100
    set_limits, ctr_lim, 'ener',    1, 100
    set_limits, ctr_lim, 'velx',    1, 100
    set_limits, ctr_lim, 'vely',    1, 100
    set_limits, ctr_lim, 'tot_vel', 1, 100
    set_limits, ctr_lim, 'snd_spd', 1, 100
    set_limits, ctr_lim, 'gamc',    1, 100
    set_limits, ctr_lim, 'game',    1, 100
    set_limits, ctr_lim, 'enuc',    1, 100

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 4
    vector.yskip          = 4
    vector.typical_vector = 10.
    vector.min_vector     = 1.
    vector.max_vector     = 1000.

    iswp = 2
  end

  ; *** Mach Step ***
  'Mach Step': begin

    set_limits, ctr_lim, 'temp',    1, 100
    set_limits, ctr_lim, 'dens',    1, 100
    set_limits, ctr_lim, 'pres',    1, 100
    set_limits, ctr_lim, 'ener',    1, 100
    set_limits, ctr_lim, 'velx',    1, 100
    set_limits, ctr_lim, 'vely',    1, 100
    set_limits, ctr_lim, 'tot_vel', 1, 100
    set_limits, ctr_lim, 'snd_spd', 1, 100
    set_limits, ctr_lim, 'gamc',    1, 100
    set_limits, ctr_lim, 'game',    1, 100
    set_limits, ctr_lim, 'enuc',    1, 100

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 4
    vector.yskip          = 4
    vector.typical_vector = 10.
    vector.min_vector     = 1.
    vector.max_vector     = 1000.

    iswp = 0
  end

  ; *** detonation ***
  'Detonation': begin

    ; set the contour limits
    set_limits, ctr_lim, 'temp',    1.e8,   1e10
    set_limits, ctr_lim, 'dens',    1.0e7,  4.0e7
    set_limits, ctr_lim, 'pres',    8.7e23, 1.3e25
    set_limits, ctr_lim, 'ener',    1.e15,  1.e20
    set_limits, ctr_lim, 'velx',    1.e4,   1.e8
    set_limits, ctr_lim, 'vely',    1.e4,   1.e8
    set_limits, ctr_lim, 'tot_vel', 1.e4,   1.e8
    set_limits, ctr_lim, 'snd_spd', 1.e7,   1.e9
    set_limits, ctr_lim, 'gamc',    4./3.,  5./3.
    set_limits, ctr_lim, 'game',    4./3.,  5./3.
    set_limits, ctr_lim, 'enuc',    1.e22,  1.e25

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 16
    vector.yskip          = 16
    vector.typical_vector = 1.e8
    vector.min_vector     = 1.e4
    vector.max_vector     = 3.e10

    iswp = 1
  end

  ; *** Bondi Accretion ***
  'Bondi Accretion': begin

    ; set the contour limits
    set_limits, ctr_lim, 'temp',    .9e6,   6.e6
    set_limits, ctr_lim, 'tion',    .9e6,   6.e6
    set_limits, ctr_lim, 'tele',    .9e6,   6.e6
    set_limits, ctr_lim, 'trad',    4900.,  26000.
    set_limits, ctr_lim, 'pres',    .8,     550.
    set_limits, ctr_lim, 'pele',    1.6e-4, .5e-3
    set_limits, ctr_lim, 'pele',    1.6e-4, .5e-3
    set_limits, ctr_lim, 'prad',    .8,     550.
    set_limits, ctr_lim, 'dena',    .6e-18,  12.e-18
    set_limits, ctr_lim, 'dens',    .6e-18,  12.e-18
    set_limits, ctr_lim, 'ener',    2.5e18,  1.6e21
    set_limits, ctr_lim, 'eint',    2.5e18,  1.6e21
    set_limits, ctr_lim, 'eion',    0.,      1.8e17
    set_limits, ctr_lim, 'eele',    1.6e17,  1.8e17
    set_limits, ctr_lim, 'erad',    2.5e18,  1.6e21
    set_limits, ctr_lim, 'ugtm',    2.4,     1650.
    set_limits, ctr_lim, 'urda',    2.4,     1650.
    set_limits, ctr_lim, 'urdd',    -10.,   0.1
    set_limits, ctr_lim, 'urdq',    0.97,   1.01
    set_limits, ctr_lim, 'fllm',    .5e-6,  24.e-6
    set_limits, ctr_lim, 'flxl',    .5e-6,  24.e-6
    set_limits, ctr_lim, 'flla',    .5e-6,  24.e-6
    set_limits, ctr_lim, 'flad',    -6.e-5, 6.e-6
    set_limits, ctr_lim, 'flaq',    .1,     1.4
    set_limits, ctr_lim, 'eddi',    .5e-6,  24.e-6
    set_limits, ctr_lim, 'mar0',    -6.e15, 6.e15
    set_limits, ctr_lim, 'mar1',    -6.e15, 6.e15
    set_limits, ctr_lim, 'vela',    -3.e7,  100.
    set_limits, ctr_lim, 'velx',    -3.e7,  100.
    set_limits, ctr_lim, 'vely',    -3.e7,  100.
    set_limits, ctr_lim, 'tot_vel', 0.,     3.e7
    set_limits, ctr_lim, 'snd_spd', 1.e7,   1.e9
    set_limits, ctr_lim, 'gamc',    1.0,    5./3.
    set_limits, ctr_lim, 'game',    1.0,    5./3.

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 16
    vector.yskip          = 16
    vector.typical_vector = 1.e6
    vector.min_vector     = 1.e4
    vector.max_vector     = 3.e10

    iswp = 0
  end

  ; *** generic ***
  'Generic': begin

    ; set the contour limits
    set_limits, ctr_lim, 'temp',    1.e5,   1e10
    set_limits, ctr_lim, 'dens',    1.e2,   1.e8
    set_limits, ctr_lim, 'pres',    1.e22,  1.e28
    set_limits, ctr_lim, 'ener',    1.e10,  1.e20
    set_limits, ctr_lim, 'velx',    1.e3,   1.e9
    set_limits, ctr_lim, 'vely',    1.e3,   1.e9
    set_limits, ctr_lim, 'tot_vel', 1.e3,   1.e9
    set_limits, ctr_lim, 'snd_spd', 1.e3,   1.e9
    set_limits, ctr_lim, 'gamc',    4./3.,  5./3.
    set_limits, ctr_lim, 'game',    4./3.,  5./3.
    set_limits, ctr_lim, 'enuc',    1.e10,  1.e28

    ; set the default vector arrow density to plot, typical vector,
    ; and minimum and maximum velocities to plot
    vector.enabled        = 0
    vector.xcomp          = "velx"
    vector.ycomp          = "vely"
    vector.xskip          = 16
    vector.yskip          = 16
    vector.typical_vector = 1.e6
    vector.min_vector     = 1.e4
    vector.max_vector     = 3.e10

    iswp = 0
  end

endcase

return
end
