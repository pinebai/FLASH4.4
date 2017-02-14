;===========================================
;  plot particles read with read_particles
;===========================================
PRO plot_particles,particles,part_names,X=xname,Y=yname, $
                   VERBOSE=verbose

; check arguments
IF (n_elements(verbose) EQ 0) then verbose = 0
IF (n_elements(xname) EQ 0) THEN xname='posx'
IF (n_elements(yname) EQ 0) THEN yname='posy'

IF (verbose) THEN print, 'plotting ',yname,' versus ',xname

xindex = where(part_names EQ xname)
yindex = where(part_names EQ yname)

x = particles(xindex,*)
y = particles(yindex,*)

plot,x,y, psym=1



END ; plot_particles.pro
