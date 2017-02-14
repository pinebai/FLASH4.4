;================
; read dump files created by Particle_dump
;=============================================

PRO particles_dump, filename, PARTICLES=particles, DT=dt
         
openr, file_id, filename, /get_lun

; get processor
fmt_p = '(a3)'
processor=''
readf, file_id, processor,format=fmt_p

; begin loop over data
fmt_info='(I2,I2,F,I2)'
info={block:0L, step:0L, time:0.0, dt:0.0, npart:0L}

; dimension by huge array for particles
maxparticles=1000
maxtimes = 100
numParticles = 0
numTimes = 0L
single = {index:0L, tag:0L, block:0L, posx:0.0, posy:0.0, posz:0.0, $
             velx:0.0, vely:0.0, velz:0.0}
particles = replicate(single,maxparticles,maxtimes)
dt = fltarr(maxtimes)



;  read until end of file
loop = -1
WHILE (eof(file_id) NE 1) DO  BEGIN

    readf,file_id,info
    loop = loop + 1

    ; store delta time
    dt(loop) = info.dt

    ; time step and error checking
    time_index = info.step - 1
    IF (time_index eq maxtimes) THEN BEGIN
        free_lun, file_id
        message, 'Maximum time steps reached; please redimension'
    ENDIF 
    numTimes = numTimes + 1

    nParticles = info.npart
    IF (nParticles GE maxparticles) THEN BEGIN
        free_lun, file_id
        message, 'Number of particles too big; please redimension'
    ENDIF 
    ; find the maximum number of particles
    IF (nParticles GT numParticles) THEN numParticles=nParticles

    ; now read the particles
    FOR i=0l, nParticles-1 DO BEGIN 
        readf,file_id,single
        particles[time_index,i] = single
        ; assign the time step as the first value in the array
        particles[time_index,i].index = info.step
    ENDFOR 
ENDWHILE 


; close the file
free_lun,file_id
   
; find all possible tags
tagIndex = particles[*,*].tag
; sort and find unique tags, per p. 80 Gumley
tagIndex = tagIndex[sort(tagIndex)]
tagIndex = tagIndex[uniq(tagIndex)]
tagIndex = tagIndex[where(tagIndex GT 0)] ; remove possible zero data
tagIndex = reform(tagIndex)

; cut down data to size 
numParticles=n_elements(tagIndex)
dt = dt[0:numTimes-1]

nonZero = where(particles.tag GT 0)
smallParticles = particles(nonZero)
; small particles is now a 1 x (numParticles * numTimes) array
print,'smallParticles'
print,smallParticles

; clear out old definition
undefine, particles
; make new definition, where the time step is in the position of the
; previous index
single = {step:0L, tag:0L, block:0L, posx:0.0, posy:0.0, posz:0.0, $
             velx:0.0, vely:0.0, velz:0.0}
particles = replicate(single,numTimes,numParticles)

; sort by tag and assign to new array
numSmall = n_elements(smallParticles)
FOR l=0L, numSmall-1 DO BEGIN 
    iTag = smallParticles[l].tag
    index = where(iTag EQ tagIndex)
    iTime = smallParticles[l].index
    particles[iTime-1,index]=smallParticles[l]
ENDFOR 

; now a function
; return, particles, dt


END ; particles_dump.pro
