;---------------------
; script commands to work with output from particles_dump.pro

; get the data from the dump file
particles_dump,"/home/lynnreid/3/tEuler/test_ParticlesDump_0000",PARTICLES=euler,DT=dt
particles_dump,"/home/lynnreid/3/tPredictorCorrector/test_ParticlesDump_0000",PARTICLES=predcor,DT=dt

; work with particle number 1 (index 0)
xeuler = euler[*,0].posx
yeuler = euler[*,0].posy
xpredcor = predcor[*,0].posx
ypredcor = predcor[*,0].posy
