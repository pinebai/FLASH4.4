function amr_mach, rho, vx, vy, vz, cs

; finds the mean mach number of an amr flow

if (n_elements(cs) eq 0) then cs=1.0

mass=amr_sum(rho)
px=amr_multiply(rho, vx)
py=amr_multiply(rho, vy)
pz=amr_multiply(rho, vz)
vxbar=amr_sum(px)/mass
vybar=amr_sum(py)/mass
vzbar=amr_sum(pz)/mass
amr_free, px
amr_free, py
amr_free, pz
vxcm=amr_subtract(vx, vxbar)
vycm=amr_subtract(vy, vybar)
vzcm=amr_subtract(vz, vzbar)
vxsqr=amr_multiply(vxcm, vxcm)
vysqr=amr_multiply(vycm, vycm)
vzsqr=amr_multiply(vzcm, vzcm)
amr_free, vxcm
amr_free, vycm
amr_free, vzcm
vtmp=amr_add(vxsqr, vysqr)
vsqr=amr_add(vtmp, vzsqr)
amr_free, vtmp
twoke=amr_multiply(rho, vsqr)
amr_free, vsqr
twoketot=amr_sum(twoke)
amr_free, twoke
vdisp=sqrt(twoketot/mass)
return, vdisp/cs

end
