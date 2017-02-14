function amr_diskanalysis, filename

; Copyright Robert Fisher (2004)
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

;------------------------------------------------------------------------
; Options
;------------------------------------------------------------------------

; Output to PS

psoutput = 1

; Print titles

titles = 1

;------------------------------------------------------------------------
; Constants (cgs)
;------------------------------------------------------------------------

G = 6.637e-8
au = 1.5e13
pi = !Pi

;------------------------------------------------------------------------
; Parameters
;------------------------------------------------------------------------

; tff is the central free-fall time

  tff = 1.21e12

; rdisk is the maximum radial extent of the disk. Note that the result
; obtained is insensitive to this value, provided solely it is large enough
; to encompass the high density cutoff region

  rdisk = 100. * au 

; central density for GWWT04 problem in cgs units

  rhoc = 3.e-18

; Density factor (below central denity) cutoff for disk

  densefac = 1.e6

; Isothermal soundspeed and EOS parameters

  rhostiff = 3.33e4 * rhoc
  ciso = 2.e4 
  gammastiff = 5./3.
  kstiff = ciso * ciso * rhostiff^(1. - gammastiff)

; Number of radial bins to average into.

  nbin = 20 

; Maximum dimensionless variance to allow in least-significant direction.
; For a thin sheetlike configuration, this should be small (< 10%).

  maxvariance = .20
 
;--------------------------------------------------------------------------
; Read in and flatten AMR data onto a single level for analysis.
;--------------------------------------------------------------------------

print, "Reading in AMR data..."

rho = get_amrcomponent (filename, "density")
vx  = get_amrcomponent (filename, "xvel")
vy  = get_amrcomponent (filename, "yvel")
vz  = get_amrcomponent (filename, "zvel")

if (rho.ndim ne 3) then begin
  print, 'Error : ndim must be 3.'
  print, 'ndim = ', rho.ndim
  return, -1
endif

print, "AMR data read in. ", rho.maxlevel, " levels represented."  
print, "Data taken at time = ", rho.time / tff

; Determine location of central maximum

maxrho = max_amr (rho, maxidx = maxidx)

print, "maximum density = ", max_amr (rho)

ix0 = maxidx (rho.ndim + 2:rho.ndim + 4)

ix01 = ix0 (0)
ix02 = ix0 (1)
ix03 = ix0 (2)

flatlevel = rho.maxlevel 

idxmax = rho.idxhi[*, flatlevel] 
idxmin = rho.idxlo[*, flatlevel]

xhi    = rho.boxmax
xlo    = rho.boxmin

dx = dblarr (rho.ndim)
dx = (xhi - xlo) / (idxmax - idxmin + 1)

; Define spatial location of maximum

rmax = dx * ix0 + xlo

print, 'rmax = ', rmax

; Define region to excise disk

cutmin = lonarr (rho.ndim)
cutmax = lonarr (rho.ndim)

ncells = round (rdisk / dx) 
cutmin = round (ix0 - ncells) 
cutmax = round (ix0 + ncells)

;print, 'ix0 = ', ix0
;print, 'rdisk / dx  = ', ncells
;print, 'cutmin = ', cutmin
;print, 'cutmax = ', cutmax

; Due to a limitation in amr_flatten, restrict the flattening to levels lower
; than the maximum level of refinement. If values are zero, then we know we will have
; overly restricted ourselves, and so we reset minlevel automatically.

minlevel = flatlevel 
jump1 : print, "Flattening data..."

rhoflat = amr_flatten (rho, flatlevel, idxmin = cutmin, idxmax = cutmax, minlevel = minlevel, /sample)  
vxflat =  amr_flatten (vx, flatlevel,  idxmin = cutmin, idxmax = cutmax, minlevel = minlevel) 
vyflat =  amr_flatten (vy, flatlevel,  idxmin = cutmin, idxmax = cutmax, minlevel = minlevel) 
vzflat =  amr_flatten (vz, flatlevel,  idxmin = cutmin, idxmax = cutmax, minlevel = minlevel) 

; Transform to rest frame determined by mass-weighted COM

vxmean = mean (vxflat * rhoflat) / total (rhoflat)
vymean = mean (vyflat * rhoflat) / total (rhoflat)
vzmean = mean (vzflat * rhoflat) / total (rhoflat)

vxflat = vxflat - vxmean
vyflat = vyflat - vymean
vzflat = vzflat - vzmean

if (min (rhoflat) eq 0) then begin
  print, "WARNING : Flattened density array <= 0."
;  minlevel = minlevel - 1
;  print, "Decreasing min level to ", minlevel, "."
;  dx = dx * rho.refratio
;  goto, jump1
endif

print, "Data flattened."

; Translate velocity matrices into a vector for use with rotation
; matrices.

sizevx = size (vxflat)
s1 = sizevx (1)
s2 = sizevx (2)
s3 = sizevx (3)

vxarray = reform (vxflat, s1 * s2 * s3)
vyarray = reform (vyflat, s1 * s2 * s3)
vzarray = reform (vzflat, s1 * s2 * s3)

vxt     = transpose (vxarray)
vyt     = transpose (vyarray)
vzt     = transpose (vzarray)

varray = [vxt, vyt, vzt]

;----------------------------------------------------------------------
; Use Principal Component Analysis (PCA) to analyze the density 
; distribution and determine the rotation matrix such that
; the disk is situated in the XY plane.  
;----------------------------------------------------------------------

; Sample the flattened density array above a density cutoff

dx0 = dx (0)

jump2:rhocut = max (rhoflat) / densefac
dense = where (rhoflat gt rhocut)

; convert 1d indices into vector spatial coordinates

wheretomulti, rhoflat, dense, ix, iy, iz

rx =  ix * 1.; 
ry =  iy * 1.;
rz =  iz * 1.; 

rxt = transpose (rx)
ryt = transpose (ry)
rzt = transpose (rz)

PCAarray = [rxt, ryt, rzt]

print, 'Doing PCA analysis on array size...', size (PCAarray)

result = pcomp (PCAarray, eigenvalues = eigenvalues, $
          coefficients = coefficients, variances = variances)

eigenvectors = coefficients / rebin (eigenvalues, 3, 3)

eigeninvert = invert (eigenvectors)

vprime = varray ## eigeninvert

rprime = PCAarray ## eigeninvert
rprime = rprime * dx0 

xprime = rprime (0, *)
yprime = rprime (1, *)
zprime = rprime (2, *)

vxprime = vprime (0, *)
vyprime = vprime (1, *)
vzprime = vprime (2, *)

vxflat = reform (vxprime, s1, s2, s3)
vyflat = reform (vyprime, s1, s2, s3)
vzflat = reform (vzprime, s1, s2, s3)

xmax = max (xprime)
xmin = min (xprime)

ymax = max (yprime)
ymin = min (yprime)

zmax = max (zprime)
zmin = min (zprime) 

x0 = (xmax + xmin) / 2.
y0 = (ymin + ymax) / 2.
z0 = (zmin + zmax) / 2.

D1  = (xmax - xmin)
D2  = (ymax - ymin)
D3  = (zmax - zmin)

xprime = xprime - x0
yprime = yprime - y0
zprime = zprime - z0

if (psoutput eq 1) then begin
  psopen, 'xzpos.ps'
endif

D  = max (D2, D3)

if (titles eq 1) then begin
  title = 'XZ Prime Cell Locations'
  xtitle = 'X Prime'
  ytitle = 'Z Prime' 
endif

plotPosition = aspect (1.)

plot, xprime, zprime, psym = 3, xstyle = 1, ystyle = 1, $
      xrange = [- D/2., D/2.], $
      yrange = [- D/2., D/2.], title = title, xtitle = xtitle, ytitle = ytitle, $
      position = plotPosition

D = max (D1, D2)

if (psoutput eq 1) then begin
  psclose
  psopen, 'xypos.ps'
endif

if (titles eq 1) then begin
  title = 'XY Prime Cell Locations'
  xtitle = 'X Prime'
  ytitle = 'Y Prime'
endif

plot, xprime, yprime, psym = 3, xstyle = 1, ystyle = 1, $
      xrange = [- D/2.,  D/2.],$
      yrange = [- D/2.,  D/2.], title = title, xtitle = xtitle, ytitle = ytitle, $
      position = plotPosition  

print, 'Z-Direction Variance = ', variances (2)
  
if ( variances (2) gt maxvariance) then begin
  densefac = densefac / 2.
  print, 'Variance in 3rd direction large. Adjusting to density cutoff = ', densefac
  if (densefac lt 100) then begin
     print, 'Variance fit not converged. Configuration likely not geometrically 2D thin.'
     stop
  endif
  goto, jump2
endif

r = xprime (*) * xprime (*) + yprime (*) * yprime (*)
r = sqrt (r)
r = (r + 2.e-4)  

phi = atan (yprime (*) , xprime (*) )
;phi = phi + pi / 2.
sinphi = sin (phi)
cosphi = cos (phi)

print, 'max r = ', max (r)
print, 'min r = ', min (r)

print, 'Effective resolution = ', max (r) / dx0, ' fine cells'
 
; create a grid to bin particles
   
dr      = rdisk / nbin

vs      = fltarr (nbin + 1)
height  = fltarr (nbin + 1)
Q       = fltarr (nbin + 1)
sigma   = fltarr (nbin + 1)
torbit  = fltarr (nbin + 1)
tinfall = fltarr (nbin + 1)
vrbar    = fltarr (nbin + 1)
vphibar  = fltarr (nbin + 1)
vphi2bar = fltarr (nbin + 1)
sigmavphi = fltarr (nbin + 1)
sigmabar= fltarr (nbin + 1)
omega   = fltarr (nbin + 1)
omegabar= fltarr (nbin + 1)
massdisk = fltarr (nbin + 1)
udisk   = fltarr (nbin + 1)
vdisp   = fltarr (nbin + 1)
gammadisk = fltarr (nbin + 1)
pressuredisk = fltarr (nbin + 1)
rhodisk = fltarr (nbin + 1)
ndisk   = lonarr (nbin + 1)
radial  = fltarr (nbin + 1)
kappa   = fltarr (nbin + 1)

s = size (dense)
Ngas = s (1)

rind    = lonarr (Ngas)
gamma   = fltarr (Ngas)
kstiff  = fltarr (Ngas)
pressure = fltarr (Ngas)
vr      = fltarr (Ngas)
vphi    = fltarr (Ngas)

gamma (*)   = 1.

rind = round (r / dx0)

print, 'Binning...', s (1), ' cells.'

print, 'Averaging over bins...'

rhofit = rhoflat (dense)
vxfit  = vxflat  (dense)
vyfit  = vyflat  (dense)
vzfit  = vzflat  (dense)

v      = sqrt (vxfit * vxfit + vyfit * vyfit)
vz     = sqrt (vzfit * vzfit)

vr     =   vxfit * cosphi + vyfit * sinphi
vphi   = - vxfit * sinphi + vyfit * cosphi

gamma (*) = 1.001

stiff = where (rhofit gt rhostiff)
if (size (stiff, /n_elements) gt 1) then begin
   gamma (stiff) = 5./3.
endif

kstiff = ciso * ciso * rhostiff^(1. - 5./3.)

print, 'min gamma = ', min (gamma)
print, 'max gamma = ', max (gamma)

radial (0) = au

dr = rdisk / nbin

mass = 0.

rasterxmax = 25 
rasterymax = 25 
array = fltarr (rasterxmax, rasterymax )
array (*, *) = 0.

dx = 2. * rdisk / rasterxmax
dy = 2. * rdisk / rasterymax

for i = 0, rasterxmax  - 1  do begin
  for j = 0, rasterymax - 1 do begin

    x0 = -rdisk + i * dx
    y0 = -rdisk + j * dy
  
    indices = where ( (xprime lt x0 + dx) and (xprime gt x0) and $
                      (yprime lt y0 + dy) and (yprime gt y0) )

    rho = min (rhofit)

    if (indices (0) ne -1) then begin
      rho = total (rhofit (indices) )
    endif

    if (rho ne 0.) then begin
      array (i, j) = alog (rho) / alog (10.)  
    endif

  endfor
endfor     

for i = 1, nbin do begin
   ringr = i * dr + 2.e-4 

   indices = where ( (r lt ringr + dx0) and (r gt ringr) )

   radial (i) = ringr

; Sum cell properties within this annulus in preparation
; for angle-averaging.

   if (indices (0) ne -1) then begin
      num = size (indices)
      ndisk (i) = num (1)

      rhodisk (i) = total ( rhofit (indices) )
      vdisp   (i) = total ( v      (indices) )
      vrbar   (i) = total ( vr     (indices) )
      vphibar (i) = total ( vphi   (indices) )
      vphi2bar (i) = total (vphi (indices) * vphi (indices) )
      gammadisk (i) = total (gamma (indices ) )
      num = size (indices)
      ndisk (i) = num (1)

      sigma   (i) = rhodisk (i) * dx0
      sigma   (i) = sigma (i) * dx0^2. / (2. * pi * ringr * dr)
      gammadisk (i) = gammadisk (i) / ndisk (i)
      maxi = i
   endif else begin ; no cells in this bin
     rhodisk (i) = 0.
     vdisp (i)   = 0.
     vrbar (i)   = 0.
     vphibar (i) = 0.
     vphi2bar (i) = 0.
     gammadisk (i) = 0.
     num = 0.
     ndisk (i) = 1.    ; avoid 0 / 0 error
   endelse

   rhodisk (i) = rhodisk (i) / ndisk (i)
   vrbar   (i) = vrbar (i) / ndisk (i)
   vphibar (i) = vphibar (i) / ndisk (i)
   vphi2bar (i) = vphi2bar (i) / ndisk (i)
   omega (i) = abs (vphi (i)) / ringr
;   omega (i) = vdisp (i) / ringr

   sigmavphi (i) = sqrt (vphi2bar (i)  - vphibar (i) * vphibar (i) )
   vs (i) = sqrt (ciso^2 + gammastiff * kstiff * rhodisk (i)^(gammastiff - 1.) )
   sigma (i) = vs (i) * sqrt (2. * rhodisk (i) / G)
   height (i) = vs (i)^2 / (G * sigma (i) ) / ringr   ; height / radius

   mass = mass + 2. * pi * ringr * dr * sigma (i)

endfor

for i = 1, nbin do begin

      ringr = i * dr + 2.e-4

      if (i lt 2) then begin
        domegadlnbin = ringr * (omega (i + 1) - omega (i)) / dr
      endif else if ( (i gt 2) and (i lt nbin - 1)) then begin
        d1           = ringr * (omega (i)     - omega (i - 1)) / dr
        d2           = ringr * (omega (i + 1) - omega (i) ) / dr
        d3           = ringr * (omega (i + 1) - omega (i - 1) ) / (2. * dr)
        if (d1 * d2  gt 0.) then begin  
          slopes = [abs (d1), abs (d2), abs (d3)]
          minslope  = min (slopes)
          sgn = d1  / abs (d1)
          domegadlnbin = minslope * sgn
        endif else begin
          domegadlnbin = 0.
        endelse
      endif else begin
        domegadlnbin = ringr * (omega (i) - omega (i - 1) ) / dr 
      endelse

      if (domegadlnbin lt - 4. * omega (i) ) then begin
        domegadlnbin = 0.
        kappa (i) = omega (i)
      endif else begin
        kappa (i) = sqrt (4. * omega (i)^2 + omega (i) * domegadlnbin)
      endelse

      if ( (sigma (i) gt 0) and (i lt maxi) ) then begin
        Q (i)     = vs (i) *  kappa (i) / (pi * G * sigma (i) )
      endif else begin
        Q (i) = 0.
      endelse 

endfor

print, "minimum rhofit ", min (rhofit)
print, "maximum rhofit ", max (rhofit)
print, "mass disk = ", mass

; Convert radial scale to AU

;radial = radial / au
;dr     = dr     / au

; Convert radial scale to Dimensionless values

radial = radial / (5000. * au)
dr     =     dr / (5000. * au)
 
logr    = alog (radial) / alog (10.)
logrmin = alog (2. * dr) / alog (10.)
logrmax = alog (dr * maxi) / alog (10.)

ymax = alog (max (vs) ) / alog (10.)
ymin = alog (ciso) / alog (10.)

plotPosition = aspect (1.)

;=====================
; Log rho vs. log r
;=====================

if (psoutput eq 1) then begin
  psclose
  psopen, "logrho.ps"
endif

if (titles eq 1) then begin
  title = "Log rho (gm / cm^3) vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Log rho (gm/cm^3)"
endif

plot, logr,  alog (rhodisk) / alog (10.), xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle

;===========================
; Scale height / r vs. log r
;=========================== 

if (psoutput eq 1) then begin
  psclose
  psopen, "height.ps"
endif

if (titles eq 1) then begin
  title = "height / radius vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "height / radius"
endif

plot, logr, height, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle, yrange = [0., 1.]
  
;=====================
; Log Sigma vs. log r
;=====================

if (psoutput eq 1) then begin
  psclose
  psopen, "logsigma.ps"
endif

if (titles eq 1) then begin
  title = "Log Sigma (gm / cm^2) vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Log Sigma (gm/cm^2)"
endif

plot, logr,  alog (sigma) / alog (10.), xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle 

;==================
; Log cs vs. log r
;==================

if (psoutput eq 1) then begin
  psclose
  psopen, "logvs.ps"
endif

if (titles eq 1) then begin
  title = "Log cs (cm/s) vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Log cs (cm/s)"
endif

ymax = alog (max (vs) ) / alog (10.)
ymin = alog (ciso) / alog (10.)

plot, logr,  alog (vs) / alog (10.), xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle, yrange = [ymin, ymax]

;====================
; v / cs vs. log r
;====================

if (psoutput eq 1) then begin
  psclose
  psopen, "machvphi.ps"
endif

if (titles eq 1) then begin
  title = "v_phi / cs vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "v_phi / cs"
endif

plot, logr, vphibar / vs, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle  

;====================
; v / cs vs. log r
;====================

if (psoutput eq 1) then begin
  psclose
  psopen, "sigmavphi.ps"
endif

if (titles eq 1) then begin
  title = "Sigma v_phi / vphi vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Sigma v_phi / vphi"
endif

plot, logr, abs (sigmavphi / vphibar), xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle, /ylog 

;====================
; vr / cs vs. log r
;====================

if (psoutput eq 1) then begin
  psclose
  psopen, "machvr.ps"
endif

if (titles eq 1) then begin
  title = "v_r / cs vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "v_r / cs"
endif

plot, logr, vrbar / vs, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle

;====================
; t_infall / t_orbit
;====================

if (psoutput eq 1) then begin
  psclose
  psopen, "tinfall.ps"
endif

if (titles eq 1) then begin
  title = "t_infall / t_orbit vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "t_infall / t_orbit"
endif

plot, logr, - vphibar / vrbar, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle

;====================
; v / cs vs. log r
;====================

if (psoutput eq 1) then begin
  psclose
  psopen, "machvz.ps"
endif

if (titles eq 1) then begin
  title = "v_z / cs vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "v_z / cs"
endif

plot, logr, vz / vs, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle

;====================
; Log omega vs. log r
;====================

if (psoutput eq 1) then begin
  psclose
  psopen, "logomega.ps"
endif

if (titles eq 1) then begin
  title = "Log Omega (1/s) vs. Log Radius (AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Log Omega (1/s)"
endif

plot, logr,  alog (omega) / alog (10.), xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle

;===========================
; Log kappa /omega vs. log r
;==========================

if (psoutput eq 1) then begin
  psclose
  psopen, "logkappa.ps"
endif

if (titles eq 1) then begin
  title = "Kappa / Omega vs. Log Radius(AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Kappa / Omega"
endif

plot, logr,  kappa / omega, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle, yrange = [1, 2]

;======================
; Q vs log r
;======================

if (psoutput eq 1) then begin
  psclose
  psopen, "Q.ps"
endif

if (titles eq 1) then begin
  title = "Q vs. Log Radius(AU)"
  xtitle = "Log Radius (AU)"
  ytitle = "Q"
endif

plot, logr, Q, xrange = [logrmin, logrmax], position = plotPosition, title = title, xtitle = xtitle, ytitle = ytitle, yrange = [0, 5.]

if (psoutput eq 1) then begin
  psclose
endif

;==========================================
; Return 1D arrays in a struct for plotting
;==========================================

plotStruct = {logr : fltarr (nbin + 1), Q : fltarr (nbin + 1), Omega : fltarr (nbin + 1), kappa : fltarr (nbin + 1), vrbar : fltarr (nbin + 1), vphibar : fltarr (nbin + 1), logsigma : fltarr (nbin + 1), rhodisk : fltarr (nbin + 1), vs : fltarr (nbin + 1), array : fltarr (rasterxmax, rasterymax) }
 
plotStruct.logr = logr
plotStruct.Q    = Q
plotStruct.Omega = Omega
plotStruct.kappa = kappa   
plotStruct.vrbar = vrbar
plotStruct.vphibar = vphibar
plotStruct.logsigma   = alog (sigma) / alog (10.)
plotStruct.rhodisk = rhodisk
plotStruct.vs      = vs
plotStruct.array   = array

return, plotStruct

end
