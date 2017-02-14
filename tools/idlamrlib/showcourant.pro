pro showcourant, nsteps, step, initstep 

; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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

; Routine to display the running effective Courant number of an AMR  
;   simulation. That is, it computes the CFL factor abs (v) + c_s
;   along all coordinate directions, and evaluates the coarse level CFL  
;
;		dt_0 = CFL dx_0 / [max (abs (v) + c_s)] 
;
;   The ratio of the actual timestep taken to dt_0 is then tabulated.
;
; c_s is computed for either an ideal EOS (default) or a barotropic
; EOS (with some modification).

base="/scratch/scratchdirs/bobf/amr/gravity/_hyperclaw.mf.rich/run/_becloud/"

;===============================================================
; Problem parameters 
;===============================================================

Courant = 0.5

;================================================================
; Barotropic EOS parameters.
; 
; Assumes an EOS of the form P (rho) = rho ciso^2 + K rho^gamma,
; where K = ciso^2 / rhocrit^(2/3).
;
;================================================================

ciso    = 1.89e4
rhoedge = 5.e-20
gamma   = 5./3.
rhocrit = 1.e4 * rhoedge
K       = ciso^2 / rhocrit^(2./3.) * gamma

tff     = 1.e13

oldtime = 0.

step = 1

; loop over steps

for num = initstep, nsteps, step do begin
   exts='000'
   exts=exts+strcompress(string(num),/remove_all)
   exts=strmid(exts,strlen(exts)-4,4)
   f=base+"plt"+exts
   f=strcompress(f,/remove_all)
   pltfile=f

; Read velocity and density fields.

  xvel    = get_amrcomponent (pltfile, "xvel")
  yvel    = get_amrcomponent (pltfile, "yvel")
  zvel    = get_amrcomponent (pltfile, "zvel")
  rho     = get_amrcomponent (pltfile, "density")

; soundspeed for an ideal EOS
  ieng = get_amrcomponent ("pltfile", "ieng")
  cs2 = amr_multiply (ieng, gamma - 1.)
  cs = amr_power (cs2, 0.5)

; soundspeed for a barotropic equation of state. 
; Uncomment following two lines if using a barotropic EOS.

;  cs2     = amr_add (amr_multiply (amr_pow (rho, 2./3.), K), ciso^2)
;  cs      = amr_pow (cs2, 0.5)

  absvx   = amr_abs (xvel)
  courfacx= amr_add (absvx, cs)

  absvy   = amr_abs (yvel)
  courfacy= amr_add (absvy, cs)

  absvz   = amr_abs (zvel)
  courfacz= amr_add (absvz, cs)

  maxcourx= max_amr (courfacx)
  maxcoury= max_amr (courfacy)
  maxcourz= max_amr (courfacz)

  maxcour = max (maxcourx, maxcoury)
  maxcour = max (maxcour, maxcourz)

  dx0     = rho.gridspacing (0, 0)
  dt      = min (Courant * dx0 / maxcour)

  time    = rho.time
  dtrun   = (time - oldtime) / step

  print,  num, ' ', time / tff, ' ', dtrun / dt

  oldtime = time
 
endfor

end
