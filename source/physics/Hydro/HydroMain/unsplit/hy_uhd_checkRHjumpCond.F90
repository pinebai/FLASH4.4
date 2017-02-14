!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_checkRHjumpCond
!!
!! NAME
!!
!!  hy_uhd_checkRHjumpCond
!!
!! SYNOPSIS
!!
!!  hy_uhd_checkRHjumpCond(integer(IN)   :: dir,
!!                         real(IN)      :: dens,
!!                         real(IN)      :: velocity(MDIM),
!!                         real(IN)      :: pres,
!!                         real(IN)      :: gamc,
!!                         real(IN)      :: Wp(HY_DENS:HY_EINT),
!!                         real(IN)      :: Wn(HY_DENS:HY_EINT),
!!                         logical(OUT)  :: SWp,SWn)
!!
!! ARGUMENTS
!!
!!  dir      - direction under reconstruction
!!  dens     - density
!!  velocity - velocity fields
!!  pres     - pressure
!!  gamc     - gamc
!!  Wp,Wn    - reconstructed right(plus) and left(negative) Riemann states
!!  SWp,SWn  - logical switch to turn on/off order reduction
!!
!!
!! DESCRIPTION
!!
!!  behind shock |  ahead of shock
!!  (region 2)   |   (region 1)
!!               | 
!!     u2        |     u1 (= shock speed, that is the flow of gas comes to the shock with the shock speed)
!!   <-----      |  <------   
!!     rho2      |    rho1
!!     p2        |     p1
!!               |
!! (shocked gas) |  (undisturbed gas)
!!               |
!!             shock
!!
!!        shock frame (in this shock frame, shock speed = 0)
!!
!! Across the hydro shock, the Rankine-Hugoniot jump relations should satisfy:
!!
!!  (1) rho2/rho1 = [(gamc+1)*M1^2]   / [2+(gamc-1)*M1^2]
!!  (2) u2/u1     = [2+(gamc-1)*M1^2] / [(gamc+1)*M1^2]
!!  (3) p2/p1     = [2*gamc*M1^2 - (gamc-1)] / [gamc+1]
!!
!! Consequences:
!!  (1) M1>=1 (the shock speed (u1) exceeds the sound speed ahead of the shock
!!  (2) u2 <= soundSpeed in region 2 (subsonic behind the shock; supersonic ahead of it)
!!  (3) p2 >= p1 and rho2 >= rho1 (the shock is compressive)
!!  (4) u2 <= u1 and temp2 >= temp1  (the shock slows down the gas and heats it up)
!!  (5) 1<= rho2/rho1 < [gamc+1]/[gamc-1] (the maximum density ratio is [gamc+1]/[gamc-1];
!!                                         the pressure ratio increases with M1^2)
!! In the presence of magnetic fields, the relationship is lightly modified: letting X= rho2/rho1,
!!  (6) u2/u1 = 1/X
!!  (7) the shock is compressive with X>= 1
!!  (8) the effect of magnetic fields is to recude X below its hydro limit
!!  (9) the shock speed (u1) must exceeds the fast magneosonic speed sqrt(C1^2+V_alfven1^2), 
!!      where C1 is sound speed in region1
!!  (10) 1<B2/B1 < [gamc+1]/[gamc-1]
!!
!! REFERENCES
!!
!!  * Kirk, Melrose, Priest, Plasma Astrophysics
!!  * Toro
!!
!!***


Subroutine hy_uhd_checkRHjumpCond(dir,dens,velocity,pres,gamc,Wp,Wn,SWp,SWn)

  implicit none

#include "Flash.h"
#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer, intent(IN)  :: dir
  real, dimension(MDIM), intent(IN) :: velocity
  real, intent(IN) :: dens, pres, gamc
  real, dimension(HY_DENS:HY_EINT), intent(IN) :: Wp,Wn
  logical, intent(OUT) :: SWp,SWn
  !!------------------------------------------------------------------------

  real :: Mach2, gammaRatio, densRatio,gamp,gamm

  SWp = .false.
  SWn = .false.


  gamp = gamc+1.
  gamm = gamc-1.
  gammaRatio = gamp/gamm

  ! Mach2 = Mach^2
  Mach2 = dot_product(velocity(1:NDIM),velocity(1:NDIM))
  Mach2 = Mach2*dens/(gamc*pres)

  if (Mach2 > 1.e-4) then
     if (Wn(HY_DENS) >    gammaRatio*dens) SWn = .true.
     if (   gammaRatio*Wn(HY_DENS) < dens) SWn = .true.
     if (Wp(HY_DENS) >    gammaRatio*dens) SWp = .true.
     if (   gammaRatio*Wp(HY_DENS) < dens) SWp = .true.
  endif

End Subroutine hy_uhd_checkRHjumpCond
