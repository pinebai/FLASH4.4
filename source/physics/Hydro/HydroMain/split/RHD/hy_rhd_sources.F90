!!****if* source/physics/Hydro/HydroMain/split/RHD/hy_rhd_sources
!!
!! NAME
!!
!!   hy_rhd_sources
!!
!! SYNOPSIS
!!   hy_rhd_sources( real   (IN) :: Vc(NUNK_VARS,n),
!!                   real   (IN) :: h(n),
!!                   real   (IN) :: cs2(n),
!!                   real  (OUT) :: src(NUNK_VARS,n),
!!                   real   (IN) :: r(n),
!!                   real   (IN) :: dr(n),
!!                   integer(IN) :: ibeg,
!!                   integer(IN) :: iend,
!!                   integer(IN) :: n,
!!                   integer(IN) :: dir,
!!                   integer(IN) :: eqnType)
!!
!!
!! DESCRIPTION
!!
!!   Compute source terms for the  relativistic equations.
!!   The argument eqnType determine whether the equations are in
!!   conservative or primitive form
!!
!! ARGUMENTS
!! 
!!   Vc      -    Array of centred values
!!   h       -    specific enthalpy   (needed in primitive form)
!!   cs2     -    sound speed squared (needed in primitive form) 
!!   src     -    source term
!!   r       -    locations of cell centroids.    
!!   dr      -    mesh spacing 
!!   ibeg    -    initial point for source term evaluation
!!   iend    -    final point for source term evaluation
!!                routine that advances species in Lagrangian fashion
!!   n       -    Size of arrays in the sweep direction
!!   dir     -    Sweep direction IAXIS, JAXIS, or KAXIS
!!   eqnType -    indicates whether equations are primitive or consevered form
!!                choices HY_CONSERVE, HY_PRIMITIVE
!!
!!***

subroutine hy_rhd_sources(Vc, h, cs2, src, r, dr, ibeg, iend, n, dir, eqnType)

  use Hydro_data, ONLY : hy_meshGeom

  implicit none
#include "constants.h"
#include "Flash.h"
#include "RHD.h"  

  !! Argument list -------------------------------------
  integer, INTENT(IN) :: ibeg, iend, n, dir,eqnType
  real, DIMENSION(NUNK_VARS,n), INTENT(OUT) :: src
  real, DIMENSION(NUNK_VARS,n), INTENT(IN)  :: Vc
  real, DIMENSION(n), INTENT(IN)       :: cs2, h, r, dr
  integer          :: i,ivn
  real             :: vel2, delta_2, scrh_1, scrh_2
  !! ---------------------------------------------------

  src = 0.d0

  if(dir==IAXIS) then
     
     if (hy_meshGeom == CYLINDRICAL) then
        if(eqnType==HY_CONSERVE) then
           
           ivn = VELX_VAR + dir - IAXIS
           do i = ibeg, iend
              src(ivn, i) = Vc(PRES_VAR, i)/r(i)
           end do

        else  !! HY_PRIMITIVE
           do i = ibeg, iend
              vel2 = Vc(VELX_VAR, i)*Vc(VELX_VAR,i) + & 
                     Vc(VELY_VAR, i)*Vc(VELY_VAR,i) + &
                     Vc(VELZ_VAR, i)*Vc(VELZ_VAR,i)
              delta_2 = 1.0/(1.0 - vel2*cs2(i))
              
              scrh_1 = Vc(VELX_VAR,i)*delta_2/r(i)
              scrh_2 = scrh_1*cs2(i)*(1.0 - vel2)
              src(DENS_VAR, i) = -Vc(DENS_VAR,i)*scrh_1
              src(VELX_VAR, i) =  Vc(VELX_VAR,i)*scrh_2
              src(VELY_VAR, i) =  Vc(VELY_VAR,i)*scrh_2
              src(VELZ_VAR, i) =  Vc(VELZ_VAR,i)*scrh_2
              src(PRES_VAR, i) = -Vc(DENS_VAR,i)*h(i)*cs2(i)*scrh_1
           end do
        end if
        
     else if (hy_meshGeom == SPHERICAL ) then
        if(eqnType==HY_CONSERVE) then        
           ivn = VELX_VAR + dir - IAXIS
           do i = ibeg, iend
              src(ivn, i) = 2.0*Vc(PRES_VAR, i)/r(i)
           end do
        end if
        !! No activity for HY_PRIMITIVE with SPHERICAL
     end if
  end if

end subroutine hy_rhd_sources
