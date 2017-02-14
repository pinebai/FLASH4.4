!!***if* source/physics/Hydro/explicit/split/RHD/hy_rhd_soundSpeed2 
!! NAME
!!
!!  hy_rhd_soundSpeed2
!!
!! SYNOPSIS
!!
!!  hy_rhd_soundSpeed2(real (IN)    :: u(NUNK_VARS,n),
!!                     real (IN)    :: h(n),
!!                     real (OUT)   :: a2(n), 
!!                     integer (IN) :: ibeg,
!!                     integer (IN) :: iend, 
!!                     integer (IN) :: n)
!!
!! DESCRIPTION
!!
!!  Compute sound speed squared 
!!  for an ideal (constant gamma) gas
!!
!!
!! ARGUMENTS
!!
!!  u    - 1-D array of primitive quantities
!!  h    - 1-D array of enthalpies
!!  a2   - 1-D array of squared sound speeds
!!  ibeg - begining of array
!!  iend - end of array
!!  n    - number of points
!!
!!***

subroutine hy_rhd_soundSpeed2(u, h, a2, ibeg, iend, n)

  use Hydro_data , ONLY : hy_gamma

  implicit none

#include "constants.h"
#include "Flash.h"
#include "RHD.h"

  !! Argument list ------------------------------------
  integer, INTENT(in)                :: ibeg, iend, n  
  real, DIMENSION(NUNK_VARS,n), INTENT(in) :: u
  real, DIMENSION(n),  INTENT(in)     :: h
  real, DIMENSION(n),  INTENT(out)    :: a2
  !! --------------------------------------------------

  integer   :: i

  do i = ibeg, iend
    a2(i) = hy_gamma*u(PRES_VAR,i)/(h(i)*u(DENS_VAR,i))
  end do 

end subroutine hy_rhd_soundSpeed2
