!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonsTrace/pi_traceODEfunctionData
!!
!! NAME
!!
!!  pi_traceODEfunctionData
!!
!! SYNOPSIS
!!
!!  use pi_traceODEfunctionData
!!  
!! DESCRIPTION
!!
!!  Data module serving as data transmitter to the proton tracing ODE functions
!!  ---------------------------------------------------------------------------
!!   
!!   B(x,y,z)       : magnetic flux density components in (x,y,z)-direction
!!   B(x,y,z)InvC   : magnetic flux density / speed of light components in (x,y,z)-direction
!!   CurlB(x,y,z)   : curl of magnetic flux density components in (x,y,z)-direction
!!   E(x,y,z)       : electric field components in (x,y,z)-direction
!!   Qm             : proton charge per unit mass
!!   (x,y,z)maxCell : the highest (x,y,z)-coordinate of the cell
!!   (x,y,z)minCell : the lowest (x,y,z)-coordinate of the cell
!!  
!! NOTES
!!
!!  All of the variables in here need to have the 'threadprivate' omp attribute!
!!
!!***

Module pi_traceODEfunctionData

  implicit none

  real,    save :: Bx, By, Bz
  real,    save :: BxInvC, ByInvC, BzInvC
  real,    save :: CurlBx, CurlBy, CurlBz
  real,    save :: Ex, Ey, Ez
  real,    save :: Qm
  real,    save :: xmaxCell, ymaxCell, zmaxCell
  real,    save :: xminCell, yminCell, zminCell

  !$omp threadprivate (Bx, By, Bz,                    &
  !$omp                BxInvC, ByInvC, BzInvC,        &
  !$omp                CurlBx, CurlBy, CurlBz,        &
  !$omp                Ex, Ey, Ez,                    &
  !$omp                Qm,                            &
  !$omp                xmaxCell, ymaxCell, zmaxCell,  &
  !$omp                xminCell, yminCell, zminCell   )

end Module pi_traceODEfunctionData
