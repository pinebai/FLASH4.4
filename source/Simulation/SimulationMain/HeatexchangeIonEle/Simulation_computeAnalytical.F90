!!****if* source/Simulation/SimulationMain/HeatexchangeIonEle/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_computeAnalytical
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                       real(IN) :: tcurr)
!!
!!
!!
!! DESCRIPTION
!!
!!  Compute an analytical solution to compare numerical solution with.
!!
!!  This implementation solves simple ion-electron temperature relaxation
!!  for the case of constant ion-electron temperature coupling.
!!
!!  The intial componeent temperatures (for t=0) are assumed to be stored in
!!  UNK variables T0IO and T0EL (and T0R).
!!
!!  The solution is for
!!  
!!   Cv{ion} d/td Tion(r,t) =  - c12 (Tion(r,t) - Tele(r,t))
!!   Cv{ele} d/td Tele(r,t) =    c12 (Tion(r,t) - Tele(r,t))
!!
!!  with initial conditions Tion(r,0), Tele(r,0) given,
!!  where Cv{ion}, Cv{ele} are heat capacities and are assumed constant.
!!
!!  Define
!!   Tm(r,t) := (Cv{ion} Tion(r,t) + Cv{ele} Tele(r,t)) / (Cv{ion} + Cv{ele})
!! (a capacity-weighted avarage temperature), then
!!
!!   d/dt Tm(r,t) = 0
!!
!! so
!!
!!   Tm(r,t) = Tm(r,0) =: Tm(r)
!!
!! Define
!!   y1(r,t) = Tion(r,t) - Tm(r)
!!   y2(r,t) = Tele(r,t) - Tm(r)
!!
!! Then both y1 and y2 have to satisfy the same equation:
!!   d/dt y1(r,t) = - c12 (Cv{ion}+Cv{ele})/(Cv{ion}Cv{ele}) y1(r,t)
!!   d/dt y2(r,t) = - c12 (Cv{ion}+Cv{ele})/(Cv{ion}Cv{ele}) y2(r,t)
!!
!! so the solution is
!!   y1(r,t) = y1(r,0) * f(t),
!!   y2(r,t) = y2(r,0) * f(t),
!! with
!!   f(t) = exp(-c12 (Cv{ion}+Cv{ele})/(Cv{ion}Cv{ele}) t)
!!
!! Substituting this back, we get the analytical solution in the original temperature
!! variables as
!!   Tion(r,t) = Tm(r) + y1(r,t) = Tm(r) + (Tion(0,t)-Tm(r)) * f(t),
!!   Tele(r,t) = Tm(r) + y2(r,t) = Tm(r) + (Tele(0,t)-Tm(r)) * f(t).
!!
!! 
!! The analytical Tion(r,t) is saved in TANA_VAR.
!!
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time
!!
!!
!!***

!!REORDER(4) solnData

subroutine Simulation_computeAnalytical(blockId, tcurr)

  use Simulation_data, ONLY : sim_xCenter, sim_yCenter, &
                              sim_zCenter, sim_rhoInit, &
                              sim_toffset, sim_CvIon, &
                              sim_CvEle, sim_meshMe
  use Heatexchange_data, ONLY : hx_c12
  use Conductivity_interface, ONLY : Conductivity
  use Eos_interface, ONLY : Eos
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  real,    intent(in) :: tcurr

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  real :: xx, yy, zz, x0, y0, z0, r2, toff, pi
  integer :: ivelx, ively, ivelz

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  real,pointer :: solnData(:,:,:,:)

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ,vecLen
  integer, dimension(MDIM) :: axis

  real, dimension(:), allocatable :: eosData
  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone, eintZone, tempZone, &
       temp1Zone, temp2Zone, temp3Zone, &
       temp1ZoneIni, temp2ZoneIni
  real :: CvIon,CvEle, ft, Tm, y20, y2, y10, y1
  real, parameter :: bogusUnusedTemperature = -1e20
  real, save :: bogusUnusedMassFrac(NSPECIES)
  
  logical :: gcell = .true.



  pi = 4. * atan (1.d0)

  ! WE COULD dump some output to stdout listing the parameters.
  ! But currently we don't.
  if (sim_meshMe == MASTER_PE) then

1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)

  endif


  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
  vecLen = 1
  allocate(eosData(vecLen * EOS_NUM))

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCenter, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getBlkPtr(blockId, solnData)
!------------------------------------------------------------------------------

! Coordinates of initial point source

  x0 = sim_xCenter
  y0 = sim_yCenter
  z0 = sim_zCenter

  CvIon = sim_CvIon
  CvEle = sim_CvEle

  ft = exp(-hx_c12 * (CvIon+CvEle)/(CvIon*CvEle) * tcurr)

! Loop over cells in the block.  For each, compute the physical position of 
!  the state using an exact solution.

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     ! get the coordinates of the cell center in the z-direction
     zz = zCenter(k)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        ! get the coordinates of the cell center in the y-direction
        yy = yCenter(j)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx  = xCenter(i)



           toff = sim_toffset  + tcurr
           rhoZone = sim_rhoInit
 
           temp1ZoneIni = solnData(T0IO_VAR,i,j,k)
           temp2ZoneIni = solnData(T0EL_VAR,i,j,k)
           Tm = (CvIon * temp1ZoneIni + CvEle * temp2ZoneIni) / (CvIon+CvEle)
!!$           print*,'Tm is',Tm
           y20 = temp2ZoneIni - Tm
           y2  = y20 * ft
           temp2Zone = y2 + Tm
!!$           print*,'y20,y2 is',y20,y2,temp2Zone
           y10 = temp1ZoneIni - Tm
           y1  = y10 * ft
           temp1Zone = y1 + Tm


           solnData(T1AN_VAR,i,j,k) = temp1Zone
           solnData(T2AN_VAR,i,j,k) = temp2Zone

        enddo
     enddo
  enddo
  deallocate(eosData)
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  call Grid_releaseBlkPtr(blockId, solnData)
     
  return

end subroutine Simulation_computeAnalytical
