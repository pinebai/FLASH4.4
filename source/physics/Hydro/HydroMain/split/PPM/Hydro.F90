!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro
!!
!! NAME
!!
!!  Hydro
!!
!! SYNOPSIS
!!
!!  Hydro( integer (IN):: blockCount, 
!!         integer (IN):: blockList(blockCount), 
!!         real    (IN):: timeEndAdv, 
!!         real    (IN):: dt, 
!!         real    (IN):: dtOld, 
!!         integer (IN):: sweepOrder )
!!
!! DESCRIPTION
!!   
!!  Perform a 1, 2, or 3D hydro update.  This version handles
!!  directionally split hydro schemes.  The input
!!  parameter sweepOrder determines the ordering of sweep
!!  directions.  For example, in 3d, SWEEP_XYZ means to perform 
!!  1d sweeps first in the x direction, then the y, then z direction, 
!!  while SWEEP_ZYX means to performs the 1d sweeps in the opposite order.
!!
!!  The blockList and blockCount arguments tell this routine on 
!!  which blocks and on how many to operate.  blockList is an 
!!  integer array of size blockCount that contains the local 
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.  
!!  
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList -   array holding local IDs of blocks on which to advance
!!  timeEndAdv -  end time
!!  dt -          timestep
!!  dtOld -       old timestep
!!  sweepOrder -  direction of hydro sweep, can be: SWEEP_XYZ or SWEEP_ZYX
!!                as defined in  constants.h
!!
!!
!!***

subroutine Hydro(         &
                   blockCount, blockList, &
                   timeEndAdv, dt, dtOld, &
                   sweepOrder )
  use hy_ppm_interface, ONLY : hy_ppm_sweep
  use Hydro_data, ONLY:  hy_useHydro,&
                         hy_gravMass, hy_gravMassXYZ, hy_gravMassZYX,&
                         hy_gravMassXZY, hy_gravMassYZX,&
                         hy_gravMassYXZ, hy_gravMassZXY

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld
  integer, INTENT(IN) :: sweepOrder

  if (.NOT. hy_useHydro) return

  !! Comput diffusive grid velocities
  !! We should keep this tag, but, these routines don't really work.
  !! call grdvel(blockID)
  !! 

  !!--------------------------------------------------------------------------
  !!  Perform three sets of one-dimensional hydro sweeps
  !!--------------------------------------------------------------------------

  hy_gravMass = 0.0   ! zero out this step's accumulation of gravitationalaccel*mass

  select case (sweepOrder)

  case (SWEEP_XYZ)

     call hy_ppm_sweep ( &
                        blockCount, blockList,  &
                        timeEndAdv, dt, dtOld,  &
                        SWEEP_X )



     if (NDIM >= 2) call  hy_ppm_sweep(      &
                                       blockCount, blockList, &
                                       timeEndAdv, dt, dtOld, &
                                       SWEEP_Y )
     
     if (NDIM == 3) call  hy_ppm_sweep(      &
                                       blockCount, blockList, &
                                       timeEndAdv, dt, dtOld, &
                                       SWEEP_Z )
     hy_gravMassXYZ = hy_gravMass

  case (SWEEP_XZY)

     call hy_ppm_sweep (      &
                        blockCount, blockList, &
                        timeEndAdv, dt, dtOld, &
                        SWEEP_X )



     if (NDIM == 3) call  hy_ppm_sweep(      &
                                       blockCount, blockList, &
                                       timeEndAdv, dt, dtOld, &
                                       SWEEP_Z )
     if (NDIM >= 2) call  hy_ppm_sweep(      &
                                       blockCount, blockList, &
                                       timeEndAdv, dt, dtOld, &
                                       SWEEP_Y )
     
     hy_gravMassXZY = hy_gravMass
     

  case (SWEEP_ZYX)
     
     if (NDIM == 3) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Z )

     if (NDIM >= 2) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Y )

     call hy_ppm_sweep(      &
                       blockCount, blockList, &
                       timeEndAdv, dt, dtOld, &
                       SWEEP_X )

     hy_gravMassZYX = hy_gravMass

  case (SWEEP_ZXY)
     
     if (NDIM == 3) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Z )

     call hy_ppm_sweep(      &
                       blockCount, blockList, &
                       timeEndAdv, dt, dtOld, &
                       SWEEP_X )

     if (NDIM >= 2) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Y )


     hy_gravMassZXY = hy_gravMass

  case (SWEEP_YXZ)
     

     if (NDIM >= 2) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Y )

     call hy_ppm_sweep(      &
                       blockCount, blockList, &
                       timeEndAdv, dt, dtOld, &
                       SWEEP_X )

     if (NDIM == 3) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Z )
     hy_gravMassYXZ = hy_gravMass

  case (SWEEP_YZX)
     
     if (NDIM >= 2) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Y )

     if (NDIM == 3) call hy_ppm_sweep(      &
                                      blockCount, blockList, &
                                      timeEndAdv, dt, dtOld, &
                                      SWEEP_Z )

     call hy_ppm_sweep(      &
                       blockCount, blockList, &
                       timeEndAdv, dt, dtOld, &
                       SWEEP_X )


     hy_gravMassYZX = hy_gravMass
  end select

end subroutine Hydro
