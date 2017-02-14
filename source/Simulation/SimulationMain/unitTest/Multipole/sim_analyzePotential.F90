!!****if* source/Simulation/SimulationMain/unitTest/Multipole/sim_analyzePotential
!!
!! NAME
!!
!!  sim_analyzePotential
!! 
!! SYNOPSIS
!!
!!  call sim_analyzePotential (logical (out) :: perfect)
!!
!! DESCRIPTION
!!
!!  This function analyzes the gravitational MacLaurin potential obtained using
!!  the grid multipole solver. The analytical and numerical gravitational potentials
!!  have been calculated and several additional data is determined and stored away
!!  for eventual visualization.
!!
!!  The following data is generated for each cell:
!!
!!    1) The difference     :    analytical - numerical
!!    2) The ratio          :    analytical / numerical
!!    3) The absolute error : | (analytical - numerical) / analytical |
!!
!!  The maximum of the absolute error is determined on the current processor and
!!  used as a criterion for success or failure. The current error tolerance is
!!  based on errors obtained from previous runs for the current refinement level.
!!  This can be extended to include several refinement levels, each of which would
!!  then have a specific error tolerance which must be determined previously.
!!
!! ARGUMENTS 
!! 
!!  perfect : test status indicator (if true, test ran without errors).
!!
!!***

!!REORDER(4): solnData

subroutine sim_analyzePotential (perfect)

  use Driver_data,       ONLY : dr_globalMe

  use Grid_interface,    ONLY : Grid_getBlkIndexLimits, &
                                Grid_getBlkPtr,         &
                                Grid_getListOfBlocks,   &
                                Grid_releaseBlkPtr

  use Simulation_data,   ONLY:  sim_passTolerance

  implicit none

# include "constants.h"
# include "Flash.h"

  logical, intent(out) :: perfect

  integer :: ib,ie,jb,je,kb,ke,i
  integer :: blkCount
  integer :: blockID
  real    :: potError, potErrorMax
  real    :: factorMin, factorMax
  real    :: apotAbsMax, gpotAbsMax

  integer :: blkLimits   (1:2,1:MDIM)
  integer :: blkLimitsGC (1:2,1:MDIM)
  integer :: blkList     (1:MAXBLOCKS)

  real, pointer :: solnData (:,:,:,:)


  integer :: ii,jj,kk


!
!
!   ...Get all leaf blocks and test the obtained potential.
!
!  
  call Grid_getListOfBlocks (LEAF, blkList, blkCount)

  potErrorMax =  tiny (0.0)
  apotAbsMax  =  tiny (0.0)
  gpotAbsMax  =  tiny (0.0)
  factorMin   =  1.0E10
  factorMax   = -1.0E10

  do i = 1,blkCount

     blockID = blkList (i)

     call Grid_getBlkPtr         (blockID, solnData)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC)

     ib = blkLimits (LOW, IAXIS)
     ie = blkLimits (HIGH,IAXIS)
     jb = blkLimits (LOW, JAXIS)
     je = blkLimits (HIGH,JAXIS)
     kb = blkLimits (LOW, KAXIS)
     ke = blkLimits (HIGH,KAXIS)
!
!
!   ...Compute error variables: ERRM = difference, ERRD = division, ERRN = normalized
!
!  
     solnData (ERRM_VAR,ib:ie,jb:je,kb:ke) =   solnData (APOT_VAR,ib:ie,jb:je,kb:ke) &
                                             - solnData (GPOT_VAR,ib:ie,jb:je,kb:ke)

     solnData (ERRD_VAR,ib:ie,jb:je,kb:ke) =   solnData (APOT_VAR,ib:ie,jb:je,kb:ke) &
                                             / solnData (GPOT_VAR,ib:ie,jb:je,kb:ke)

     solnData (ERRN_VAR,ib:ie,jb:je,kb:ke) = abs (  solnData (ERRM_VAR,ib:ie,jb:je,kb:ke) &
                                                  / solnData (APOT_VAR,ib:ie,jb:je,kb:ke) )


!    if (maxval (abs (solnData (ERRN_VAR,ib:ie,jb:je,kb:ke))) > 0.0115) then
!        write (*,*) ' blockID = ',blockID
!    end if


!     if (blockID == 74) then
!         do ii = ib,ie
!         do jj = jb,je
!         do kk = kb,ke
!            if (abs (solnData (ERRN_VAR,ii,jj,kk)) > 0.0115) then
!                write (*,*) ' i,j,k,P = ',ii,jj,kk,solnData (GPOT_VAR,ii,jj,kk)
!            end if
!         end do
!         end do
!         end do
!     end if


     poterrorMax = max (potErrorMax, maxval (abs (solnData (ERRN_VAR,ib:ie,jb:je,kb:ke))))
     apotAbsMax  = max (apotAbsMax , maxval (abs (solnData (APOT_VAR,ib:ie,jb:je,kb:ke))))
     gpotAbsMax  = max (gpotAbsMax , maxval (abs (solnData (GPOT_VAR,ib:ie,jb:je,kb:ke))))
     factorMax   = max (factorMax  , maxval (     solnData (ERRD_VAR,ib:ie,jb:je,kb:ke)) )
     factorMin   = min (factorMin  , minval (     solnData (ERRD_VAR,ib:ie,jb:je,kb:ke)) )

     call Grid_releaseBlkPtr (blockID, solnData)

  end do
!
!
!   ...Analyze the maximum absolute error obtained.
!
!  
  if(potErrorMax < sim_passTolerance) then
     perfect=.true.
     write (*,*) 'Proc',dr_globalMe, ' Maximum potential error = ',potErrorMax,' SUCCESS! '
  else
     perfect=.false.
     write (*,*) 'Proc',dr_globalMe, ' Maximum potential error = ',potErrorMax,' FAILURE! '
  end if
!
!
!   ...Ready!
!
!
  return
end subroutine sim_analyzePotential
