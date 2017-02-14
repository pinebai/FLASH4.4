!!****if* source/Simulation/SimulationMain/unitTest/ParticlesRefine/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This routine implements the Strang splitting scheme for time
!! advancement. A single step in the this driver 
!! includes two sweeps, the first one in order XYZ, and
!! the second one in order ZYX. This driver works with directionally
!! split operators only. The routine also controls the regridding of
!! the mesh if necessary and the simulation output.
!!
!!  
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

#include "constants.h"
#include "Flash.h"


subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr,&
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  use Particles_interface, ONLY : Particles_updateGridVar, Particles_initPositions
  use Simulation_data, ONLY : sim_print, sim_minBlks

  implicit none
  
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(MAXBLOCKS)  :: blkList
  integer :: blkCount, i,j,k, b, blkID



  ! for unit test
  logical,save :: perfect = .true.
  character(len=20) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp
  real :: err1, err=0.0
  logical :: dummy, success=.false.

  temp = dr_globalMe
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe
  ! ------------ end of unitTest setup ---------------------------------------
  

  call Particles_updateGridVar(MASS_PART_PROP, PDEN_VAR)

  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  do b = 1, blkCount
     blkID=blkList(b)
     err=0.0
     call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blkID, solnData, CENTER)
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              err1 = solnData(DENS_VAR,i,j,k)-(solnData(DMPS_VAR,i,j,k)+solnData(PDEN_VAR,i,j,k))
              if(err1>err)err=err1
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blkID,solnData,CENTER)
  end do
  perfect=((err==0.0).and.(blkCount .GE. sim_minBlks))

  if (perfect) then
     write(fileUnit,*)"all results conformed with expected values."
     write(*,*)"all results conformed with expected values."
  else
     write(fileUnit,*)"Failure in particles based refine ",err
  end if


  ! ------------------------------- Gravity unitTest output
   close (fileUnit)   ! for Gravity_unitTest
  ! --------------------------------



  return
  
end subroutine Driver_evolveFlash

