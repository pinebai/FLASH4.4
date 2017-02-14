!!****if* source/Driver/DriverMain/Driver_setupParallelEnv
!!
!! NAME
!!
!!  Driver_setupParallelEnv
!!
!! SYNOPSIS
!!
!!  Driver_setupParallelEnv()
!!
!! DESCRIPTION
!!
!!  Initialize the parallel message-passing environment,
!!  including generation of needed communicators for all units
!!  based upon runtime information such as number of mesh copies
!!  or whether directional communicators are needed.
!!
!!  ARGUMENTS
!!  
!!
!!
!!***
#include "constants.h"
#include "Flash.h"
subroutine Driver_setupParallelEnv ()

  use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs, dr_globalComm,&
       dr_meshComm, dr_meshMe, dr_meshNumProcs,&
       dr_meshAcrossComm, dr_meshAcrossMe, dr_meshAcrossNumProcs, &
       dr_axisComm, dr_axisMe, dr_axisNumProcs,&
       dr_meshCopyCount
  use Driver_interface, ONLY : Driver_abortFlash

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none             

  include "Flash_mpi.h"

  integer, parameter :: nonrep_maxlocs(0:NONREP_COUNT) = NONREP_MAXLOCS
  integer, parameter :: nonrep_rpcount_start(1:NONREP_COUNT+1) = NONREP_RPCOUNT_START
  character(len=*), parameter :: nonrep_rpcount_flat = NONREP_RPCOUNT_FLAT
  integer, parameter :: nonrep_namef_start(1:NONREP_COUNT+1) = NONREP_NAMEF_START
  character(len=*), parameter :: nonrep_namef_flat = NONREP_NAMEF_FLAT_LWR
  
  integer ::   error
  integer :: color, key
  integer :: i,j
  integer :: countInComm
  integer :: n
  
  call RuntimeParameters_get("meshCopyCount",dr_meshCopyCount)

#ifdef FLASH_GRID_UG
  call RuntimeParameters_get("iProcs", dr_axisNumProcs(IAXIS))
  call RuntimeParameters_get("jProcs", dr_axisNumProcs(JAXIS))
  call RuntimeParameters_get("kProcs", dr_axisNumProcs(KAXIS))   
  if (dr_globalNumProcs .ne. &
       dr_meshCopyCount * dr_axisNumProcs(IAXIS) * &
       dr_axisNumProcs(JAXIS) * dr_axisNumProcs(KAXIS)) then
     if(dr_globalMe == MASTER_PE)print*,'the processors are ',dr_axisNumProcs,dr_globalnumProcs
     call Driver_abortFlash("[Driver_init] Must set runtime parameters iProcs, jProcs, kProcs, " // &
        "and meshCopyCount so that iProcs*jProcs*kProcs*meshCopyCount equals number of processors")
  end if
#else
  dr_axisNumProcs=1
#endif

  do i=1, NONREP_COUNT
     call RuntimeParameters_get(nonrep_rpcount_flat(nonrep_rpcount_start(i):nonrep_rpcount_start(i+1)-1), n)
     if(dr_meshCopyCount*nonrep_maxlocs(i) < n) then
        call Driver_abortFlash("["//FILE_AT_LINE//"] the parameter " // &
           nonrep_rpcount_flat(nonrep_rpcount_start(i):nonrep_rpcount_start(i+1)-1) // &
           " must not be greater than the number of local unks " // &
           nonrep_namef_flat(nonrep_namef_start(i):nonrep_namef_start(i+1)-1) // &
           " times meshCopyCount.")
     end if
  end do

  !! first make a communicator for group of processors 
  !! that have the whole computational grid
  !! The grid is duplicated on all communicators
  countInComm=dr_globalNumProcs/dr_meshCopyCount

  if((countInComm*dr_meshCopyCount) /= dr_globalNumProcs)&
       call Driver_abortFlash("when duplicating mesh, numProcs should be a multiple of meshCopyCount")
  
  color = dr_globalMe/countInComm
  key = mod(dr_globalMe,countInComm)
  call MPI_Comm_split(dr_globalComm,color,key,dr_meshComm,error)
  call MPI_Comm_split(dr_globalComm,key,color,dr_meshAcrossComm,error)
   
  call MPI_COMM_RANK(dr_meshComm,dr_meshMe, error)
  call MPI_COMM_SIZE(dr_meshComm, dr_meshNumProcs,error)

  call MPI_COMM_RANK(dr_meshAcrossComm,dr_meshAcrossMe, error)
  call MPI_COMM_SIZE(dr_meshAcrossComm, dr_meshAcrossNumProcs,error)
  
  !! Now create the communicators for each dimension
  dr_axisComm = dr_meshComm
#ifdef DEBUG_GRID
  write(6,*)'the communicator is ',dr_axisComm,dr_meshMe
  write(6,*)'the procgrid is',dr_axisNumProcs
#endif

#if NDIM == 2
  
  color = dr_meshMe/dr_axisNumProcs(1)
  key = mod(dr_meshMe,dr_axisNumProcs(1))
  call MPI_Comm_split(dr_meshComm,color,key,dr_axisComm(1),error)
  call MPI_Comm_split(dr_meshComm,key,color,dr_axisComm(2),error)

#elif NDIM==3

  color = dr_meshMe/dr_axisNumProcs(1)
  key   = mod(dr_meshMe,dr_axisNumProcs(1))
  call MPI_Comm_split(dr_meshComm,color,key,dr_axisComm(1),error)


  color = key+dr_axisNumProcs(1)*(color/dr_axisNumProcs(2))
  key = dr_meshMe/dr_axisNumProcs(1)
  key = mod(key,dr_axisNumProcs(2))
  call MPI_Comm_split(dr_meshComm,color,key,dr_axisComm(2),error)

  color = dr_axisNumProcs(1)*dr_axisNumProcs(2)
  key   = dr_meshMe/color
  color = mod(dr_meshMe,color)
  call MPI_Comm_split(dr_meshComm,color,key,dr_axisComm(3),error)
  
#endif


!! Round off the data size to the nearest multiple of number of 
!! processors along each dimension to determine the local and
!! global blkLimitsGC

  do i = 1,NDIM
     call MPI_Comm_rank(dr_axisComm(i),dr_axisMe(i),error)
  end do

  return
end subroutine Driver_setupParallelEnv
