!!****if* source/Grid/GridMain/Chombo/UG/gr_createDomain
!!
!! NAME
!!
!!  gr_createDomain
!!
!!
!! SYNOPSIS
!!
!!  gr_createDomain()
!!                  
!!
!!
!! DESCRIPTION
!!
!! Creates the Uniform Grid domain. The Uniform Grid expects to be
!! given the processor grid at runtime in the form of runtime parameters
!! iprocs,jprocs and  kprocs, where iprocs*jprocs*kprocs is the
!! number of processor in the run. If UG is running in fixed
!! blocksize mode, then the blocksizes NXB, NYB and NZB are specified
!! at compile time. One block is placed on each processor, so that
!! the global domain size is <NXB*iprocs, NYB*jprocs, NZB*kprocs>.
!! However, if UG is running in non-fixed blocksize mode, the
!! blocksize is determined at runtime. In this mode, UG expects to be
!! given the global domain size in form of runtime parameters
!! iGridSize,jGridSize and kGridSize, and the blocksize is <iGridSize
!!/iprocs,jGridSize/jprocs,kGridSize/kprocs>. As in fixedblocksize mode,
!! only one blocks is placed on each processor
!!  
!! This routine also creates directional communicators used in
!! guardcell exchanges, allocates space for storing physical
!! coordinates of grid points, and computes them.
!!
!!  
!! ARGUMENTS
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_createDomain()
  use iso_c_binding
  use RuntimeParameters_interface
  use Driver_interface, ONLY : Driver_abortFlash
  use flash_ftypes, ONLY : flash_ug_info_t
  use chombo_f_c_interface, ONLY : ch_define_uniform_grid
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs, &
       gr_guard, &
       gr_gIndexSize, &
       gr_lIndexSize, gr_kmax, gr_kmin, &
       gr_jmax, gr_jmin, gr_imax, gr_imin, gr_blkBC, &
       gr_domainBC, gr_iguard,gr_jguard,gr_kguard, &
       gr_verbosity, gr_axisNumProcs
  
  implicit none

#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"

  type(flash_ug_info_t) :: flashUGInfo
  integer :: ierr
  integer :: color, key, range_b, range_e
  real, dimension(MDIM) :: lowDomain, highDomain

  integer :: i, j, k, l, m, meshVar, meshMap
  integer, parameter :: Total_strings = &
       NUNK_VARS + &
       NFACE_VARS + &
       NFACE_VARS + &
       NFACE_VARS + &
       NSCRATCH_GRID_VARS + &
       NSCRATCH_CENTER_VARS + &
       NSCRATCH_FACEX_VARS + &
       NSCRATCH_FACEY_VARS + &
       NSCRATCH_FACEZ_VARS
  integer, parameter, dimension(MAX_GRID_DATA_STRUCT_TMP) :: Mesh_types = &
       (/ CENTER, &
       FACEX, &
       FACEY, &
       FACEZ, &
       SCRATCH, & 
       SCRATCH_CTR, & 
       SCRATCH_FACEX, &
       SCRATCH_FACEY, &
       SCRATCH_FACEZ /)
  integer, parameter, dimension(MAX_GRID_DATA_STRUCT_TMP) :: Mesh_variables = &
       (/ NUNK_VARS, &
       NFACE_VARS, &
       NFACE_VARS, &
       NFACE_VARS, &
       NSCRATCH_GRID_VARS, &
       NSCRATCH_CENTER_VARS, &
       NSCRATCH_FACEX_VARS, &
       NSCRATCH_FACEY_VARS, &
       NSCRATCH_FACEZ_VARS /)
  integer, parameter, dimension(MAX_GRID_DATA_STRUCT_TMP) :: Mesh_map_names = &
       (/ MAPBLOCK_UNK, &
       MAPBLOCK_FACES, &
       MAPBLOCK_FACES, &
       MAPBLOCK_FACES, &
       MAPBLOCK_SCRATCH, &
       MAPBLOCK_SCRATCH_CENTER, &
       MAPBLOCK_SCRATCH_FACEX, &
       MAPBLOCK_SCRATCH_FACEY, &
       MAPBLOCK_SCRATCH_FACEZ /)
  character (len=MAX_STRING_LENGTH) :: meshString
  integer(kind=c_int), dimension(Total_strings) :: meshStringLen
  character(kind=c_char), dimension &
       (Total_strings * MAX_STRING_LENGTH) :: meshStrings
  logical :: restart
  
  call RuntimeParameters_get("restart", restart)

  !Pack all strings into meshStrings (interoperable with C) and also
  !keep track of the length of each mesh string inside meshStrings.
  meshStrings(:) = C_NULL_CHAR
  k = 1; m = 1 !Used to access unit-based arrays.
  do i = 1, MAX_GRID_DATA_STRUCT_TMP
     meshVar = Mesh_variables(i)
     meshMap = Mesh_map_names(i)

     do j = 1, meshVar
        call Simulation_mapIntToStr(j, meshString, meshMap)
        meshStringLen(k) = len(trim(meshString))
        do l = 1, meshStringLen(k)
           meshStrings(m) = meshString(l:l)
           m = m + 1
        end do
        k = k + 1
     end do
  end do


  !store local index size for each block
  gr_lIndexSize = gr_gIndexSize/gr_axisNumProcs


  gr_iguard = gr_guard(IAXIS)
  if(NDIM>1)then
     gr_jguard = gr_guard(JAXIS)
  else
     gr_jguard = 0
  end if
  if(NDIM>2)then
     gr_kguard = gr_guard(KAXIS)
  else
     gr_kguard = 0
  endif
  


  gr_blkBC = gr_domainBC
  if(gr_axisMe(IAXIS)/=0)gr_blkBC(LOW,IAXIS)=NOT_BOUNDARY
  if(gr_axisMe(IAXIS)/=(gr_axisNumProcs(IAXIS)-1))&
       gr_blkBC(HIGH,IAXIS)=NOT_BOUNDARY
  
  if(NDIM > 1) then
     if(gr_axisMe(JAXIS)/=0)gr_blkBC(LOW,JAXIS)=NOT_BOUNDARY
     if(gr_axisMe(JAXIS)/=(gr_axisNumProcs(JAXIS)-1))&
          gr_blkBC(HIGH,JAXIS)=NOT_BOUNDARY
  end if
  
  if(NDIM > 2) then
     if(gr_axisMe(KAXIS)/=0)gr_blkBC(LOW,KAXIS)=NOT_BOUNDARY
     if(gr_axisMe(KAXIS)/=(gr_axisNumProcs(KAXIS)-1))&
          gr_blkBC(HIGH,KAXIS)=NOT_BOUNDARY
  end if



  do i = 1, MDIM
     if (gr_domainBC(LOW,i) == PERIODIC) then
        if (gr_domainBC(LOW,i) /= gr_domainBC(HIGH,i)) then
           call Driver_abortFlash("Chombo requires periodicity on both sides")
        end if
     end if
  end do

  !Pack mesh info into flash_ug_info_t structure.
  flashUGInfo % lowDomain(IAXIS) = gr_imin
  flashUGInfo % lowDomain(JAXIS) = gr_jmin
  flashUGInfo % lowDomain(KAXIS) = gr_kmin
  flashUGInfo % highDomain(IAXIS) = gr_imax
  flashUGInfo % highDomain(JAXIS) = gr_jmax
  flashUGInfo % highDomain(KAXIS) = gr_kmax
  flashUGInfo % procGrid(1:MDIM) = gr_axisNumProcs(1:MDIM)
  flashUGInfo % baseDomainSize(1:MDIM) = gr_gIndexSize(1:MDIM)
  flashUGInfo % guardCells(1:MDIM) = gr_guard(1:MDIM)
  flashUGInfo % domainBC(1:MDIM) = gr_domainBC(LOW,1:MDIM)
  flashUGInfo % meshTypes(1:MAX_GRID_DATA_STRUCT_TMP) = &
       Mesh_Types(1:MAX_GRID_DATA_STRUCT_TMP)
  flashUGInfo % meshNumVars(1:MAX_GRID_DATA_STRUCT_TMP) = &
       Mesh_variables(1:MAX_GRID_DATA_STRUCT_TMP)
  flashUGInfo % verbosity = gr_verbosity

  call ch_define_uniform_grid(flashUGInfo, meshStringLen, meshStrings, logical(restart,kind=c_bool))

end subroutine gr_createDomain
