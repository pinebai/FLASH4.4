!!****if* source/Grid/GridMain/Chombo/AMR/gr_createDomain
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
!! DESCRIPTION
!!
!! Defines the initial domain
!!  
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "flash_bool.h"
#include "Flash.h"

subroutine gr_createDomain()
  use iso_c_binding
  use Driver_interface, ONLY : Driver_abortFlash
  use flash_ftypes, ONLY : flash_amr_info_t
  use chombo_f_c_interface, ONLY : ch_define_adaptive_grid
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use RuntimeParameters_interface
  use Grid_data, ONLY : gr_guard, &
       gr_kmax, gr_kmin, &
       gr_jmax, gr_jmin, gr_imax, gr_imin, &
       gr_domainBC, gr_iguard,gr_jguard,gr_kguard, &
       gr_gIndexSize, gr_maxBlockSize, lrefine_max, &
       gr_BRMeshRefineFillRatio, gr_BRMeshRefineBufferSize, &
       gr_BRMeshRefineBlockFactor, gr_verbosity, gr_useQuadCFInterp, &
       gr_useFluxCorrect, gr_refRatio, gr_restrictBeforeGhostExchange, &
       gr_scaleFineFluxes
  implicit none
  type(flash_amr_info_t) :: flashAMRInfo
  integer, dimension(MDIM) :: maxBlockSize
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
  integer(c_int), dimension(Total_strings) :: meshStringLen
  character(kind=c_char), dimension &
       (Total_strings * MAX_STRING_LENGTH) :: meshStrings
  integer(c_int) :: quadCFInterp, fluxCorrect, restrictBeforeGhostExchange, &
       scaleFineFluxes
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

  gr_iguard = gr_guard(IAXIS)
  if (NDIM>1) then
     gr_jguard = gr_guard(JAXIS)
  else
     gr_jguard = 0
  end if
  if (NDIM>2) then
     gr_kguard = gr_guard(KAXIS)
  else
     gr_kguard = 0
  endif
  
  do i = 1, MDIM
     if (gr_domainBC(LOW,i) == PERIODIC) then
        if (gr_domainBC(LOW,i) /= gr_domainBC(HIGH,i)) then
           call Driver_abortFlash("Chombo requires periodicity on both sides")
        end if
     end if
  end do

  maxBlockSize(1:MDIM) = 1
  maxBlockSize(1:NDIM) = gr_maxBlockSize

  if (gr_useQuadCFInterp) then
     quadCFInterp = FLASH_TRUE
  else
     quadCFInterp = FLASH_FALSE
  end if

  if (gr_useFluxCorrect) then
     fluxCorrect = FLASH_TRUE
  else
     fluxCorrect = FLASH_FALSE
  end if

  if (gr_restrictBeforeGhostExchange) then
     restrictBeforeGhostExchange = FLASH_TRUE
  else
     restrictBeforeGhostExchange = FLASH_FALSE
  end if

  if (gr_scaleFineFluxes) then
     scaleFineFluxes = FLASH_TRUE
  else
     scaleFineFluxes = FLASH_FALSE
  end if


  !Pack mesh info into flash_amr_info_t structure.
  flashAMRInfo % lowDomain(IAXIS) = gr_imin
  flashAMRInfo % lowDomain(JAXIS) = gr_jmin
  flashAMRInfo % lowDomain(KAXIS) = gr_kmin
  flashAMRInfo % highDomain(IAXIS) = gr_imax
  flashAMRInfo % highDomain(JAXIS) = gr_jmax
  flashAMRInfo % highDomain(KAXIS) = gr_kmax
  flashAMRInfo % BRMeshRefineFillRatio = gr_BRMeshRefineFillRatio
  flashAMRInfo % BRMeshRefineBufferSize = gr_BRMeshRefineBufferSize
  flashAMRInfo % BRMeshRefineBlockFactor = gr_BRMeshRefineBlockFactor
  flashAMRInfo % maxBlockSize(1:MDIM) = maxBlockSize(1:MDIM)
  flashAMRInfo % baseDomainSize(1:MDIM) = gr_gIndexSize(1:MDIM)
  flashAMRInfo % guardCells(1:MDIM) = gr_guard(1:MDIM)
  flashAMRInfo % domainBC(1:MDIM) = gr_domainBC(LOW,1:MDIM)
  flashAMRInfo % meshTypes(1:MAX_GRID_DATA_STRUCT_TMP) = &
       Mesh_Types(1:MAX_GRID_DATA_STRUCT_TMP)
  flashAMRInfo % meshNumVars(1:MAX_GRID_DATA_STRUCT_TMP) = &
       Mesh_variables(1:MAX_GRID_DATA_STRUCT_TMP)
  flashAMRInfo % maxRefineLevel = lrefine_max - 1
  flashAMRInfo % verbosity = gr_verbosity
  flashAMRInfo % quadCFInterp = quadCFInterp
  flashAMRInfo % fluxCorrect = fluxCorrect
  flashAMRInfo % refRatio = gr_refRatio
  flashAMRInfo % restrictBeforeGhostExchange = restrictBeforeGhostExchange;
  flashAMRInfo % scaleFineFluxes = scaleFineFluxes

  call ch_define_adaptive_grid(flashAMRInfo, meshStringLen, meshStrings, logical(restart,kind=c_bool))

end subroutine gr_createDomain
