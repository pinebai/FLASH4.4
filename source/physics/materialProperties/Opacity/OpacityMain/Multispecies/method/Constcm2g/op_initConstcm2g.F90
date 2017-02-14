!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constcm2g/op_initConstcm2g
!!
!! NAME
!!
!!  op_initConstcm2g
!!
!! SYNOPSIS
!!
!!  call op_initConstcm2g ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the constant-cm2g section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initConstcm2g ()

  use Driver_interface,   ONLY : Driver_abortFlash

  use Opacity_data,       ONLY : op_totalSpecies,         &
                                 op_absorptionKind,       &
                                 op_emissionKind,         &
                                 op_transportKind,        &
                                 op_writeOpacityInfo

  use op_constcm2gData,    ONLY : op_initializedConstcm2g,  &
                                  op_absorptionConstcm2g,   &
                                  op_emissionConstcm2g,     &
                                  op_transportConstcm2g

  use op_interface,       ONLY : op_writeConstcm2g

  use op_numericsData,    ONLY : zero
  
  use Simulation_interface, ONLY : Simulation_mapIntToStr

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt

  implicit none

#include "Flash.h"
#include "Opacity.h"
#include "constants.h"

  character (len=20) :: spec_str
  character (len=MAX_STRING_LENGTH) :: str
  character (len=MAX_STRING_LENGTH) :: rtpar

  integer :: species
  integer :: status
  integer :: absorptionKind
  integer :: emissionKind
  integer :: transportKind
!
!
!    ...Safety net.
!
!
  if (op_totalSpecies < 1) then
      call Driver_abortFlash ('[op_initConstcm2g] ERROR: No species present!')
  end if
!
!
!    ...Allocate the arrays holding the constant opacities (if wanted) for each species.
!
!
  allocate (op_absorptionConstcm2g (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initConstcm2g] ERROR: op_absorptionConstcm2g allocation failed')
  end if

  allocate (op_emissionConstcm2g   (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initConstcm2g] ERROR: op_emissionConstcm2g allocation failed')
  end if

  allocate (op_transportConstcm2g  (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initConstcm2g] ERROR: op_transportConstcm2g allocation failed')
  end if
!
!
!    ...Initialize the arrays with zeros.
!
!
  op_absorptionConstcm2g = zero
  op_emissionConstcm2g   = zero
  op_transportConstcm2g  = zero
!
!
!    ...Read in the needed info from the runtime parameters and set
!       the handle arrays.
!
!
  do species = 1,op_totalSpecies
     call Simulation_mapIntToStr(species+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)

     write(rtpar,'(3a)') "op_", trim(spec_str), "Absorb"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,absorptionKind)

     write(rtpar,'(3a)') "op_", trim(spec_str), "Emiss"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,emissionKind)

     write(rtpar,'(3a)') "op_", trim(spec_str), "Trans"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,transportKind)

     if (absorptionKind == OP_CONSTCM2G) then
         op_absorptionKind     (species) = OP_CONSTCM2G
         write(rtpar,'(3a)') "op_", trim(spec_str), "AbsorbConstant"
         call RuntimeParameters_get(rtpar, op_absorptionConstcm2g(species))
         if(op_absorptionConstcm2g(species) < 0.0) &
              call Driver_abortFlash("[op_initConstcm2g] Error: constant absorption opacity not set")
     end if

     if (emissionKind == OP_CONSTCM2G) then
         op_emissionKind     (species) = OP_CONSTCM2G
         write(rtpar,'(3a)') "op_", trim(spec_str), "EmissConstant"
         call RuntimeParameters_get(rtpar, op_emissionConstcm2g(species))
         if(op_emissionConstcm2g(species) < 0.0) &
              call Driver_abortFlash("[op_initConstcm2g] Error: constant emission opacity not set")
     end if

     if (transportKind == OP_CONSTCM2G) then
         op_transportKind     (species) = OP_CONSTCM2G
         write(rtpar,'(3a)') "op_", trim(spec_str), "TransConstant"
         call RuntimeParameters_get(rtpar, op_transportConstcm2g(species))
         if(op_transportConstcm2g(species) < 0.0) &
              call Driver_abortFlash("[op_initConstcm2g] Error: constant transport opacity not set")
     end if

  end do
!
!
!    ...Write out the Opacity constants to file (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeConstcm2g ()
  end if
!
!
!    ...Set initialization status.
!
!
  op_initializedConstcm2g = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initConstcm2g
