!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constant/op_initConstant
!!
!! NAME
!!
!!  op_initConstant
!!
!! SYNOPSIS
!!
!!  call op_initConstant ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the constant section of the opacity unit.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initConstant ()

  use Driver_interface,     ONLY : Driver_abortFlash

  use Opacity_data,         ONLY : op_totalSpecies,         &
                                   op_absorptionKind,       &
                                   op_emissionKind,         &
                                   op_transportKind,        &
                                   op_writeOpacityInfo

  use op_constantData,      ONLY : op_initializedConstant,  &
                                   op_absorptionConstant,   &
                                   op_emissionConstant,     &
                                   op_transportConstant

  use op_interface,         ONLY : op_writeConstants

  use op_numericsData,      ONLY : zero

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
      call Driver_abortFlash ('[op_initConstant] ERROR: No species present!')
  end if
!
!
!    ...Allocate the arrays holding the constant opacities (if wanted) for each species.
!
!
  allocate (op_absorptionConstant (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initConstant] ERROR: op_absorptionConstant allocation failed')
  end if

  allocate (op_emissionConstant   (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initConstant] ERROR: op_emissionConstant allocation failed')
  end if

  allocate (op_transportConstant  (1:op_totalSpecies), stat = status)

  if (status > 0) then
      call Driver_abortFlash ('[op_initConstant] ERROR: op_transportConstant allocation failed')
  end if
!
!
!    ...Initialize the arrays with zeros.
!
!
  op_absorptionConstant = zero
  op_emissionConstant   = zero
  op_transportConstant  = zero
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

     if (absorptionKind == OP_CONSTANT) then
         op_absorptionKind     (species) = OP_CONSTANT
         write(rtpar,'(3a)') "op_", trim(spec_str), "AbsorbConstant"
         call RuntimeParameters_get(rtpar, op_absorptionConstant(species))
         if(op_absorptionConstant(species) < 0.0) &
              call Driver_abortFlash("[op_initConstant] Error: constant absorption opacity not set")
     end if

     if (emissionKind == OP_CONSTANT) then
         op_emissionKind     (species) = OP_CONSTANT
         write(rtpar,'(3a)') "op_", trim(spec_str), "EmissConstant"
         call RuntimeParameters_get(rtpar, op_emissionConstant(species))
         if(op_emissionConstant(species) < 0.0) &
              call Driver_abortFlash("[op_initConstant] Error: constant emission opacity not set")
     end if

     if (transportKind == OP_CONSTANT) then
         op_transportKind     (species) = OP_CONSTANT
         write(rtpar,'(3a)') "op_", trim(spec_str), "TransConstant"
         call RuntimeParameters_get(rtpar, op_transportConstant(species))
         if(op_transportConstant(species) < 0.0) &
              call Driver_abortFlash("[op_initConstant] Error: constant transport opacity not set")
     end if

  end do
!
!
!    ...Write out the Opacity constants to file (if requested).
!
!
  if (op_writeOpacityInfo) then
      call op_writeConstants ()
  end if
!
!
!    ...Set initialization status.
!
!
  op_initializedConstant = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initConstant
