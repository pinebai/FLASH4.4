!!****if* source/Multispecies/MultispeciesMain/Multispecies_init
!!
!!
!! NAME
!!
!!  Multispecies_init
!!
!! SYNOPSIS
!!
!!  Multispecies_init() &
!!                    
!!                    
!!
!! DESCRIPTION
!!  
!!  Initializes data structures of the Multispecies unit.
!!  The implementation subroutine first sets property values to UNDEFINED
!!  and then calls Simulation_initSpecies, where the user can implement
!!  setting the properties to meaningful values.
!!
!! ARGUMENTS
!!
!!  
!!  
!!  
!!
!! SEE ALSO
!!
!!   Simulation_initSpecies
!!
!!***!

subroutine Multispecies_init()

  use Multispecies_data, ONLY : ms_Array, ms_globalMe

  use Simulation_interface, ONLY : Simulation_initSpecies, &
                                   Simulation_mapIntToStr

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
                                         RuntimeParameters_mapStrToInt

  use Driver_interface, ONLY : Driver_getMype

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "constants.h"

  integer                           :: i, n
  character (len=20)                :: spec_str
  character (len=MAX_STRING_LENGTH) :: rtpar
  character (len=MAX_STRING_LENGTH) :: str

!! Assign things that every init unit should know
  call Driver_getMype(GLOBAL_COMM, ms_globalMe)

  do i=1, NSPECIES
     ms_Array(i)%name           = i + NPROP_VARS
     ms_Array(i)%numTotal       = UNDEFINED_REAL
     ms_Array(i)%numPositive    = UNDEFINED_REAL
     ms_Array(i)%numNeutral     = UNDEFINED_REAL
     ms_Array(i)%numNegative    = UNDEFINED_REAL
     ms_Array(i)%bindingEnergy  = UNDEFINED_REAL
     ms_Array(i)%adiabaticIndex = UNDEFINED_REAL
     ms_Array(i)%zMin           = UNDEFINED_REAL
     ms_Array(i)%opacityLowTemp = 0.0
     ms_Array(i)%eosType        = UNDEFINED_INT
     ms_Array(i)%eosSubtype     = UNDEFINED_INT
     ms_Array(i)%eosZfreetableFile= MS_UNDEFINED_STRING
     ms_Array(i)%eosEnertableFile= MS_UNDEFINED_STRING
     ms_Array(i)%eosPrestableFile= MS_UNDEFINED_STRING
     ms_Array(i)%eosGroupName   = MS_UNDEFINED_STRING
     ms_Array(i)%numElems       = 0
     ms_Array(i)%zElems(:)      = UNDEFINED_INT
     ms_Array(i)%aElems(:)      = UNDEFINED_REAL
     ms_Array(i)%fractions(:)   = 0.0
  end do

#ifdef SPECIES_SETUPVAR
  ! SPECIES_SETUPVAR will be defined if the "species" setup variable
  ! has been set, as in
  !
  !   ./setup  ...  species=spec1,spec2,...,specN 
  !
  ! When set, runtime parameters will be created
  ! for each property in the multispecies database for each
  ! species. Thus, the user can automatically set these options using
  ! flash.par instead of having to use Simulation_initSpecies.
  ! 
  do i=1, NSPECIES
     call Simulation_mapIntToStr(i+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)     
     
     write(rtpar,'(3a)') "ms_", trim(spec_str), "A"
     call RuntimeParameters_get(rtpar, ms_Array(i)%numTotal)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "Z"
     call RuntimeParameters_get(rtpar, ms_Array(i)%numPositive)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "Neutral"
     call RuntimeParameters_get(rtpar, ms_Array(i)%numNeutral)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "Negative"
     call RuntimeParameters_get(rtpar, ms_Array(i)%numNegative)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "BindEnergy"
     call RuntimeParameters_get(rtpar, ms_Array(i)%bindingEnergy)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "Gamma"
     call RuntimeParameters_get(rtpar, ms_Array(i)%adiabaticIndex)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "ZMin"
     call RuntimeParameters_get(rtpar, ms_Array(i)%zMin)

     write(rtpar,'(3a)') "op_", trim(spec_str), "LowTemp"
     call RuntimeParameters_get(rtpar, ms_Array(i)%opacityLowTemp)

     write(rtpar,'(3a)') "ms_", trim(spec_str), "NumElems"
     call RuntimeParameters_get(rtpar, ms_Array(i)%numElems)

     do n = 1, MS_MAXELEMS
        write(rtpar,'(3a,i0)') "ms_", trim(spec_str), "ZElems_", n
        call RuntimeParameters_get(rtpar, ms_Array(i)%zElems(n))

        write(rtpar,'(3a,i0)') "ms_", trim(spec_str), "AElems_", n
        call RuntimeParameters_get(rtpar, ms_Array(i)%aElems(n))

        write(rtpar,'(3a,i0)') "ms_", trim(spec_str), "Fractions_", n
        call RuntimeParameters_get(rtpar, ms_Array(i)%fractions(n))
     end do

     write(rtpar,'(3a)') "eos_", trim(spec_str), "EosType"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,ms_Array(i)%eosType)

     write(rtpar,'(3a)') "eos_", trim(spec_str), "SubType"
     call RuntimeParameters_get(rtpar, str)
     call RuntimeParameters_mapStrToInt(str,ms_Array(i)%eosSubType)

     write(rtpar,'(3a)') "eos_", trim(spec_str), "TableFile"
     call RuntimeParameters_get(rtpar, ms_Array(i)%eosZfreetableFile)
     call RuntimeParameters_get(rtpar, ms_Array(i)%eosEnertableFile)
     call RuntimeParameters_get(rtpar, ms_Array(i)%eosPrestableFile)

     write(rtpar,'(3a)') "eos_", trim(spec_str), "GroupName"
     call RuntimeParameters_get(rtpar, ms_Array(i)%eosGroupName)

  end do
#endif



  call Simulation_initSpecies()

end subroutine Multispecies_init
