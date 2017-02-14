!!****if* source/physics/Eos/EosMain/eos_fillMapLookup
!!
!! NAME
!!
!!  eos_fillMapLookup
!!
!! SYNOPSIS
!!
!!  call eos_fillMapLookup()
!!
!! DESCRIPTION
!!
!!  Alleviates the expense of the eos_variableMap function call by storing 
!!  eos map data in a module level array.
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Eos_map.h"

subroutine eos_fillMapLookup()
  use Eos_data, ONLY : eos_mapLookup
  implicit none

  !Eventually we will re-use existing constants which will enumerate 1,2,3,4,5.
!!$  integer, parameter, dimension(5) :: StructLookup = (/ &
!!$       CENTER, &
!!$       FACEX, &
!!$       FACEY, &
!!$       FACEZ, &
!!$       SCRATCH &
!!$       /)
  integer :: i, j, k, dataStruct
  integer, external :: eos_variableMap

  eachStruct: do k = 1, MAX_GRID_DATA_STRUCT
     dataStruct = k
     eachDirection: do j = EOS_IN,EOS_OUT
        eachVariable: do i = 1, EOSMAP_NUM_ROLES

           eos_mapLookup(i, j, k) = &
                eos_variableMap(dataStruct, i, (j-1)) 

        end do eachVariable
     end do eachDirection
  end do eachStruct
end subroutine eos_fillMapLookup
