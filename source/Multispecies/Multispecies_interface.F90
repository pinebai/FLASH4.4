!!****h* source/Multispecies/Multispecies_interface
!!
!! This is the header file for the multispecies module
!! that defines its public interfaces.
!!***
Module Multispecies_interface
#include "Flash.h"
#include "constants.h"

  interface Multispecies_setProperty
    subroutine Multispecies_setProperty(name, property, value)
      integer, intent(in)           :: name, property
      real, intent(in)              :: value
    end subroutine Multispecies_setProperty
    subroutine Multispecies_setIntegerProperty(name, property, value)
      integer, intent(in)           :: name, property
      integer, intent(in)           :: value
    end subroutine Multispecies_setIntegerProperty
    subroutine Multispecies_setStringProperty(name, property, value)
      integer, intent(in)           :: name, property
      character(len=*), intent(in)  :: value
    end subroutine Multispecies_setStringProperty

    subroutine Multispecies_setRealArrProperty(name, property, value)
      integer, intent(in)           :: name, property
      real, intent(in)              :: value(:)
    end subroutine Multispecies_setRealArrProperty

    subroutine Multispecies_setIntArrProperty(name, property, value)
      integer, intent(in)           :: name, property
      integer, intent(in)           :: value(:)
    end subroutine Multispecies_setIntArrProperty
  end interface

  interface Multispecies_list
    subroutine Multispecies_list(fileUnit)
      integer, intent(in)                :: fileUnit
    end subroutine Multispecies_list
  end interface

  interface Multispecies_init
    subroutine Multispecies_init()
      
    end subroutine Multispecies_init
  end interface

  interface Multispecies_getSumSqr
    subroutine Multispecies_getSumSqr(property, value, weights, speciesMask)
      integer, intent(in)               :: property
      real, intent(out)                 :: value
      real, intent(in), optional        :: weights(:)
      integer, intent(in), optional     :: speciesMask(NSPECIES)
    end subroutine Multispecies_getSumSqr
  end interface

  interface Multispecies_getSumInv
    subroutine Multispecies_getSumInv(property, value, weights, speciesMask)
      integer, intent(in)                   :: property
      real, intent(out)                     :: value
      real, intent(in), optional            :: weights(:)
      integer, intent(in), optional         :: speciesMask(NSPECIES)
    end subroutine Multispecies_getSumInv
  end interface

  interface Multispecies_getSumFrac
    subroutine Multispecies_getSumFrac(property, value, weights, speciesMask)
      integer, intent(in)           :: property
      real, intent(out)             :: value
      real, intent(in), optional    :: weights(:)
      integer, intent(in), optional :: speciesMask(NSPECIES)
    end subroutine Multispecies_getSumFrac
  end interface

  interface Multispecies_getSum
    subroutine Multispecies_getSum(property, value, weights, speciesMask)
      integer, intent(in)               :: property
      real, intent(out)                 :: value
      real, intent(in), optional        :: weights(:)
      integer, intent(in), optional     :: speciesMask(NSPECIES)
    end subroutine Multispecies_getSum
  end interface

  interface Multispecies_getProperty
    subroutine Multispecies_getProperty(name, property, value)
      integer, intent(in)       :: name, property
      real, intent(out)         :: value
    end subroutine Multispecies_getProperty
    subroutine Multispecies_getIntegerProperty(name, property, value)
      integer, intent(in)       :: name, property
      integer, intent(out)      :: value
    end subroutine Multispecies_getIntegerProperty
    subroutine Multispecies_getStringProperty(name, property, value)
      integer, intent(in)       :: name, property
      character(len=*),intent(out):: value
    end subroutine Multispecies_getStringProperty

    subroutine Multispecies_getRealArrProperty(name, property, value)
      integer, intent(in) :: name, property
      real, intent(out)   :: value(:)
    end subroutine Multispecies_getRealArrProperty

    subroutine Multispecies_getIntArrProperty(name, property, value)
      integer, intent(in)  :: name, property
      integer, intent(out) :: value(:)
    end subroutine Multispecies_getIntArrProperty
  end interface

  interface
    subroutine Multispecies_getPropertyVector(property, values)
      integer, intent(in)       :: property
      real, intent(out)         :: values(:)
    end subroutine Multispecies_getPropertyVector
  end interface

  interface Multispecies_getAvg
    subroutine Multispecies_getAvg(property, value, weights, speciesMask)
      integer, intent(in)               :: property
      real, intent(out)                 :: value
      real, intent(in), optional        :: weights(:)
      integer, intent(in), optional     :: speciesMask(NSPECIES)
    end subroutine Multispecies_getAvg
  end interface

  interface Multispecies_finalize
    subroutine Multispecies_finalize()
    end subroutine Multispecies_finalize
  end interface

  interface Multispecies_unitTest
     subroutine Multispecies_unitTest(fileUnit,perfect)
       integer, intent(IN) :: fileUnit
       logical, intent(INOUT) :: perfect
     end subroutine Multispecies_unitTest
  end interface

end Module Multispecies_interface
