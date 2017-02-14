#include "constants.h"
#include "Flash.h"

module eos_c_interface
  implicit none

#ifdef USE_EOS_C_INTERFACE
  interface
     subroutine eos_tabReadhdf5 &
          (tableName, groupName, temperatures, densities, avg_ionization, ion_press, &
          ele_press, ion_energy, ele_energy, ion_heat_capacity, ele_heat_capacity) &
          bind(c)
       use iso_c_binding, only : c_char, c_long
       integer(c_long), dimension(*), intent(IN) :: temperatures, densities, avg_ionization, &
                                                    ion_press, ele_press, ion_energy, &
                                                    ele_energy, ion_heat_capacity, &
                                                    ele_heat_capacity 
       character(kind=c_char), dimension(*), intent(IN) :: tableName, groupName 
     end subroutine eos_tabReadhdf5 &
  end interface

  interface
     subroutine eos_tabBrowsehdf5 &
          (tableName, groupName, ntemp, ndens) & bind(c)
       use iso_c_binding, only : c_int, c_char
       integer(c_int), intent(IN) :: ntemp, ndens
       character(kind=c_char), dimension(*), intent(IN) :: tableName, groupName
     end subroutine eos_tabBrowsehdf5
  end interface
#endif

end module