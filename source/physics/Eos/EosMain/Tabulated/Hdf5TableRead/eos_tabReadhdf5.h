#indef EOS_TABREADHDF5_H
#define EOS_TABREADHDF5_H

#include "mangle_names.h"
#include "constants.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOX(x) x
#endif

void FTOC(eos_read_hdf5)(char * tableName,
                         char * groupName,
                         double * temperatures,
                         double * densities,
                         double * avg_ionization,
                         double * ion_press,
                         double * ele_press,
                         double * ion_energy,
                         double * ele_energy,
                         double * ion_heat_capacity,
                         double * ele_heat_capacity,
                         double * entropy)

#endif
