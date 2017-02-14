#ifndef SET_HDF5_C_H
#define SET_HDF5_C_H

#include <hdf5.h>
#include "assert.h"
#include "mangle_names.h"

void FTOC(init_hdf5_c)(int *pFileID);
		      
void FTOC(finalise_hdf5_c)(int *pFileID);

#endif
