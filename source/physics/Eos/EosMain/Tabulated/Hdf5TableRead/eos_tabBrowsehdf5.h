#ifndef EOS_TABBROWSEHDF5_H
#define EOS_TABBROWSEHDF5_H

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "mangle_names.h"
#include "constants.h"


void FTOC(eos_tabBrowsehdf5)(const char * tableName,
			     const char * groupName,
			     const int * ntemp, 
			     const int * ndens);

#endif
