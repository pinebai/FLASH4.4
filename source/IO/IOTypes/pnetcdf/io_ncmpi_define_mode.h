#ifndef IO_NCMPI_DEFINE_MODE_H
#define IO_NCMPI_DEFINE_MODE_H

#include <pnetcdf.h>
#include <assert.h>
#include "mangle_names.h"

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

void FTOC(io_ncmpi_define_mode_enddef)(const int * const pFileID);
void FTOC(io_ncmpi_define_mode_redef)(const int * const pFileID);

#endif
