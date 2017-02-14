#ifndef IO_NCMPI_TYPE_H
#define IO_NCMPI_TYPE_H

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include <pnetcdf.h>
#include <assert.h>

nc_type io_ncmpi_type_primitive(const int flashType);

#endif
