#ifndef IO_MPI_TYPE_H
#define IO_MPI_TYPE_H

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include <mpi.h>
#include <assert.h>

MPI_Datatype io_mpi_type_primitive(const int flashType);

#endif
