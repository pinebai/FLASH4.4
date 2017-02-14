#include "io_mpi_type.h"

MPI_Datatype io_mpi_type_primitive(const int flashType)
{
  MPI_Datatype mpiType;

  switch (flashType) {
  case (IO_FLASH_INT):
    mpiType = MPI_INT;
    break;
  case (IO_FLASH_DOUBLE):
    mpiType = MPI_DOUBLE;
    break;
  case (IO_FLASH_FLOAT):
    mpiType = MPI_FLOAT;
    break;
  case (IO_FLASH_CHAR):
    mpiType = MPI_CHAR;
    break;
  default:
    Driver_abortFlashC("[io_mpi_type]: unknown type");
  }
  return mpiType;
}
