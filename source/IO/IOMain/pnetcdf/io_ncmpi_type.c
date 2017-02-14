#include "io_ncmpi_type.h"

nc_type io_ncmpi_type_primitive(const int flashType)
{
  nc_type nType;

  switch (flashType) {
  case (IO_FLASH_INT):
    nType = NC_INT;
    break;
  case (IO_FLASH_DOUBLE):
    nType = NC_DOUBLE;
    break;
  case (IO_FLASH_FLOAT):
    nType = NC_FLOAT;
    break;
  case (IO_FLASH_CHAR):
    nType = NC_CHAR;
    break;
  default:
    Driver_abortFlashC("[io_ncmpi_type]: unknown type");
  }
  return nType;
}
