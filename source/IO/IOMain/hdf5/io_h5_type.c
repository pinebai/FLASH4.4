#include "io_h5_type.h"
#define N_FLASH_PRIM_TYPES 4

hid_t io_h5_type_hid_primitive(const int flashType)
{
  hid_t hType;

  switch (flashType) {
  case (IO_FLASH_INT):
    hType = H5T_NATIVE_INT;
    break;
  case (IO_FLASH_DOUBLE):
    hType = H5T_NATIVE_DOUBLE;
    break;
  case (IO_FLASH_FLOAT):
    hType = H5T_NATIVE_FLOAT;
    break;
  case (IO_FLASH_CHAR):
    hType = H5T_NATIVE_CHAR;
    break;
  default:
    Driver_abortFlashC("[io_h5_type]: unknown type");
  }
  return hType;
}

int io_h5_type_flash_primitive(const hid_t hType)
{
  hid_t hPrimitiveTypes[N_FLASH_PRIM_TYPES];
  int fPrimitiveTypes[N_FLASH_PRIM_TYPES];
  int i;

  /* Initialize as follows to avoid:
     warning: initializer element is not computable at load time */
  hPrimitiveTypes[0] = H5T_NATIVE_DOUBLE;
  hPrimitiveTypes[1] = H5T_NATIVE_FLOAT;
  hPrimitiveTypes[2] = H5T_NATIVE_INT;
  hPrimitiveTypes[3] = H5T_NATIVE_CHAR;

  fPrimitiveTypes[0] = IO_FLASH_DOUBLE;
  fPrimitiveTypes[1] = IO_FLASH_FLOAT;
  fPrimitiveTypes[2] = IO_FLASH_INT;
  fPrimitiveTypes[3] = IO_FLASH_CHAR;

  for(i=0; i<N_FLASH_PRIM_TYPES; ++i) {
    if (H5Tequal(hPrimitiveTypes[i],hType) > 0) {
      return fPrimitiveTypes[i];
    }
  }
  Driver_abortFlashC("[io_h5_type]: unknown type");
  return -1;
}

hid_t io_h5_type_create_string(const int strLen)
{
  hsize_t hStrLen = (hsize_t) strLen;
  hid_t hType;
  herr_t hErr;

  hType = H5Tcopy(H5T_C_S1);
  assert (hType >= 0);
  hErr = H5Tset_size(hType, hStrLen);
  assert (hErr >= 0);
  return hType;
}

void io_h5_type_free_string(const hid_t hType)
{
  herr_t hErr;
  hErr = H5Tclose(hType);
  assert (hErr >= 0);
}
