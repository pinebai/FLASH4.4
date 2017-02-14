#ifndef IO_H5_TYPE_H
#define IO_H5_TYPE_H

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include <hdf5.h>
#include <assert.h>

hid_t io_h5_type_hid_primitive(const int flashType);
int io_h5_type_flash_primitive(const hid_t hType);
hid_t io_h5_type_create_string(const int strLen);
void io_h5_type_free_string(const hid_t hType);

#endif
