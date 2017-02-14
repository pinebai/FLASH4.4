#ifndef IO_NCMPI_NONBLOCKING_H
#define IO_NCMPI_NONBLOCKING_H

#include <pnetcdf.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "mangle_names.h"

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

/* Test for Pnetcdf >= 1.2.0 */
#if ((PNETCDF_VERSION_MAJOR > 1) ||					\
     (PNETCDF_VERSION_MAJOR == 1 && PNETCDF_VERSION_MINOR >= 2))

typedef int FLASH_IO_Request;
#define flash_wait_all(a,b,c,d) ncmpi_wait_all(a,b,c,d)

#else

/* Version <= 1.1.1 */
typedef NCMPI_Request FLASH_IO_Request;
#define flash_wait_all(a,b,c,d) ncmpi_waitall(b,c)

#endif

FLASH_IO_Request * io_ncmpi_nonblocking_get_request();
void FTOC(io_ncmpi_nonblocking_init)();
void FTOC(io_ncmpi_nonblocking_complete_requests)();
void FTOC(io_ncmpi_nonblocking_finalize)();

#endif
