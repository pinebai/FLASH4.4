#include "io_ncmpi_nonblocking.h"

#define MAX_REQUESTS 128
FLASH_IO_Request m_requests[MAX_REQUESTS];
int m_statuses[MAX_REQUESTS];
int m_req_number, m_myPE;

#ifdef DEBUG_IO
const int debugIO = 1;
#else
const int debugIO = 0;
#endif

FLASH_IO_Request * io_ncmpi_nonblocking_get_request()
{
  if (m_req_number < MAX_REQUESTS) {
    m_req_number++;
    if (debugIO && m_myPE == MASTER_PE) {
      printf(" [%s]: Returning pointer to request %d\n",
	     __FILE__, m_req_number);
    }
    return &(m_requests[m_req_number-1]);
  } else {
    if (debugIO && m_myPE == MASTER_PE) {
      printf(" [%s]: Returning a null request pointer\n",
	     __FILE__);
    }
    return NULL;
  }
}

void FTOC(io_ncmpi_nonblocking_init)(const int * const pMyPE)
{
  m_myPE = *pMyPE;
  m_req_number = 0;
}

void FTOC(io_ncmpi_nonblocking_complete_requests)(int *ncid)
{
  int err;
  if (m_req_number == 0) return;
  if (debugIO && m_myPE == MASTER_PE) {
    printf(" [%s]: Waiting on %d non-blocking writes\n",
	   __FILE__, m_req_number);
  }
  
  /* This is a processed string that selects the approrpiate ncmpi 
     wait all function */
  err = flash_wait_all(*ncid, m_req_number, m_requests, m_statuses);

  assert(err == NC_NOERR);
  m_req_number = 0;
}

void FTOC(io_ncmpi_nonblocking_finalize)()
{
	/* clean up here, but we need a netcdf identifier to call
	 * ncmpi_wait_all now which means
	 * io_ncmpi_nonblocking_complete_requests needs one, but we're called
	 * in the finalize path so there's none to pass */
}
