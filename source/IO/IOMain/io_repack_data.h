#ifndef IO_REPACK_DATA_H
#define IO_REPACK_DATA_H

#include <stdlib.h>
#include <mpi.h>

int io_repack_data(void *inbuf,
		   int incount,
		   MPI_Datatype intype,
		   void *outbuf,
		   int outcount,
		   MPI_Datatype outtype);

#endif
