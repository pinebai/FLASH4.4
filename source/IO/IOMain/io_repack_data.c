#include "io_repack_data.h"

/* This is currently a copy of ncmpii_data_repack */

int io_repack_data(void *inbuf,
		   int incount,
		   MPI_Datatype intype,
		   void *outbuf,
		   int outcount,
		   MPI_Datatype outtype)
{
  int intypesz, outtypesz;
  int packsz;
  void *packbuf;
  int packpos;

  MPI_Type_size(intype, &intypesz);
  MPI_Type_size(outtype, &outtypesz);

  if (incount*intypesz != outcount*outtypesz) {
    /* input data amount does not match output data amount */
    /* NOTE: we ignore it for user responsibility or add error handling ? */

    /* for rescue, guarantee output data amount <= input data amount */
    if (incount*intypesz < outcount*outtypesz)
      outcount = incount*intypesz/outtypesz;
  }

  if (incount == 0)
    return 0;

  /* local pack-n-unpack, using MPI_COMM_SELF */
  MPI_Pack_size(incount, intype, MPI_COMM_SELF, &packsz);
  packbuf = (void *)malloc(packsz);
  packpos = 0;
  MPI_Pack(inbuf, incount, intype, packbuf, packsz, &packpos, MPI_COMM_SELF);
  packpos = 0;
  MPI_Unpack(packbuf, packsz, &packpos, outbuf, outcount, outtype, MPI_COMM_SELF);
  free(packbuf);

  return 0;
}
