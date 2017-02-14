#include <mpi.h>

int pfocu(int argc, char **argv, int mype, int nprocs);
void test_reader(char *filename);

int main(int argc, char **argv){
  int retValue = 0;
  int mype, nprocs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#if 0
  /* test_reader(argv[1], mype, nprocs);*/
#else
  retValue = pfocu(argc, argv, mype, nprocs);
#endif

  MPI_Finalize();

  return retValue;
}

