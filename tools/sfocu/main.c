#ifdef NEED_MPI
#include <mpi.h>
#endif

int sfocu(int argc, char **argv);
void test_reader(char *filename);

int main(int argc, char **argv){
  int retValue = 0;

#ifdef NEED_MPI
  MPI_Init(&argc, &argv);
#endif


#if 0
  /*test_reader(argv[1], mype, nprocs);*/
#else
  retValue = sfocu(argc, argv);
#endif

#ifdef NEED_MPI
  MPI_Finalize();
#endif

  return retValue;
}

