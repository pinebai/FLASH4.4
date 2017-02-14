


typedef int MPI_Comm;

#define MPI_COMM_WORLD 1
#define MPI_SUM 20
#define MPI_DOUBLE_PRECISION 27
#define MPI_INTEGER 28
#define MPI_INTEGER1 1
#define MPI_REAL 26
#define MPI_REAL4 10
#define MPI_REAL8 11
#define MPI_REAL16 0



#define MPI_INTEGER 28
#define MPI_INTEGER 28



#define MPI_MIN 101
#define MPI_MAX 102
#define MPI_SUCCESS 1001

int MPI_Init(int *argc,char ***argv);
void mpi_init(int *ierr);
int MPI_Abort(MPI_Comm comm, int errorcode);
int MPI_Comm_rank(MPI_Comm comm, int *rank);
int MPI_Allreduce(void *sendBuf, void *recvBuf, int count, int datatype, int op, MPI_Comm comm);
int MPI_Reduce(void *sendBuf, void *recvBuf, int count, int datatype, int op, MPI_Comm comm);

