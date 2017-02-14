!!****h
!!
!! Header File for Dummy MPI Library

       interface MPI_AllReduce
       
       subroutine MPI_AllReduceIntegerArray (sendBuf, recvBuf, count, 
     &  datatype ,reduceOp, mpiComm, ierr)
       INTEGER :: mpiComm
       INTEGER,dimension(:) :: sendBuf
       INTEGER :: count
       INTEGER :: datatype
       INTEGER :: reduceOp
       INTEGER,dimension(:) :: recvBuf
       INTEGER :: ierr
       end subroutine
      
      
       subroutine  MPI_AllReduceRealArray (sendBuf, recvBuf, count, 
     &  datatype, reduceOp, mpiComm, ierr)
       INTEGER :: mpiComm
       REAL,dimension(:) :: sendBuf
       INTEGER :: count
       INTEGER :: datatype
       INTEGER :: reduceOp
       REAL,dimension(:) :: recvBuf
       INTEGER :: ierr
       end subroutine

       subroutine MPI_AllReduceInteger (sendBuf, recvBuf, count, 
     &              datatype,reduceOp, mpiComm, ierr)
       INTEGER :: mpiComm
       INTEGER :: sendBuf
       INTEGER :: count
       INTEGER :: datatype
       INTEGER :: reduceOp
       INTEGER :: recvBuf
       INTEGER :: ierr
       end subroutine
      
      
       subroutine  MPI_AllReduceReal(sendBuf, recvBuf, 
     &		count, datatype,reduceOp, mpiComm, ierr)
       INTEGER :: mpiComm
       REAL :: sendBuf
       INTEGER :: count
       INTEGER :: datatype
       INTEGER :: reduceOp
       REAL :: recvBuf
       INTEGER :: ierr
       end subroutine
      
      
       end interface

       INTEGER, PARAMETER :: MPI_PROC_NULL = -1
       INTEGER, PARAMETER :: MPI_ANY_TAG = -1
       INTEGER, PARAMETER :: MPI_ANY_SOURCE = -2

       INTEGER, PARAMETER :: MPI_Comm_World = 1
       INTEGER, PARAMETER :: MPI_SOURCE = 2
       INTEGER, PARAMETER :: MPI_TAG = 3
       INTEGER, PARAMETER :: MPI_Status_Size = 4

       INTEGER, PARAMETER :: MPI_SUM = 20
       INTEGER, PARAMETER :: MPI_LOGICAL = 25
       INTEGER, PARAMETER :: MPI_REAL = 26
       INTEGER, PARAMETER :: MPI_DOUBLE_PRECISION = 27
       INTEGER, PARAMETER :: MPI_INTEGER = 28
       INTEGER, PARAMETER :: MPI_INTEGER1 = 1
       INTEGER, PARAMETER :: MPI_INTEGER2 = 4
       INTEGER, PARAMETER :: MPI_INTEGER4 = 6
       INTEGER, PARAMETER :: MPI_INTEGER8 = 13
       INTEGER, PARAMETER :: MPI_INTEGER16 = 0
       INTEGER, PARAMETER :: MPI_REAL4 = 10
       INTEGER, PARAMETER :: MPI_REAL8 = 11
       INTEGER, PARAMETER :: MPI_REAL16 = 0

       INTEGER, PARAMETER :: MPI_2DOUBLE_PRECISION = 33

       INTEGER, PARAMETER :: MPI_MIN = 101
       INTEGER, PARAMETER :: MPI_MAX = 102
       INTEGER, PARAMETER :: MPI_LOR = 106
       INTEGER, PARAMETER :: MPI_MINLOC = 110
       INTEGER, PARAMETER :: MPI_MAXLOC = 111
       INTEGER, PARAMETER :: MPI_ERRORS_RETURN = 120

      REAL MPI_WTIME



      
      
