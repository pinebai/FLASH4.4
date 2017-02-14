!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/PlaneSolver_blktri/gr_pfftPoissonDirect
!!
!! NAME
!!
!!  gr_pfftPoissonDirect
!!
!!
!! SYNOPSIS
!!
!!   gr_pfftPoissonDirect(integer(IN) :: iDirection, 
!!                      integer(IN) :: solveflag, 
!!                      integer(IN) :: inSize, 
!!                      integer(IN) :: localSize(MDIM), 
!!                      integer(IN) :: globalSize(MDIM), 
!!                      integer(IN) :: transformType(MDIM)
!!                      real(IN)    :: inArray(:)
!!                      real(OUT)   :: outArray(:))
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the 
!!   fft based method where at most one dimension can be 
!!   periodic, the Neuman solution must happen in a plane
!!
!!
!! ARGUMENTS
!!
!!   iDirection  - direction of the transform, valid values are 
!!                 PFFT_FORWARD and PFFT_INVERSE 
!!   solveflag   - Indicates the solvers to apply
!!                 solveflag==1 => Neuman in Y and Z
!!   inSize      - size of inArray and outArray
!!   localSize   - the local bounds (on myPe) of the original 3D data
!!                 to be transformed
!!   globalSize  - global size of the 3D data to be transformed
!!   transformType - The type if transform to be applied along each
!!                 - of the dimensions
!!   inArray       - single dimension array containing linearized 3D data
!!                   to be transformed
!!   outArray      - single dimension array containing transformed 
!!                   linearized 3D data
!!  
!!
!!
!!
!!***

subroutine gr_pfftPoissonDirect (iDirection, solveflag, inSize, localSize, globalSize, &
                               transformType, inArray, outArray)


  use gr_pfftData, ONLY : pfft_globalLen,pfft_inLen, pfft_outLen, &
                        pfft_midLen,pfft_t1Len, pfft_t2Len, pfft_transformType, &
                        pfft_comm, pfft_me, pfft_procGrid,&
                        pfft_work1,pfft_work2,pfft_localLimits,&
                        pfft_trigIaxis,pfft_trigJaxis,&
                        pfft_trigKaxis,pfft_scale, pfft_workSize,pfft_ndim, pfft_usableProc, pfft_myPE
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, &
                        gr_pfftTranspose, gr_pfftGetLocalLimitsAnytime


!!$  use Grid_interface, ONLY : Grid_pfftMapToInput,Grid_getGlobalIndexLimits,&
!!$       Grid_getBlkIndexLimits,Grid_pfftInit, Grid_pfft,Grid_pfftMapFromOutput, &
!!$       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_pfftFinalize
     

!  use gr_pfftInterface, ONLY : gr_pfftDerivs
!  use gr_interface, ONLY:  gr_findMean

  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe,gr_meshComm, &
       gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
!  use Driver_interface, ONLY : Driver_abortFlash

  use fish

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  
#include "Flash_mpi.h" 


  integer, intent(in)    :: iDirection, solveflag, inSize  
  integer, dimension(MDIM),intent(IN) :: localSize,globalSize,transformType
  real,  dimension(inSize),intent(IN) :: inArray
  real,  dimension(inSize),intent(OUT) :: outArray

  !--------------------------------------------------------------------------
  integer :: numVec , ierr
  real, save :: invScale
  integer,dimension(2,MDIM) :: pfftBlkLimits
  real, save :: DELX,DELY,DELZ,DELXSQ,DELYSQ,DELZSQ
  integer, dimension(MDIM) :: gridOrientation, pfft_meorientation

  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AM, BM, CM
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AN, BN, CN
  TYPE (fishworkspace), SAVE :: W
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: W2

  REAL,DIMENSION(:),ALLOCATABLE, SAVE :: BML
  REAL,DIMENSION(:,:),ALLOCATABLE, SAVE :: RHSP

  real :: DPM
!
!     LOCAL ARRAYS FOR FFT
!
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AK,AL
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: WSAVE,WSAVF

  real, dimension(:,:,:), allocatable :: tmp3DArray

  !NOTE !!!!!!!
  !We used to have routines named gr_pfft1Dto3D and gr_pfft3Dto1D 
  !which have since been deleted.  There may still be references 
  !to these subroutines in old commented out code.  I don't want
  !to remove this old commented out code in case it is of any use.
  !In place of these routines we now use a local "reshape" call.
  !Last time I tested, I found this code to be broken.  It has not 
  !been touched for a while!  Any useful code has already been 
  !incorporated into DirectSolver directory.

  integer, save :: NP,MP,N,M,L,NWK,IERROR
  integer K,I,J,JL,LL

  logical, save :: firstcall = .true.
  !=========================================================================================


  ! Initialization of solver arrays. Depends on solveflag.
  ! ....
  if (firstcall) then

     invScale = 1.0

!-----------------------------------------------------------------------
!     CELL SIZE
!-----------------------------------------------------------------------
     DELX = (gr_imax - gr_imin)/real(globalSize(1))
     DELY = (gr_jmax - gr_jmin)/real(globalSize(2))
     DELZ = (gr_kmax - gr_kmin)/real(globalSize(3))

     DELXSQ = DELX**2.
     DELYSQ = DELY**2.
     DELZSQ = DELZ**2.

!-----------------------------------------------------------------------
!     DIMENSIONS
!-----------------------------------------------------------------------

     N = globalSize(2) !NYG-2 ! Total Number of Points in Y
     M = globalSize(3) !NZ-2  ! Total Number of Points in Z
     L = globalSize(1) !NX-2  ! Total Number of Points in X
    
!     write(*,*) 'In dimensions',N,M,L,solveflag

!-----------------------------------------------------------------------
!     DIMENSION OF ARRAY FOR WORK SPACE AND SAVED ARRAYS
!-----------------------------------------------------------------------
     IF(solveflag==1) THEN  ! Allocate For Blktri Solution

        NP =1 ! Neuman in Y 
        MP =1 ! Neuman in Z

     ELSEIF(solveflag==2) THEN

        NWK=4*N+(10+INT(LOG10(REAL(N))/LOG10(2.)))*M
        ALLOCATE(W2(NWK))

     ENDIF

     ALLOCATE(AK(L))
     ALLOCATE(AM(M),BM(M),CM(M))
!-----------------------------------------------------------------------
!     MODIFIED WAVE NUMBERS
!-----------------------------------------------------------------------
       
     DO K=1, L/2
        AK(K) =  2.*PI*real(K-1)
     ENDDO

     DO K=L/2+1, L
        AK(K) = -2.*PI*real(L-K+1)
     ENDDO

     AK(1:L) = 2.*(1.-COS(AK(1:L)/REAL(L)))/DELXSQ

!-----------------------------------------------------------------------
!     CREATE MATRIX FOR BLKTRI 
!-----------------------------------------------------------------------
     AM(1:M) = 1/DELZSQ !AP(2:M+1)*AU(1:M  )
     CM(1:M) = 1/DELZSQ !AP(2:M+1)*AU(2:M+1)

     IF(MP==1) THEN  
        AM(1)=0.
        CM(M)=0.
     ENDIF

     BM = - AM - CM

     ALLOCATE(BML(M))
     ALLOCATE(RHSP(M,N))

     IF(solveflag==1) THEN

        ALLOCATE(AN(N),BN(N),CN(N))

        AN(1:N) = 1/DELYSQ !CP(2:N+1)*CW(1:N  )
        CN(1:N) = 1/DELYSQ !CP(2:N+1)*CW(2:N+1)

        IF(NP==1) THEN 
           AN(1)=0.
           CN(N)=0.
        ENDIF

        BN = - AN - CN
!
!     BLKTRI NEEDS TO BE CALLED ONCE TO INITIALIZE THE WORKING ARRAY!
!


        IERROR = 0
        CALL BLKTRI(0,NP,N,AN,BN,CN,MP,M,AM,BM,CM,M,RHSP,IERROR,W)
        IF(IERROR/=0) WRITE(6,*) 'INIT. BLKTRI: IERROR = ',IERROR

!!$        ELSEIF(solveflag==2) THEN
!!$
!!$          AM(:)=AM(:)*DELZSQ
!!$          BM(:)=BM(:)*DELZSQ
!!$          CM(:)=CM(:)*DELZSQ
!!$          AK(:)=AK(:)*DELZSQ
!!$
!!$        ELSEIF(solveflag==3) THEN
!!$
!!$          ALLOCATE(WSAVF(2*N+15))
!!$          ALLOCATE(AL(N))
!!$
!!$!-----------------------------------------------------------------------
!!$!     MODIFIED WAVE NUMBERS
!!$!-----------------------------------------------------------------------
!!$          DO K=1,N/2
!!$            AL(K)= 2.*PI*(K-1)
!!$          ENDDO
!!$          DO K=N/2+1,N
!!$            AL(K)=-2.*PI*(N-K+1)
!!$          ENDDO
!!$          
!!$          AL(1:N)=2.*(1.-COS(AL(1:N)/REAL(N)))/DELZSQ
!!$
!!$          CALL RFFTI(N,WSAVF)

     ENDIF

  invScale = 1.0
  firstcall = .false.

  endif



  ! Forward Sweep and Solution:
  IF (iDirection .eq. PFFT_FORWARD) THEN

!!$     write(*,*) 'Before fft..',pfft_globalLen(IAXIS),pfft_scale(IAXIS)

     ! A this point I have the data distributed along pencils in the x direction
     ! Forward transform of only in x (npx = 1, npy=1, npz=np). pfft_ndim==1,  direction==PFFT_FORWARD 
     ! call Grid_pfft(PFFT_FORWARD,inArray,tranArray)
     
!!$     write(*,*) 'Dumping inarray'
!!$     write(*,*)  inArray
!!$     stop         

!!$     write(*,*)'TRANSFORMATION=', pfft_transformType(IAXIS)

     numVec = pfft_inLen(JAXIS)*pfft_inLen(KAXIS)

     call gr_pfftDcftForward(inArray,pfft_work1,pfft_trigIaxis,&
                          pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                          numVec,pfft_transformType(IAXIS),pfft_scale(IAXIS))


        ! Check pz == 1
        ! If pz = 1  (px is also 1 always) and solver BLKTRI => Use Block Tridiagonal
        ! -- ------  --------------------  --- -------------------------------------- 
        ! Need to pass the data to 3D arrays. Go from 1D pfft_work2, 
!!$        call gr_pfft1Dto3D(pfft_work1, 'xyz', pfft_inLen, .false., p3DArray)

!!$        do I = 1 , pfft_inLen(IAXIS)
!!$           write(*,*) transpose(p3DArray(I,:,:))  !transpose
!!$           pause
!!$        enddo
!!$        deallocate(p3DArray)
!!$        nullify(p3DArray)        



  
!!$     do i = 1,36
!!$         write(*,'(I4,6g16.8)') i,pfft_work1((i-1)*6+1:(i)*6)
!!$         pause
!!$     enddo
!!$
!!$     numvec = pfft_inLen(JAXIS)*pfft_inLen(KAXIS)*pfft_inLen(IAXIS)
!!$     write(*,*) 'NUMVEC=',numvec
!!$     write(*,*) 'Sum of Reals =',sum(abs(pfft_work1(1:numvec:2)))
!!$     write(*,*) 'Sum of Imaginary =',sum(abs(pfft_work1(2:numvec:2)))

!!$     write(*,*) 'Dumping pfft_work1',minval(inArray),maxval(inArray)
!!$     write(*,*) 'Dumping pfft_work1',minval(pfft_work1),maxval(pfft_work1)

     ! Transpose to go from xyz to yzx (we end up with npx=np, npy=1, npz=1)
!!$     write(*,*) 'Real or Complex=',pfft_transformType(JAXIS)
!!$     write(*,*) 'T1 len =',pfft_t1Len
!!$     write(*,*) 'Midlen len =',pfft_midLen
!!$     write(*,*) 'Size Work1=',size(pfft_work1)
!!$     write(*,*) 'Size Work2=',size(pfft_work2) 
  
     call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work1,&           !pfft_transformType(JAXIS)
                           pfft_work2,pfft_inLen,pfft_midLen,&
                           pfft_procGrid(JAXIS),pfft_comm(JAXIS))

!!$     write(*,*) 'Sum of Reals =',sum(abs(pfft_work2(1:numvec:2)))
!!$     write(*,*) 'Sum of Imaginary =',sum(abs(pfft_work2(2:numvec:2)))
!!$
!!$     stop

!!$     do i = 1,6
!!$         write(*,*) pfft_work2((i-1)*36+1:(i)*36)
!!$         pause
!!$     enddo

!!$        ! Check pz == 1
!!$        ! If pz = 1  (px is also 1 always) and solver BLKTRI => Use Block Tridiagonal
!!$        ! -- ------  --------------------  --- -------------------------------------- 
!!$        ! Need to pass the data to 3D arrays. Go from 1D pfft_work2, 
!!$        call gr_pfft1Dto3D(pfft_work2, 'yzx', pfft_midLen, .false., p3DArray)
!!$
!!$        do I = 1 , pfft_midLen(IAXIS)
!!$           write(*,*) transpose(p3DArray(:,:,I))  !transpose
!!$           pause
!!$        enddo
!!$        deallocate(p3DArray)
!!$        nullify(p3DArray)    






     call MPI_BARRIER(gr_meshComm,ierr)

     if (solveflag .eq. 1) then   ! isplane 

        ! Check pz == 1
        ! If pz = 1  (px is also 1 always) and solver BLKTRI => Use Block Tridiagonal
        ! -- ------  --------------------  --- -------------------------------------- 
        ! Need to pass the data to 3D arrays. Go from 1D pfft_work2, 
        allocate(tmp3DArray &
             (pfft_midLen(IAXIS), pfft_midLen(JAXIS), pfft_midLen(KAXIS)))
        tmp3DArray = reshape(pfft_work2, pfft_midLen)

        gridOrientation(IAXIS) = JAXIS
        gridOrientation(JAXIS) = KAXIS     
        gridOrientation(KAXIS) = IAXIS

        pfft_meorientation(IAXIS) = IAXIS
        pfft_meorientation(JAXIS) = KAXIS
        pfft_meorientation(KAXIS) = JAXIS
        do i = 1, NDIM
        call gr_pfftGetLocalLimitsAnytime(i,gridOrientation(i),pfft_meorientation(i),pfft_midLen, &
                                            PFFT_PCLDATA_REAL,pfftBlkLimits)
        end do

!!$        write(*,*) 'After Transform and Trnaspose, Processor', pfft_myPE, &
!!$          "has coordinates:", pfft_me, "and grid orienation:", gridOrientation, &
!!$          "the shape of the grid:", pfft_midLen, &
!!$          'portion of grid (IAXIS):', pfftBlkLimits(:,IAXIS), &
!!$          'portion of grid (JAXIS):', pfftBlkLimits(:,JAXIS), &
!!$          'portion of grid (KAXIS):', pfftBlkLimits(:,KAXIS) !, &



        ! BLKTRI BY YZ planes
        ! The solution of this solve by planes is dumped back to tmp3DArray
        LL = pfft_midLen(3)  ! Number of Planes in X as tmp3DArray(y, z, x)
        DO JL = 1, LL        ! LOCAL INDEX

           J= pfftBlkLimits(LOW,KAXIS) + JL      !LL*MYRANK+JL    ! GLOBAL INDEX   !!!!!!! Warning, fix thisss !!!!
  

!!$           write(*,*) ' '
!!$           write(*,*) 'Value of Plane J=',pfft_myPe,J
!!$ 
!!$           call MPI_BARRIER(gr_meshComm,ierr)
!!$           call MPI_FINALIZE(ierr)
!!$           stop    


           !-----------------------------------------------------------------------
           !    CHANGE THE DIAGONAL COEFFICIENTS
           !-----------------------------------------------------------------------
           BML(1:M) = BM(1:M) - AK(J/2+1) !**2.

           !    ASSIGN THE RHS
           DO K=1,N ! Y Direction
              DO I=1,M ! Z Direction
                 RHSP(I,K) = tmp3DArray(K,I,JL)   ! DUMMYT(I+(JL-1)*M+(K-1)*M*LL) (K,I,JL)
              ENDDO
           ENDDO

!!$           write(*,*) 'RHSP 1st  plane, I (Z dir)='
!!$           WRITE(*,*) RHSP(:,:)
!!$           !stop
!!$           pause 
          

!!$           write(*,*) 'NP,MP=',NP,MP
!!$           write(*,*) 'N, M=',N,M
!!$           write(*,*) 'AN =',AN
!!$           write(*,*) 'BN =',BN
!!$           write(*,*) 'CN =',CN
!!$           write(*,*) 'AM =',AM
!!$           write(*,*) 'BML =',BML
!!$           write(*,*) 'CM =',CM        

           !-----------------------------------------------------------------------
           !    CALL THE SOLVER
           !-----------------------------------------------------------------------
           CALL BLKTRI(1,NP,N,AN,BN,CN,MP,M,AM,BML,CM,M,RHSP,IERROR,W)

           IF(IERROR/=0) WRITE(6,*) ' J,IERROR = ',J,IERROR

           !-----------------------------------------------------------------------
           !    FIX THE SUM OF PRESSURE CORRECTION OF WHOLE FIELD TO ZERO
           !-----------------------------------------------------------------------
           IF(J==1) THEN
              DPM=SUM(RHSP)/REAL(M*N)
              RHSP(:,:)=RHSP(:,:)-DPM
           ENDIF

           !    DUMP THE SOLUTION INTO tmp3DArray
           DO K=1,N
              DO I=1,M
                 tmp3DArray(K,I,JL) = RHSP(I,K) 
                 !DUMMYT(I+(JL-1)*M+(K-1)*M*LL) = RHSP(I,K)
                 !DUMMY2(I,JL,K) = RHSP(I,K)
              ENDDO
           ENDDO

!!$           write(*,*) 'RHSP 1st  plane, I (Z dir)='
!!$           WRITE(*,*) tmp3DArray(:,:,JL)
!!$           pause 


        ENDDO

        ! Need to pass the data from 3D array to 1D arry pfft_work2, 
        pfft_work2(1:product(pfft_midLen)) = &
             reshape(tmp3DArray, (/product(pfft_midLen)/))
        deallocate(tmp3DArray)

        ! Transpose to go from yzx to xyz (we end up with npx=1, npy=np, npz=1)
        call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,pfft_work2,&           !pfft_transformType(JAXIS)
                              pfft_work1,pfft_midLen,pfft_inLen,&
                              pfft_procGrid(JAXIS),pfft_comm(JAXIS))


!!$        ! Check pz == 1
!!$        ! If pz = 1  (px is also 1 always) and solver BLKTRI => Use Block Tridiagonal
!!$        ! -- ------  --------------------  --- -------------------------------------- 
!!$        ! Need to pass the data to 3D arrays. Go from 1D pfft_work2, 
!!$        call gr_pfft1Dto3D(pfft_work1, 'xyz', pfft_inLen, .false., p3DArray)

!!$        
!!$        do I = 1 , pfft_inLen(KAXIS)
!!$           write(*,*) 'After inverse transpose=',I
!!$           write(*,*) transpose(p3DArray(:,:,I))  !transpose
!!$           pause
!!$        enddo
!!$        deallocate(p3DArray)
!!$        nullify(p3DArray)        




     else 

        ! Do other solves


     endif

  ELSE     ! DO INVERSE

     if (solveflag .ne. 1) then ! notisplane
  
        ! Do other inverse transforms (related to other solvers)

     endif

     ! Inverse transform of only in x. pfft_ndim==1,  direction==PFFT_INVERSE 
     numVec=pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     call gr_pfftDcftInverse(pfft_work1,outArray,pfft_trigIaxis,&
                             pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                             numVec,pfft_transformType(IAXIS),invScale)
!     outArray=outArray/REAL(L)

  ENDIF


  return
end subroutine gr_pfftPoissonDirect
