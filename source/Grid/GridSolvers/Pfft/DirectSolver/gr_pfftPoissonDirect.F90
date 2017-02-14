!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/gr_pfftPoissonDirect
!!
!! NAME
!!
!!  gr_pfftPoissonDirect
!!
!!
!! SYNOPSIS
!!
!!  call gr_pfftPoissonDirect(integer(IN)                :: iDirection,
!!                            integer(IN)                :: solveflag,
!!                            integer(IN)                :: inSize,
!!                            integer(IN)                :: localSize(MDIM),
!!                            integer(IN)                :: globalSize(MDIM),
!!                            integer(IN)                :: transformType(MDIM),
!!                            real(IN) ,dimension(inSize):: inArray(:),
!!                            real(OUT),dimension(inSize):: outArray(:))
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the 
!!   fft based method for periodic and dirichlet problems
!!   on a uniform grid.  
!!   Isolated problems are not supported.
!!
!!
!! ARGUMENTS
!!
!!   iDirection  - direction of the transform, valid values are 
!!                 PFFT_FORWARD and PFFT_INVERSE 
!!   solveflag   - Indicates the boundary conditions and therefore
!!                 the solvers to be applied.
!!                 solveflag==1 => Periodic in X, Neumann in Y and Z
!!                 solveflag==2 => Periodic in X and Y, Neumann in Z
!!                 solveflag==3 => Periodic in X, Y and Z
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


subroutine gr_pfftPoissonDirect (iDirection, solveflag, &
     inSize, localSize, globalSize, &
     transformType, inArray, outArray)
  
  
  use gr_pfftData, ONLY : pfft_globalLen,pfft_inLen, pfft_outLen, &
       pfft_midLen,pfft_t1Len, pfft_t2Len, pfft_transformType, &
       pfft_comm, pfft_me, pfft_procGrid,&
       pfft_work1,pfft_work2,&
       pfft_trigIaxis,pfft_trigJaxis,&
       pfft_trigKaxis,pfft_scale, pfft_workSize,pfft_ndim, &
       pfft_usableProc, pfft_myPE
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, &
       gr_pfftTranspose, gr_pfftGetLocalLimitsAnytime
  
  use Grid_data, ONLY : gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
  use Timers_interface, ONLY : Timers_start, Timers_stop
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
  integer, save, dimension(2,MDIM) :: pfftBlkLimits
  real, save :: DELX,DELY,DELZ,DELXSQ,DELYSQ,DELZSQ
  integer, dimension(MDIM) :: gridOrientation, pfft_meorientation
  
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AM, BM, CM
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AN, BN, CN

  TYPE (fishworkspace), SAVE :: W
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: W2

  REAL,DIMENSION(:),ALLOCATABLE, SAVE :: BML,BMM,CMM,DPZ
  REAL,DIMENSION(:,:),ALLOCATABLE, SAVE :: RHSP

  real :: DPM,DPMaux
!
!     LOCAL ARRAYS FOR FFT
!
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AK,AL,AZ
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: WSAVE,WSAVF

  real, dimension(:,:,:), allocatable :: tmp3DArray

  integer, save :: LP,NP,MP,N,M,L,NWK,IERROR
  integer K,I,J,IL,JL,NL,LL,Ydir

  logical, save :: firstcall = .true.

!!$  logical, parameter :: USE_CYCLIC = .false.
  !=========================================================================================

  call Timers_start("PFFT")

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
     IF(solveflag==1) THEN  
        LP =0 ! Periodic in X
        NP =1 ! Neumann in Y 
        MP =1 ! Neumann in Z
     ELSEIF(solveflag==2) THEN
        LP =0 ! Periodic in X
        NP =0 ! Periodic in Y 
        MP =1 ! Neumann in Z
     ELSEIF(solveflag==3) THEN
        LP =0 ! Periodic in X
        NP =0 ! Periodic in Y 
        MP =0 ! Preiodic in Z
     ENDIF
     
     ALLOCATE(AK(L))
     ALLOCATE(AM(M),BM(M),CM(M))
!-----------------------------------------------------------------------
!     MODIFIED WAVE NUMBERS
!-----------------------------------------------------------------------

     ! In X direction : pfft_wave(:,IAXIS)       
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

        !Please note: Blktri makes extensive use of assumed size arrays.
        !Compilers (absoft) are not able to keep track of allocated memory
        !bounds.  Please compile without array bounds checking!
        CALL BLKTRI(0,NP,N,AN,BN,CN,MP,M,AM,BM,CM,M,RHSP,IERROR,W)
        IF(IERROR/=0) WRITE(6,*) 'INIT. BLKTRI: IERROR = ',IERROR

        ELSE !IF(solveflag==2 .or. solveflag==3) THEN

          ALLOCATE(BMM(M),CMM(M),DPZ(M)) !,DPY(N)
          ALLOCATE(AL(N),AZ(M))
!-----------------------------------------------------------------------
!     MODIFIED WAVE NUMBERS
!-----------------------------------------------------------------------

          ! In y direction, pfft_wave(:,JAXIS)
          DO K=1,N/2
            AL(K)= 2.*PI*(K-1)
          ENDDO
          DO K=N/2+1,N
            AL(K)=-2.*PI*(N-K+1)
          ENDDO          
          AL(1:N)=2.*(1.-COS(AL(1:N)/REAL(N)))/DELYSQ


          ! 3 fft solution using Greens Function:
          ! In z direction, pfft_wave(:,KAXIS)
          DO K=1,M/2
            AZ(K)= 2.*PI*(K-1)
          ENDDO
          DO K=M/2+1,M
            AZ(K)=-2.*PI*(M-K+1)
          ENDDO
          AZ(1:M)=2.*(1.-COS(AZ(1:M)/REAL(M)))/DELZSQ          

          ALLOCATE(WSAVF(2*M+15))  
          CALL RFFTI(M,WSAVF)          

     ENDIF

  invScale = 1.0
  firstcall = .false.

  endif



  ! Forward Sweep and Solution:
  IF (iDirection .eq. PFFT_FORWARD) THEN

     ! A this point I have the data distributed along pencils in the x direction
     ! Forward transform of only in x (npx = 1, npy=np, npz=1). pfft_ndim==1,  direction==PFFT_FORWARD 
     numVec = pfft_inLen(JAXIS)*pfft_inLen(KAXIS)

     call gr_pfftDcftForward(inArray,pfft_work1,pfft_trigIaxis,&
                          pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                          numVec,pfft_transformType(IAXIS),pfft_scale(IAXIS)) ! Consecutive RFFTF

     call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work1,&                    
                           pfft_work2,pfft_inLen,pfft_midLen,&
                           pfft_procGrid(JAXIS),pfft_comm(JAXIS))


     if (solveflag .eq. 1) then   ! isplane 

        ! Check pz == 1
        ! If pz = 1  (px is also 1 always) and solver BLKTRI => Use Block Tridiagonal
        ! -- ------  --------------------  --- --------------------------------------
        !It is easier to work with 3D arrays.  Grid is orientated as YZX.
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

           J= pfftBlkLimits(LOW,KAXIS) + JL         ! GLOBAL INDEX   
           if (J-1 > pfftBlkLimits(HIGH,KAXIS)) cycle
  
           !-----------------------------------------------------------------------
           !    CHANGE THE DIAGONAL COEFFICIENTS
           !-----------------------------------------------------------------------
           BML(1:M) = BM(1:M) - AK(J/2+1) 

           !    ASSIGN THE RHS
           DO K=1,N ! Y Direction
              DO I=1,M ! Z Direction
                 RHSP(I,K) = tmp3DArray(K,I,JL)   
              ENDDO
           ENDDO
       
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
              ENDDO
           ENDDO

        ENDDO

        ! Need to pass the data from 3D array to 1D array pfft_work2, 
        pfft_work2(1:product(pfft_midLen)) = &
             reshape(tmp3DArray, (/product(pfft_midLen)/))
        deallocate(tmp3DArray)


        ! Transpose to go from yzx to xyz (we end up with npx=1, npy=np, npz=1)
        call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,pfft_work2,&                    
                              pfft_work1,pfft_midLen,pfft_inLen,&
                              pfft_procGrid(JAXIS),pfft_comm(JAXIS))
   

     else 


        ! Do Fourier Transform in Y direction:
        numVec=pfft_midLen(JAXIS)*pfft_midLen(KAXIS)
        call gr_pfftDcftForward(pfft_work2,pfft_work1,pfft_trigJaxis,  &                  
                              pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
                              numVec,pfft_transformType(JAXIS),pfft_scale(JAXIS)) ! Consecutive RFFTF

#if NDIM == 3  
        call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work1,  &
                            pfft_work2,pfft_midLen,pfft_outLen,&
                            pfft_procGrid(KAXIS),pfft_comm(KAXIS))        

        ! Need to pass the data to 3D arrays. Go from 1D pfft_work2. 
        ! Read to 2 periodic zxy sequence or yxz sequence (2D) 
        ! Grid is orientated as ZXY.
        allocate(tmp3DArray &
             (pfft_outLen(IAXIS), pfft_outLen(JAXIS), pfft_outLen(KAXIS)))
        tmp3DArray = reshape(pfft_work2, pfft_outLen)

      
        ! Orientation of the grid in tmp3DArray Z,X,Y
        gridOrientation(IAXIS) = KAXIS
        gridOrientation(JAXIS) = IAXIS     
        gridOrientation(KAXIS) = JAXIS

        ! Processor distribution on the Z,X,Y domain (1 proc,py procs,pz procs)
        pfft_meorientation(IAXIS) = IAXIS ! Refers to 1 proc
        pfft_meorientation(JAXIS) = JAXIS ! Refers to py
        pfft_meorientation(KAXIS) = KAXIS ! Refers to pz
        do i = 1, NDIM
           call gr_pfftGetLocalLimitsAnytime(i,gridOrientation(i),pfft_meorientation(i),pfft_outLen, &
                                               PFFT_PCLDATA_REAL,pfftBlkLimits)
        end do

        LL = pfft_outLen(JAXIS)  ! Number of points in X as tmp3DArray(z, x, y)
        NL = pfft_outLen(KAXIS)  ! Number of points in Y as tmp3DArray(z, x, y)

        Ydir = KAXIS

#else

        ! Need to pass the data to 3D arrays. Go from 1D pfft_work1. 
        ! Read to 2 periodic yxz sequence (2D) 
        allocate(tmp3DArray &
             (pfft_midLen(IAXIS), pfft_midLen(JAXIS), pfft_midLen(KAXIS)))
        tmp3DArray = reshape(pfft_work1, pfft_midLen)


        ! Orientation of the grid in tmp3DArray Y, X, Z
        gridOrientation(IAXIS) = JAXIS
        gridOrientation(JAXIS) = IAXIS     
        gridOrientation(KAXIS) = KAXIS

        ! Processor distribution on the Y,X,Z domain (1 proc,py procs,1)
        pfft_meorientation(IAXIS) = IAXIS
        pfft_meorientation(JAXIS) = JAXIS
        pfft_meorientation(KAXIS) = KAXIS
        do i = 1, NDIM
           call gr_pfftGetLocalLimitsAnytime(i,gridOrientation(i),pfft_meorientation(i),pfft_midLen, &
                                               PFFT_PCLDATA_REAL,pfftBlkLimits)
        end do

        LL = pfft_midLen(JAXIS)  ! Number of points in X as tmp3DArray(y, x, z)
        NL = pfft_midLen(IAXIS)  ! Number of points in Y as tmp3DArray(y, x, z)

        Ydir = IAXIS
#endif

        ! Operate on the Z direction
        DO JL=1,NL  ! Y Direction

           J= pfftBlkLimits(LOW,Ydir) + JL         ! GLOBAL INDEX Y DIRECTION   
!!$           if (J .LE. pfftBlkLimits(HIGH,Ydir)) then
!!$              print*,pfft_myPE,' J=',J,pfftBlkLimits(:,Ydir),' ok'
!!$           else
!!$              print*,pfft_myPE,' J=',J,pfftBlkLimits(:,Ydir),' CYCLE J!'
!!$           end if
           if (J-1 > pfftBlkLimits(HIGH,Ydir)) cycle
           
           !-----------------------------------------------------------------------
           !    CHANGE THE DIAGONAL COEFFICIENTS
           !-----------------------------------------------------------------------
           BML(1:M) = BM(1:M) - AL(J/2+1) 
           
           DO IL=1,LL

              I= pfftBlkLimits(LOW,JAXIS) + IL      ! GLOBAL INDEX X DIRECTION  
              if (I-1 > pfftBlkLimits(HIGH,JAXIS)) cycle

              !-----------------------------------------------------------------------
              !    CHANGE THE DIAGONAL COEFFICIENTS
              !-----------------------------------------------------------------------
              BMM(1:M) = BML(1:M) - AK(I/2+1) 

#if NDIM == 3
              DPZ(1:M) = tmp3DArray(1:M,IL,JL)

              IF (solveflag .eq. 2) THEN
                 CMM = CM
                 IF( J .eq. 1 .AND. I .eq. 1 ) THEN 
                    CMM(1)=0.
                    DPZ(1)=0.
                 ENDIF
                 call TRIDAG(AM,BMM,CMM,DPZ,DPZ,M)

              ELSE  ! solveflag .eq. 3              
                 
                 ! 3D fft solution
                 ! Fourier transform in z direction:
                 CALL RFFTF(M,DPZ,WSAVF)                 

                 ! Green's Function
                 K = 1
                 IF (I+J+K .eq. 3) THEN
                    DPZ(K) = 0.
                 ELSE
                    DPZ(K) = -DPZ(K)/(AL(J/2+1)+AK(I/2+1)+AZ(K/2+1))
                 ENDIF
                 DO K = 2,M
                    DPZ(K) = -DPZ(K)/(AL(J/2+1)+AK(I/2+1)+AZ(K/2+1))
                 ENDDO

                 ! Inverse Fourier Transform
                 DPZ = DPZ/real(M)
                 CALL RFFTB(M,DPZ,WSAVF)

              ENDIF

              tmp3DArray(1:M,IL,JL) = DPZ(1:M)
             
#else

              DPZ(1:M) = tmp3DArray(JL,IL,1:M)
              
              IF (J+I .eq. 2) then
                 DPZ(1) = 0.                           ! Wave numbers in x and y = 0.
              ELSE
                 DPZ(1) = DPZ(1)/(BMM(1)+AM(1)+CM(1))  ! 2D M==1, Rest of wavenumber combinations
              ENDIF
              tmp3DArray(JL,IL,1:M) = DPZ(1:M)
#endif

           ENDDO
        ENDDO

#if NDIM ==3
        ! Need to pass the data from 3D array to 1D array pfft_work2, 
        pfft_work2(1:product(pfft_outLen)) = &
             reshape(tmp3DArray, (/product(pfft_outLen)/))
        deallocate(tmp3DArray)


        ! Transpose to go from zxy to yzx (we end up with npx=py, npy=1, npz=pz)
        call gr_pfftTranspose(PFFT_INVERSE,PFFT_PCLDATA_REAL,pfft_work2,  &                    
                              pfft_work1,pfft_outLen,pfft_midLen,&
                              pfft_procGrid(KAXIS),pfft_comm(KAXIS))
#else

        ! Need to pass the data from 3D array to 1D array pfft_work1, 
        pfft_work1(1:product(pfft_midLen)) = &
             reshape(tmp3DArray, (/product(pfft_midLen)/))
        deallocate(tmp3DArray)

#endif
        ! Solution still needs inverse fft in Y, Transpose to xyz and ifft in X
      
     endif

  ELSE     ! DO INVERSE

     if (solveflag .ne. 1) then ! not isplane
        ! Do other inverse transforms (related to solver flags 2 or 3)
        numVec=pfft_midLen(JAXIS)*pfft_midLen(KAXIS)
        call gr_pfftDcftInverse(pfft_work1,pfft_work2,pfft_trigJaxis,  &
                              pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
                              numVec,pfft_transformType(JAXIS),invScale)


#if NDIM==3
        !-----------------------------------------------------------------------
        !    FIX THE SUM OF PRESSURE CORRECTION OF WHOLE FIELD TO ZERO
        !-----------------------------------------------------------------------
        ! Need to pass the data to 3D arrays. Go from 1D pfft_work2.
        ! Grid orientation is YZX.
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

        IF(pfftBlkLimits(LOW,KAXIS)==0) THEN  ! Wavenumber 0 is in this block          
           DPM = SUM(tmp3DArray(:,:,1))/REAL(M*N) 
        ELSE
           DPM = 0.
        ENDIF

        DPMaux = DPM
        call MPI_Allreduce(DPMaux, DPM, 1, FLASH_REAL,&
                           MPI_SUM, pfft_comm(IAXIS), ierr)

        IF(pfftBlkLimits(LOW,KAXIS)==0) THEN  ! Wavenumber 0 is in this block          
           tmp3DArray(:,:,1) = tmp3DArray(:,:,1) - DPM
        ENDIF

        ! Need to pass the data from 3D array to 1D arry pfft_work2, 
        pfft_work2(1:product(pfft_midLen)) = &
             reshape(tmp3DArray, (/product(pfft_midLen)/))
        deallocate(tmp3DArray)
#endif

        !-----------------------------------------------------------------------
        call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work2,   &
                              pfft_work1,pfft_midLen,pfft_inLen,&
                              pfft_procGrid(JAXIS),pfft_comm(JAXIS))


     endif

     ! Inverse transform of only in x. pfft_ndim==1,  direction==PFFT_INVERSE 
     numVec=pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     call gr_pfftDcftInverse(pfft_work1,outArray,pfft_trigIaxis,     &
                             pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                             numVec,pfft_transformType(IAXIS),invScale)


  ENDIF

  call Timers_stop("PFFT")

  return
end subroutine gr_pfftPoissonDirect


      SUBROUTINE tridag(a,b,c,r,u,n)
      use Driver_interface, ONLY : Driver_abortFlash 
      INTEGER n,NMAX
      REAL a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=3000)
      INTEGER j
      REAL bet,gam(NMAX)
      if(b(1).eq.0.) call Driver_abortFlash('tridag: rewrite equations')
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)call Driver_abortFlash('tridag failed')
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
