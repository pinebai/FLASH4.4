!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/Generic_Direct/gr_pfftPoissonDirect
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


subroutine gr_pfftPoissonDirect (iDirection, &
              solveflag, inSize, localSize, globalSize, &
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
  
  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe,gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
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
  
  TYPE (fishworkspace), SAVE :: W
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: W2

  REAL,DIMENSION(:),ALLOCATABLE, SAVE :: DPZ

  real :: DPM,DPMaux
!
!     LOCAL ARRAYS FOR FFT
!
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: AK,AL,AZ

  real, dimension(:,:,:), allocatable :: tmp3DArray

  integer, save :: LP,NP,MP,N,M,L,NWK,IERROR
  integer K,I,J,IL,JL,NL,LL,Ydir,n2mh,ib,jb,kb

  logical, save :: firstcall = .true.

  !=========================================================================================

  call Timers_start("PFFT")

  ! Initialization of solver arrays.
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
    
!-----------------------------------------------------------------------
!     DIMENSION OF ARRAY FOR WORK SPACE AND SAVED ARRAYS
!-----------------------------------------------------------------------

     write(*,*) 'transformtype(IAXIS)',transformtype(IAXIS)
     write(*,*) 'transformtype(JAXIS)',transformtype(JAXIS)
     write(*,*) 'transformtype(KAXIS)',transformtype(KAXIS)
     write(*,*) 'L,N,M=',L,N,M
    
     ALLOCATE(AK(L),AL(N),AZ(M),DPZ(M))

!-----------------------------------------------------------------------
!     MODIFIED WAVE NUMBERS
!-----------------------------------------------------------------------

     ! In X direction : pfft_wave(:,IAXIS)   
     ! X wavenumbers:
     select case(transformtype(IAXIS))
     case(PFFT_SINQ) ! Dirichlet BCs
         do k=1 , L 
            AK(k) = cos(PI*real(k)/real(L)); 
            AK(k) = 2.*(1.-AK(k));
         enddo

    case(PFFT_COSQ) ! Neuman BCs  
         AK(1) = 0.0            ! initialize - KW
         do k=1 , L-1  
            AK(k+1) = cos(PI*real(k)/real(L)); 
            AK(k+1) = 2.*(1.-AK(k+1));
         enddo

    case(PFFT_REAL)
         n2mh = L/2;
         do k=1,n2mh
             AK(k) = 2.*PI*real(k-1)
         enddo
         do k=n2mh+1,L
             AK(k) = -2.*PI*real(L+1-k)
         enddo
         AK(1:L)=2.*(1.-cos(AK(1:L)/real(L)))

    end select             

    AK(1:L) = AK(1:L)/DELXSQ

    
     ! In Y direction : pfft_wave(:,JAXIS)   
     ! Y wavenumbers:
     select case(transformtype(JAXIS))
     case(PFFT_SINQ) ! Dirichlet BCs
         do k=1 , N 
            AL(k) = cos(PI*real(k)/real(N)); 
            AL(k) = 2.*(1.-AL(k));
         enddo

    case(PFFT_COSQ) ! Neuman BCs  
         AL(1) = 0.0            ! initialize - KW
         do k=1 , N-1  
            AL(k+1) = cos(PI*real(k)/real(N)); 
            AL(k+1) = 2.*(1.-AL(k+1));
         enddo

    case(PFFT_REAL)
         n2mh = N/2;
         do k=1,n2mh
             AL(k) = 2.*PI*real(k-1)
         enddo
         do k=n2mh+1,N
             AL(k) = -2.*PI*real(N+1-k)
         enddo
         AL(1:N)=2.*(1.-cos(AL(1:N)/real(N)))

    end select             

    AL(1:N) = AL(1:N)/DELYSQ


#if NDIM == 3
     ! In Z direction : pfft_wave(:,KAXIS)   
     ! Z wavenumbers:
     select case(transformtype(KAXIS))
     case(PFFT_SINQ) ! Dirichlet BCs
         do k=1 , M 
            AZ(k) = cos(PI*real(k)/real(M)); 
            AZ(k) = 2.*(1.-AZ(k));
         enddo

    case(PFFT_COSQ) ! Neuman BCs  
         AZ(1) = 0.0            ! initialize - KW
         do k=1 , M-1  
            AZ(k+1) = cos(PI*real(k)/real(M)); 
            AZ(k+1) = 2.*(1.-AZ(k+1));
         enddo

    case(PFFT_REAL)
         n2mh = M/2;
         do k=1,n2mh
             AZ(k) = 2.*PI*real(k-1)
         enddo
         do k=n2mh+1,M
             AZ(k) = -2.*PI*real(M+1-k)
         enddo
         AZ(1:M)=2.*(1.-cos(AZ(1:M)/real(M)))

    end select             

    AZ(1:M) = AZ(1:M)/DELZSQ 
#endif

    invScale = 1.0
    firstcall = .false.

  endif



  ! Forward Sweep and Solution:
  IF (iDirection .eq. PFFT_FORWARD) THEN


     ! A this point I have the data distributed along pencils in the x direction
     ! Forward transform of only in x (npx = 1, npy, npz). pfft_ndim==1,  direction==PFFT_FORWARD 
     numVec = pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     call gr_pfftDcftForward(inArray,pfft_work1,pfft_trigIaxis,&
                          pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                          numVec,pfft_transformType(IAXIS),pfft_scale(IAXIS)) 

     call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work1,&                    
                           pfft_work2,pfft_inLen,pfft_midLen,&
                           pfft_procGrid(JAXIS),pfft_comm(JAXIS))

     ! Do Transform Transform in Y direction:
     numVec=pfft_midLen(JAXIS)*pfft_midLen(KAXIS)
     call gr_pfftDcftForward(pfft_work2,pfft_work1,pfft_trigJaxis,  &                  
                             pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
                             numVec,pfft_transformType(JAXIS),pfft_scale(JAXIS)) 

#if NDIM == 3  
        call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work1,  &
                            pfft_work2,pfft_midLen,pfft_outLen,&
                            pfft_procGrid(KAXIS),pfft_comm(KAXIS))


        ! Do Transform in Z direction:
        numVec=pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
        call gr_pfftDcftForward(pfft_work2,pfft_work1,pfft_trigKaxis,  &                  
                             pfft_globalLen(KAXIS),pfft_outLen(IAXIS),&
                             numVec,pfft_transformType(KAXIS),pfft_scale(KAXIS)) 


        ! Need to pass the data to 3D arrays. Go from 1D pfft_work1. 
        ! Read to zxy sequence (3D) 
        ! Grid is orientated as ZXY.
        allocate(tmp3DArray &
             (pfft_outLen(IAXIS), pfft_outLen(JAXIS), pfft_outLen(KAXIS)))
        tmp3DArray = reshape(pfft_work1, pfft_outLen)

      
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
           if (J-1 > pfftBlkLimits(HIGH,Ydir)) cycle

           ! Define wavenumber index:
           select case(transformtype(JAXIS))
           case(PFFT_SINQ)
              jb = J
           case(PFFT_COSQ)
              jb = J
           case(PFFT_REAL)
              jb = (J)/2+1
           end select


           DO IL=1,LL

              I= pfftBlkLimits(LOW,JAXIS) + IL      ! GLOBAL INDEX X DIRECTION  
              if (I-1 > pfftBlkLimits(HIGH,JAXIS)) cycle

              ! Define wavenumber index:
              select case(transformtype(IAXIS))
              case(PFFT_SINQ)
                 ib = I
              case(PFFT_COSQ)
                 ib = I
              case(PFFT_REAL)
                 ib = (I)/2+1
              end select

#if NDIM == 3

              DPZ(1:M) = tmp3DArray(1:M,IL,JL)          
                 
              ! Green's Function Solution
              K = 1
              ! Define wavenumber index:
              select case(transformtype(KAXIS))
              case(PFFT_SINQ)
                 kb = K
              case(PFFT_COSQ)
                 kb = K
              case(PFFT_REAL)
                 kb = (K)/2+1
              end select
              ! Set to zero constant coeff, if no Dirichlet BC:
              if((transformtype(IAXIS) .ne. PFFT_SINQ) .and. &
                 (transformtype(JAXIS) .ne. PFFT_SINQ) .and. &
                 (transformtype(KAXIS) .ne. PFFT_SINQ) .and. &
                 (I+J+K .eq. 3)) then
                 DPZ(K) = 0.
              else
                 DPZ(K) = -DPZ(K)/(AL(jb)+AK(ib)+AZ(kb))
              endif

              DO K = 2,M
                 ! Define wavenumber index:
                 select case(transformtype(KAXIS))
                 case(PFFT_SINQ)
                    kb = K
                 case(PFFT_COSQ)
                    kb = K
                 case(PFFT_REAL)
                    kb = (K)/2+1
                 end select
                 DPZ(K) = -DPZ(K)/(AL(jb)+AK(ib)+AZ(kb))
              ENDDO

              tmp3DArray(1:M,IL,JL) = DPZ(1:M)
             
#else

              DPZ(1:M) = tmp3DArray(JL,IL,1:M)
              
              K = 1
              ! Set to zero constant coeff, if no Dirichlet BC:
              if((transformtype(IAXIS) .ne. PFFT_SINQ) .and. &
                 (transformtype(JAXIS) .ne. PFFT_SINQ) .and. &
                 (J+I .eq. 2)) then
                 DPZ(K) = 0. 
              else
                 DPZ(K) = -DPZ(K)/(AL(jb)+AK(ib))  ! 2D M==1, Rest of wavenumber combinations
              endif

              tmp3DArray(JL,IL,1:M) = DPZ(1:M)
#endif

           ENDDO
        ENDDO

#if NDIM ==3
        ! Need to pass the data from 3D array to 1D array pfft_work1, 
        pfft_work1(1:product(pfft_outLen)) = &
             reshape(tmp3DArray, (/product(pfft_outLen)/))
        deallocate(tmp3DArray)

#else

        ! Need to pass the data from 3D array to 1D array pfft_work1, 
        pfft_work1(1:product(pfft_midLen)) = &
             reshape(tmp3DArray, (/product(pfft_midLen)/))
        deallocate(tmp3DArray)

#endif
        ! Solution still needs inverse fft in Y, Transpose to xyz and ifft in X
      

  ELSE     ! DO INVERSE

#if NDIM ==3

     ! Inverse transforms in Z direction:
     numVec=pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
     call gr_pfftDcftInverse(pfft_work1,pfft_work2,pfft_trigKaxis,  &
                             pfft_globalLen(KAXIS),pfft_outLen(IAXIS),&
                             numVec,pfft_transformType(KAXIS),invScale)

     ! Transpose to go from zxy to yzx (we end up with npx=py, npy=1, npz=pz)
     call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work2,  &                    
                           pfft_work1,pfft_outLen,pfft_midLen,&
                           pfft_procGrid(KAXIS),pfft_comm(KAXIS))
#endif


     ! Do other inverse transforms 
     ! Inverse transforms in Y direction:
     numVec=pfft_midLen(JAXIS)*pfft_midLen(KAXIS)
     call gr_pfftDcftInverse(pfft_work1,pfft_work2,pfft_trigJaxis,  &
                             pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
                             numVec,pfft_transformType(JAXIS),invScale)

     !-----------------------------------------------------------------------
     call gr_pfftTranspose(idirection,PFFT_PCLDATA_REAL,pfft_work2,   &
                           pfft_work1,pfft_midLen,pfft_inLen,&
                           pfft_procGrid(JAXIS),pfft_comm(JAXIS))


     ! Inverse transform of only in x. pfft_ndim==1,  direction==PFFT_INVERSE 
     numVec=pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
     call gr_pfftDcftInverse(pfft_work1,outArray,pfft_trigIaxis,     &
                             pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                             numVec,pfft_transformType(IAXIS),invScale)

  ENDIF

  call Timers_stop("PFFT")

  return
end subroutine gr_pfftPoissonDirect
