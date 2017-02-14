!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/rt_init
!!
!!  NAME 
!!
!!  rt_init
!!
!!  SYNOPSIS
!!
!!  call rt_init()
!!
!!  DESCRIPTION 
!!    Initialize neutrino leakage
!!
!!  NOTES
!!      This unit implements ray-by-ray multispecies neutrino leakage.
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in rt_calcLeak.F90 and 
!!      rt_calcTau.F90 are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  stellarcollapse.org/codes.html.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M., & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine rt_init
  use rt_data
  use RadTrans_data, ONLY : rt_useRadTrans
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm, Driver_abortFlash, &
                               Driver_getNStep
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  Use Logfile_interface, ONLY : Logfile_stamp

#include "Flash.h"
#include "constants.h"
  
  implicit none

  include "Flash_mpi.h"
  
  integer :: i,j,k
  real :: dTht, dPhi

  integer :: nconstant, nlog, n1, ghosts1
  real :: dx, xmin, dxfactor
  real, allocatable :: x1(:), x1i(:)
  integer :: error, raysPerProc, count
  logical :: threadBlockListBuild, threadWithinBlockBuild

  integer :: numRaysPerProc, istart, iend, lnum, num

  call Driver_getMype(GLOBAL_COMM,rt_globalMe)
  call Driver_getComm(MESH_COMM, rt_meshComm)
  call Driver_getNStep(rt_nstepStart)

  call RuntimeParameters_get('leak_radLog', rt_leakRadLog)
  call RuntimeParameters_get('leak_dx', rt_leakDx)
  call RuntimeParameters_get('leak_radMax', rt_leakRadMax)
  call RuntimeParameters_get('leak_thtMax', rt_leakThtMax) ! in units of Pi
  call RuntimeParameters_get('leak_phiMax', rt_leakPhiMax) ! in units of Pi
  call RuntimeParameters_get('leak_doHeat', rt_leakDoHeat)
  call RuntimeParameters_get('leak_heatFac', rt_leakHeatFac)

#ifndef LEAK_STATIC
  !! These are initialized in rt_data.F90, if LEAK_STATIC
  call RuntimeParameters_get('leak_numRad', rt_leakNumRad)
  call RuntimeParameters_get('leak_numTht', rt_leakNumTht)
  call RuntimeParameters_get('leak_numPhi', rt_leakNumPhi)
#endif

  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadLeakBlockList", rt_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadLeakWithinBlock", rt_threadWithinBlock)
  call RuntimeParameters_get("leak_subCommSize", rt_subCommSize)

  call RuntimeParameters_get("leak_reducedTime", rt_reducedTime)
  call RuntimeParameters_get("leak_reducedSteps", rt_reducedSteps)

  if (rt_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Leak_init]')
     rt_threadBlockList = .false.
  end if
  if (rt_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Leak_init]')
     rt_threadWithinBlock = .false.
  end if

#ifndef LEAK_STATIC
  if (NDIM == 1 .AND. (rt_leakNumTht > 1 .OR. rt_leakNumPhi > 1)) then
     rt_leakNumTht = 1
     rt_leakNumPhi = 1
     if (rt_globalMe == MASTER_PE) write(*,*) 'Setting rt_leakNumTht/Phi=1 because NDIM=1'
  else if (NDIM == 2 .AND. rt_leakNumPhi > 1) then
     rt_leakNumPhi = 1
     if (rt_globalMe == MASTER_PE) write(*,*) 'Setting rt_leakNumPhi=1 because NDIM=2'
  end if
#endif 

  if (.NOT. rt_useRadTrans) return

  ghosts1 = 0
  n1 = rt_leakNumRad + 2*ghosts1
#ifndef LEAK_STATIC
  rt_leakNumRays = rt_leakNumTht*rt_leakNumPhi
#endif
  rt_arraySize = 5*3*n1*rt_leakNumRays
  rt_arraySmall = 3*n1*rt_leakNumRays

#ifndef LEAK_STATIC
  allocate(rt_leakRadii(n1))
  allocate(rt_dr(n1))
  allocate(rt_leakTheta(rt_leakNumTht))
  allocate(rt_leakPhi(rt_leakNumPhi))
  allocate(rt_leakX(n1,rt_leakNumTht,rt_leakNumPhi))
  allocate(rt_leakY(n1,rt_leakNumTht,rt_leakNumPhi))
  allocate(rt_leakZ(n1,rt_leakNumTht,rt_leakNumPhi))
#endif

  ! Indexing of rt_rayData(1:5,...):
  ! 1 - tauRuff
  ! 2 - chiRoss
  ! 3 - heatFlux
  ! 4 - heatErms
  ! 5 - heatEave
#ifndef LEAK_STATIC
  allocate(rt_rayData(5,n1,3,rt_leakNumRays))
  rt_rayData = 0.
#endif

  ! Set up pointers
  rt_tauRuff => rt_rayData(1,:,:,:)
  rt_chiRoss => rt_rayData(2,:,:,:)
  rt_heatFlx => rt_rayData(3,:,:,:)
  rt_heatErms => rt_rayData(4,1,:,:)
  rt_heatEave => rt_rayData(5,1,:,:)

  ! Allocate variables in leak_module
  allocate(x1(n1))
  allocate(x1i(n1))

  ! Now setup the radial spacing of the leakage rays
  ! have constant spacing up to leakRadLog, then logarithmic spacing
  ! put equal number of zones in each region, unless leakRadLog = 0,
  ! in which case the whole ray will be logarithmic, or 
  ! leakRadLog = leakRadMax, in which case the whole ray will be 
  ! logarithmically-spaced.
  x1i = 0.
  x1 = 0.
  nconstant = nint(rt_leakRadLog/rt_leakDx)
  if (nconstant.eq.0) then
     stop "grid_custom parameters wrong, or use log"
  endif
  nlog=n1-ghosts1*2-nconstant
  xmin = real(nconstant)*rt_leakDx
  call series2(nlog,xmin,rt_leakRadMax,rt_leakDx,dxfactor)
  x1(ghosts1+1) = rt_leakDx/2.0d0
  x1i(ghosts1+2) = rt_leakDx

  do i=ghosts1+3,ghosts1+nconstant
     dx = rt_leakDx
     x1i(i) = x1i(i-1) + dx
  enddo
     
  do i=ghosts1+nconstant+1,n1
     dx = dx*dxfactor
     x1i(i) = x1i(i-1)+dx
  enddo
     
  do i=ghosts1+2,n1-1
     x1(i) = 0.5d0*(x1i(i+1)+x1i(i))
  enddo
     
  x1(n1) = 2.0d0*x1i(n1) - x1(n1-1)
  
  ! setup inner ghost zones !! Hold-over!  there are no ghost zones here.
  j=ghosts1+1
  do i=ghosts1,1,-1
     x1(i) = -x1(j)
     x1i(i) = -x1i(j+1)
     j=j+1
  enddo
 
  rt_leakRadii = x1
  do i=1,rt_leakNumRad-1
     rt_dr(i) = x1i(i+1) - x1i(i)
  end do
  rt_dr(rt_leakNumRad) = rt_dr(rt_leakNumRad-1)

  dTht = 0.
  rt_leakTheta = PI/2.
  if (rt_leakNumTht > 1) then 
     dTht = PI*rt_leakThtMax/(rt_leakNumTht-1)
     do i=1,rt_leakNumTht
        rt_leakTheta(i) = (i-1)*dTht
     end do
  end if
  rt_dTht = dTht

  dPhi = 0.
  rt_leakPhi = 0.
  if (rt_leakNumPhi > 1) then
     dPhi = PI*rt_leakPhiMax/(rt_leakNumPhi-1)
     do i=1,rt_leakNumPhi
        rt_leakPhi(i) = (i-1)*dPhi
     end do
  end if
  rt_dPhi = dPhi

  ! Initial arrays
  rt_leakX = 0.0
  rt_leakY = 0.0
  rt_leakZ = 0.0
  do k=1,rt_leakNumPhi
     do j=1,rt_leakNumTht
        do i=1,n1
           rt_leakX(i,j,k) = rt_leakRadii(i)*sin(rt_leakTheta(j))*cos(rt_leakPhi(k))
#if NDIM > 1
           rt_leakY(i,j,k) = rt_leakRadii(i)*cos(rt_leakTheta(j))
#if NDIM == 3
           rt_leakZ(i,j,k) = rt_leakRadii(i)*sin(rt_leakTheta(j))*sin(rt_leakPhi(k))
#endif
#endif
        end do
     end do
  end do

  !initialize PNS coords
  rt_pnsCoord = 0.

  ! Now setup the arrays for leakage data
#ifndef LEAK_STATIC
  allocate(rt_leakArr(3,n1,rt_leakNumRays))
  allocate(rt_leakSrc(4,n1,rt_leakNumRays))
#endif
  rt_leakArr = 0.
  rt_leakSrc = 0.

  call MPI_COMM_SIZE(rt_meshComm, rt_meshNumProcs,error)

  if (mod(rt_meshNumProcs,rt_subCommSize) /= 0 .OR. rt_subCommSize == -1) then
     if (rt_globalMe == MASTER_PE) write(*,*) "RadTrans: setting leak_subCommSize to meshNumProcs"
     rt_subCommSize = rt_meshNumProcs
  end if

  rt_subMeshMe =  mod(rt_globalMe, rt_subCommSize)
  call MPI_COMM_SPLIT(rt_meshComm, rt_globalMe/rt_subCommSize, rt_subMeshMe, rt_subMeshComm, error)

  ! Now let's parse up the rays among different processors.
  allocate(rt_recvCnt(rt_subCommSize))
  allocate(rt_dsplCnt(rt_subCommSize))
  rt_recvCnt = 0
  rt_dsplCnt = 0

  numRaysPerProc = rt_leakNumRays / rt_subCommSize
  num = mod(rt_leakNumRays, rt_subCommSize)
  iend = 0
  do i=1,rt_subCommSize
     if (i <= num) then
        lnum = numRaysPerProc+1
     else if (i <= rt_leakNumRays) then
        lnum = numRaysPerProc
     else
        lnum = 0
     end if
     istart = iend + 1
     iend = istart + lnum - 1
     rt_recvCnt(i) = lnum * 3*5*n1
     rt_dsplCnt(i) = (istart-1) *3*5*n1
     if (rt_subMeshMe == i-1) then
        rt_istart = istart
        rt_iend = iend
     end if
  end do

  rt_arraySizeLocal = 3*n1*(rt_iend-rt_istart+1)


end subroutine rt_init

subroutine series2(nzones,xmin,xmax,mindx,dxfac)
! This routine is motivated by the "series2" subroutine of
! Cala resp. Prometheus.
!
! It solves for a factor dxfac by which each dx is slightly
! larger than the preceding dx.

  implicit none
  
  real dxfac
  integer nzones
  real xmin,xmax,mindx
  
  ! internal vars
  real tol
  real al,aold,ferror,sum,F,dsum,dFda,anew
  integer k,i,itermax
  
  tol = 1.0d-6
  itermax = 100
  
  ! solve for dxfac
  dxfac=0.0d0

  ! estimate
  al = log( (xmax-xmin)/mindx )/ ( real(nzones-2) )
  aold = exp(al)
  k = 1
  ferror = 1.0d0

  !-------------------------------------------------
  ! Solve: F = (xmax-xmin)/mindx - (Sum[ a^j],j=0,N)
  ! let x = a, y(x) = F
  !-------------------------------------------------

  ! evaluate F
  do while( (ferror.gt.tol).and.(k.lt.itermax))
     sum = 0.0d0
     do i=1,nzones-1
        sum = sum + aold**(i-1)
     enddo
     F = ( (xmax-xmin)/mindx) - sum
     
     ! evaluate dFDa
     dsum = 1.0d0
     do i=4,nzones
        dsum = dsum + (i-2)*(aold**(i-3))
     enddo
     dFda = -1.0d0*dsum
     ! next root
     anew = aold - F/dFda
     ferror = abs(anew-aold)/aold
     k = k + 1
     aold = anew
  enddo
  dxfac = anew
  
end subroutine  series2
     
