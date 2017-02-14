!!****if* source/Grid/GridSolvers/Multigrid/gr_hgInit
!!
!! NAME
!!  gr_hgInit
!!
!! SYNOPSIS
!!
!!  gr_hgInit()
!!
!! DESCRIPTION
!!  This routine initializes the data necessary for HG Multigrid
!!
!! PARAMETERS
!!
!!  igeomx, igeomy, igeomz, geometry, quadrant, nblockx, nblocky, nblockz
!!  mg_maxcorrections, mg_maxresidualnorm, mg_printnorm
!!
!!  NOTE
!!    Although there is code here for the non-Cartesian case,
!!    it is suspicious (e.g. loads of assumptions about where 0 is located)
!!    and has never been tested.  Therefore, only Cartesian geometry is allowed.
!!
!!***

#include "Flash.h"

subroutine gr_hgInit

  !==================================================================
  use gr_hgData, ONLY: hg_myPE, gr_hgMeshRefineMin, gr_hgMeshRefineMax, &
       gr_hgMaxCorrections, gr_hgPrintNorm, gr_hgMaxResidualNorm, &
       gr_hgQuadrant, gr_hgGeometry, gr_hgSaveNodetype, gr_hgSaveNewchild, &
       Px,Py,Pz,Pew,Pns,Pud,&
       hg_cx, hg_cy, hg_cz, &
       hg_ili,hg_ile,hg_iui,hg_iue,hg_jli,hg_jle,hg_jui,hg_jue,hg_kli,hg_kle,hg_kui,hg_kue, &
       n1off, n2off, n3off, nmax1, nmax2, &
       nxl1,nxl2,nxr1,nxr2,nyl1,nyl2,nyr1,nyr2,nzl1,nzl2,nzr1,nzr2, &
       hg_restrict_n1, hg_restrict_n2, hg_restrict_n3, &
       send_restrict_data, send_prolong_data, recv_restrict_data, recv_prolong_data, &
       send_prolong_req, send_restrict_req,  nbbuf_prolong, nbbuf_restrict, &
       gr_hgUsingPfft

#ifdef FLASH_GRID_PARAMESH2
  use tree, ONLY : nguard, maxblocks
#else
  use paramesh_dimensions, ONLY : nguard, maxblocks
#endif

  use Grid_data, ONLY : gr_meshComm, gr_meshMe, gr_meshNumProcs, &
       gr_nblockx, gr_nblocky, gr_nblockz
  use tree, ONLY : lrefine, nodetype, newchild, nchild, nfaces, lnblocks, maxblocks_tr
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY:  Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Multigrid.h"
#include "Flash_mpi.h"


  integer                         :: lb, ierr, ierr1, ierr2, ierr3, c
  integer                            :: mylrefmin, mylrefmax
  integer                         :: nbr_blks(2*MDIM)
  integer                         :: i,j,k

  integer                         :: igeomx, igeomy, igeomz

  character(len=MAX_STRING_LENGTH) :: strGeometry,eosModeString
  integer                         :: geometry, refineLevel, blkCount

  !=======================================================================

  allocate(gr_hgSaveNodetype(maxblocks_tr),stat=ierr1)
  allocate(gr_hgSaveNewchild(maxblocks_tr),stat=ierr2)
  if ((ierr1 /= 0) .or. (ierr2 /= 0)) &
       call Driver_abortFlash("Allocation error in hgInit for gr_hgSaveNodetype or gr_hgSaveNewchild")
  

  hg_myPE = gr_meshMe
  do lb = 1, lnblocks
     gr_hgSaveNodetype(lb) = nodetype(lb)
     gr_hgSaveNewchild(lb) = newchild(lb)
  end do

  ! Determine minimum and maximum levels of refinement in the mesh.

  mylrefmin = HUGE(0)
  mylrefmax = TINY(0.0)

  do lb = 1, lnblocks
     mylrefmin = min(mylrefmin,lrefine(lb))
     mylrefmax = max(mylrefmax,lrefine(lb))
  end do
  
  call mpi_allreduce(mylrefmin, gr_hgMeshRefineMin, 1, MPI_INTEGER, MPI_MIN, &
       gr_meshComm, ierr)
  call mpi_allreduce(mylrefmax, gr_hgMeshRefineMax, 1, MPI_INTEGER, MPI_MAX, &
       gr_meshComm, ierr)
  
  !  mylrefmin = minval(lrefine(1:lnblocks))
  !  call mpi_allreduce(mylrefmin, gr_hgMeshRefineMin, 1, MPI_INTEGER, MPI_MIN, &
  !       gr_meshComm, ierr)
  
  !  mylrefmax = maxval(lrefine(1:lnblocks))
  !  call mpi_allreduce(mylrefmax, gr_hgMeshRefineMax, 1, MPI_INTEGER, MPI_MAX, &
  !       gr_meshComm, ierr)
  
  call Logfile_stamp(gr_hgMeshRefineMax,"[gr_hgInit] max refine level = ")
  
  hg_ili = 1 + NGUARD
  hg_iui = NXB + NGUARD
  hg_jli = 1 + K2D*NGUARD
  hg_jui = NYB + K2D*NGUARD
  hg_kli = 1 + K3D*NGUARD
  hg_kui = NZB + K3D*NGUARD
  
  ! Determine index ranges for exterior zones.
  
  hg_ile = 1
  hg_iue = NXB + 2*NGUARD
  hg_jle = 1
  hg_jue = NYB + 2*NGUARD*K2D
  hg_kle = 1
  hg_kue = NZB + 2*NGUARD*K3D
  
  ! Determine mesh geometry and decide whether we support it (only on the first
  ! call).
  
  call RuntimeParameters_get ("geometry", strGeometry)
  call RuntimeParameters_mapStrToInt(strGeometry, geometry)
  call RuntimeParameters_get("quadrant", gr_hgQuadrant)


!! NOTE -- this routine used to work by asking about geometry in all three directions
!! However, since ONLY cartesian is supported, LBR eliminated the options
!! If they should be reimplemented, look at Hydro_init to see how to get them into 
!! a sensible structure.  Or possible Grid_initGeometry
!  call RuntimeParameters_get("igeomx", igeomx)
!  call RuntimeParameters_get("igeomy", igeomy)
!  call RuntimeParameters_get("igeomz", igeomz)

  select case (NDIM)
  case (1)
     if (geometry == CARTESIAN) then
        gr_hgGeometry = MG_GEOM_1DCARTESIAN
        !    elseif (igeomx == geom_cylrad) then
        !      gr_hgGeometry = MG_GEOM_1DCYLINDRICAL
        !    elseif (igeomx == geom_sphrad) then
        !      gr_hgGeometry = MG_GEOM_1DSPHERICAL
     else
        gr_hgGeometry = MG_GEOM_INVALID
     endif
  case (2)
     if (geometry == CARTESIAN) then
        gr_hgGeometry = MG_GEOM_2DCARTESIAN
        !    elseif ((igeomx == geom_cylrad) .and. &
        !            (igeomy == geom_planar)) then
        !      gr_hgGeometry = MG_GEOM_2DCYLAXISYM
        !    elseif ((igeomx == geom_cylrad) .and. &
        !            (igeomy == geom_cylang)) then
        !      gr_hgGeometry = MG_GEOM_2DCYLPOLAR
        !    elseif ((igeomx == geom_sphrad) .and. &
        !            (igeomy == geom_sphtheta)) then
        !      gr_hgGeometry = MG_GEOM_2DSPHAXISYM
     else
        gr_hgGeometry = MG_GEOM_INVALID
     endif
  case (3)
     if (geometry == CARTESIAN) then
        gr_hgGeometry = MG_GEOM_3DCARTESIAN
        !    elseif ((igeomx == geom_cylrad) .and. &
        !            (igeomy == geom_cylang) .and. &
        !            (igeomz == geom_planar)) then
        !      gr_hgGeometry = MG_GEOM_3DCYLINDRICAL
        !    elseif ((igeomx == geom_sphrad) .and. &
        !            (igeomy == geom_sphtheta) .and. &
        !            (igeomz == geom_sphphi)) then
        !      gr_hgGeometry = MG_GEOM_3DSPHERICAL
     else
        gr_hgGeometry = MG_GEOM_INVALID
     endif
  end select
  if (gr_hgGeometry == MG_GEOM_INVALID) then
     call Driver_abortFlash ('[gr_hgInit]  unsupported grid geometry, only Cartesian allowed')
  endif

  ! Make sure we only have one mesh block on the coarsest level.
  if ((gr_nblockx /= 1) .or. (gr_nblocky /= 1) .or. (gr_nblockz /= 1)) then
     if (.not.gr_hgUsingPfft) then
        call Driver_abortFlash("[gr_hgInit] only one block allowed on coarsest level")
     end if
  endif

  
  if (NDIM == 1) then
     nmax1 = 1
     nmax2 = 1
  elseif (NDIM == 2) then
     nmax1 = max(NXB, NYB)
     nmax2 = 1
  else ! NDIM == 3
     nmax1 = max(NXB, NYB)
     nmax2 = max(NYB, NZB)
  endif
  
  allocate(Px(-2:2,NXB,nchild), stat=ierr1)
  allocate(Py(-2:2,NYB,nchild), stat=ierr2)
  allocate(Pz(-2:2,NZB,nchild), stat=ierr3)
  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
       call Driver_abortFlash("Allocation error for the P<dim> arrays")

  ! It helps if you initialize variables before you allocate with them LBR  
  nbbuf_restrict = max(maxblocks, 1)
  nbbuf_prolong = max(maxblocks/8, 1)
  
  !initialize the prolongation data structures

  allocate(send_prolong_data(nmax1,nmax2,nfaces,nchild,nbbuf_prolong), &
           stat=ierr1)
  allocate(recv_prolong_data(nmax1,nmax2,nfaces), stat=ierr2)
  allocate(send_prolong_req(nchild*nbbuf_prolong), stat=ierr3)
  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
       call Driver_abortFlash("Allocation error in hgInit for prolongation")
  
  hg_restrict_n1 = max(NXB/2, 1)
  hg_restrict_n2 = max(NYB/2, 1)
  hg_restrict_n3 = max(NZB/2, 1)

  allocate(send_restrict_data(hg_restrict_n1,hg_restrict_n2,&
           hg_restrict_n3,nbbuf_restrict), stat=ierr1)
  allocate(recv_restrict_data(hg_restrict_n1, &
           hg_restrict_n2,hg_restrict_n3), stat=ierr2)
  allocate(send_restrict_req(nbbuf_restrict), stat=ierr3)
  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
       call Driver_abortFlash("Allocation error in hgInit for restriction")
  
  n1off(1) = 0
  n2off(1) = 0
  n3off(1) = 0
  n1off(2) = NXB/2
  n2off(2) = 0
  n3off(2) = 0
  n1off(3) = 0
  n2off(3) = NYB/2
  n3off(3) = 0
  n1off(4) = NXB/2
  n2off(4) = NYB/2
  n3off(4) = 0
  n1off(5) = 0
  n2off(5) = 0
  n3off(5) = NZB/2
  n1off(6) = NXB/2
  n2off(6) = 0
  n3off(6) = NZB/2
  n1off(7) = 0
  n2off(7) = NYB/2
  n3off(7) = NZB/2
  n1off(8) = NXB/2
  n2off(8) = NYB/2
  n3off(8) = NZB/2
  
  if (NDIM == 1) then
     
     Pew(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
     Pew(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
     
  endif
  
  if (NDIM == 2) then
     
     do c = 1, nchild
        do i = 1, NXB-1, 2
           Px(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
           Px(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
        enddo
        do i = 1, NYB-1, 2
           Py(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
           Py(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
        enddo
     enddo
     
     Pns(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
     Pns(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
     Pew(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
     Pew(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
  endif
  
  if (NDIM == 3) then
     
     do c = 1, nchild
        do i = 1, NXB-1, 2
           Px(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
           Px(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
        enddo
        do i = 1, NYB-1, 2
           Py(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
           Py(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
        enddo
        do i = 1, NZB-1, 2
           Pz(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
           Pz(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
        enddo
     enddo
     
     Pns(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
     Pns(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
     Pew(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
     Pew(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
     Pud(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
     Pud(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
     
  endif
  !======================================================================

  do i = 1, NXB
     hg_cx(:,i) = (/ 0., 1., -2., 1., 0. /)
     !!    hg_cx(:,i) = (/ -1./12., 4./3., -5./2., 4./3., -1./12. /)
  enddo
  !  hg_cx(:,1  ) = (/ 0., 4., -19./3., 8./3., -1./3. /)
  !  hg_cx(:,NXB) = (/ -5./3., 22./3., -41./3., 8., 0. /)
  do j = 1, NYB
     hg_cy(:,j) = (/ 0., 1., -2., 1., 0. /)
     !!    hg_cy(:,j) = (/ -1./12., 4./3., -5./2., 4./3., -1./12. /)
  enddo
  !  hg_cy(:,1  ) = (/ 0., 4., -19./3., 8./3., -1./3. /)
  !  hg_cy(:,NYB) = (/ -5./3., 22./3., -41./3., 8., 0. /)
  do k = 1, NZB
     hg_cz(:,k) = (/ 0., 1., -2., 1., 0. /)
     !!    hg_cz(:,k) = (/ -1./12., 4./3., -5./2., 4./3., -1./12. /)
  enddo
  !  hg_cz(:,1  ) = (/ 0., 4., -19./3., 8./3., -1./3. /)
  !  hg_cz(:,NZB) = (/ -5./3., 22./3., -41./3., 8., 0. /)


  call RuntimeParameters_get("mg_maxCorrections", gr_hgMaxCorrections)
  call RuntimeParameters_get("mg_maxResidualNorm", gr_hgMaxResidualNorm)
  call RuntimeParameters_get("mg_printNorm", gr_hgPrintNorm)
  

  nxl1 = 1
  nxl2 = max(NXB/2, 1)
  nxr1 = max(NXB/2+1, 1)
  nxr2 = NXB
  nyl1 = 1
  nyl2 = max(NYB/2, 1)
  nyr1 = max(NYB/2+1, 1)
  nyr2 = NYB
  nzl1 = 1
  nzl2 = max(NZB/2, 1)
  nzr1 = max(NZB/2+1, 1)
  nzr2 = NZB
  call Logfile_stamp("done","[gr_hgInit]")

  return
end subroutine gr_hgInit
