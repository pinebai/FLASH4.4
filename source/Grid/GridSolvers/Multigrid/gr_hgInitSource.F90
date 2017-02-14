!!****if* source/Grid/GridSolvers/Multigrid/gr_hgInitSource
!!
!! NAME
!!  gr_hgInitSource
!! 
!! SYNOPSIS
!!  gr_hgInitSource(integer, intent(in) :: isrc
!!                  integer, intent(in) :: isoln
!!
!! DESCRIPTION
!!  Initializes the source variable of the equation to be solved by
!!   the HG multigrid solver.  In this process it centers the variable
!!  if the boundary conditions are neumann or periodic, and applies
!!  external boundary values in the case of nonhomogenous dirichlet
!!  boundaries.  It also restricts the source function to all levels
!!
!!  Note that these operations are done on the source variable, so
!!  when using periodic or given-value boundaries a temporary variable
!!  should be used instead of the actual Grav_iSource
!!
!! ARGUMENTS
!!  isrc - the source variable of the equation
!!  isoln - the solution variable, which may contain given-value BCs.
!!
!! SIDE EFFECTS
!!
!!  Initializes module variable gr_hgAvgSource.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgInitSource (isrc, isoln)

  use tree, ONLY :lnblocks,neigh,nodetype,lrefine,bsize
  use gr_hgData, ONLY :gr_hgMeshRefineMax,gr_hgAvgSource
  use physicaldata, ONLY:  unk
  use gr_hgData, ONLY: hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui, &
       gr_hgBndTypes, gr_hgGeometry, gr_hgQuadrant
  use Grid_data, ONLY: gr_meshComm, gr_delta

  use Timers_interface, ONLY: Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgRestrict

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multigrid.h"
#include "Flash_mpi.h"

  integer, intent(in) :: isrc, isoln

  real,parameter:: coeff2 = 2.0 !13./6.
  integer       :: lb, i, j, k, ierr
  real          :: lvol, vol, lsum, bsum, sum
  real          :: nbinv, cvol, bvol
  real,dimension(MDIM) :: deltas


  !==========================================================================


  ! Need to restrict the source so it is valid at all levels.


  call Timers_start("gr_restrictTree")
  ! Note indices are different but call is the same
  do i = gr_hgMeshRefineMax-1, 1, -1
     call gr_hgRestrict(i+1, isrc, isrc)
  enddo
  call Timers_stop("gr_restrictTree")


     ! For given-value boundary conditions, we must subtract exterior boundary
     ! values (assumed to be stored in the first layer of boundary zones for
     ! the solution variable).
     
  if (ANY(gr_hgBndTypes(1:2*NDIM) == MG_BND_GIVENVAL)) then

     !! NOTE that results are wrong if the following loop gets applied to leaf blocks only. - KW
     do lb = 1, lnblocks
        deltas = gr_delta(1:MDIM,lrefine(lb))
        if ((gr_hgGeometry == MG_GEOM_1DCARTESIAN) .or. &
            (gr_hgGeometry == MG_GEOM_2DCARTESIAN) .or. &
            (gr_hgGeometry == MG_GEOM_3DCARTESIAN)) then

           if (gr_hgBndTypes(1)==MG_BND_GIVENVAL) then
              if (neigh(1,ILO_FACE,lb)<=PARAMESH_PHYSICAL_BOUNDARY) then
                 unk(isrc,NGUARD+1,:,:,lb) = &
                      unk(isrc,NGUARD+1,:,:,lb) - &
                      coeff2*unk(isoln,NGUARD,:,:,lb)/deltas(IAXIS)**2
              endif
           endif

        endif

        if ((gr_hgGeometry /= MG_GEOM_2DCYLAXISYM) .or. &
             (.not. gr_hgQuadrant)) then

           if ((NDIM >= 2) .and. (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              if (gr_hgBndTypes(3)==MG_BND_GIVENVAL) then
                 unk(isrc,:,(NGUARD*K2D)+1,:,lb) = &
                      unk(isrc,:,(NGUARD*K2D)+1,:,lb) - &
                      coeff2*unk(isoln,:,(NGUARD-1)*K2D+1,:,lb)/deltas(JAXIS)**2
              endif
           endif

        endif

        if ((NDIM == 3) .and. (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
           if (gr_hgBndTypes(5)==MG_BND_GIVENVAL) then
              unk(isrc,:,:,(NGUARD*K3D)+1,lb) = &
                   unk(isrc,:,:,(NGUARD*K3D)+1,lb) - &
                   coeff2*unk(isoln,:,:,((NGUARD-1)*K3D)+1,lb)/deltas(KAXIS)**2
           endif
        endif

        if (gr_hgBndTypes(2)==MG_BND_GIVENVAL) then
           if (neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then
              unk(isrc,NGUARD+NXB,:,:,lb) = &
                   unk(isrc,NGUARD+NXB,:,:,lb) - &
                   coeff2*unk(isoln,NGUARD+NXB+1,:,:,lb)/deltas(IAXIS)**2
           endif
        endif

        if ((NDIM >= 2) .and. (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
           if (gr_hgBndTypes(4)==MG_BND_GIVENVAL) then
              unk(isrc,:,(NGUARD*K2D)+NYB,:,lb) = &
                   unk(isrc,:,(NGUARD*K2D)+NYB,:,lb) - &
                   coeff2*unk(isoln,:,(NGUARD+NYB)*K2D+1,:,lb)/deltas(JAXIS)**2
           endif
        endif

        if ((NDIM == 3) .and. (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
           if (gr_hgBndTypes(6)==MG_BND_GIVENVAL) then
              unk(isrc,:,:,(NGUARD*K3D)+NZB,lb) = &
                   unk(isrc,:,:,(NGUARD*K3D)+NZB,lb) - &
                   coeff2*unk(isoln,:,:,((NGUARD+NZB)*K3D)+1,lb)/deltas(KAXIS)**2
           endif
        endif

     enddo

  endif

  ! For periodic or Neumann boundary conditions, compute the average value of
  ! the source and subtract it from the source.

  gr_hgAvgSource = 0.

  if (ALL((gr_hgBndTypes(1:2*NDIM) == MG_BND_PERIODIC) .or. &
          (gr_hgBndTypes(1:2*NDIM) == MG_BND_NEUMANN))) then
     
     lvol = 0.
     lsum = 0.
     nbinv = 1. / real(NXB)
     if (NDIM >= 2) nbinv = nbinv / real(NYB)
     if (NDIM == 3) nbinv = nbinv / real(NZB)

     do lb = 1, lnblocks
        if (lrefine(lb) == 1) then 
           bvol = bsize(IAXIS,lb)
           if (NDIM >= 2) bvol = bvol * bsize(JAXIS,lb)
           if (NDIM == 3) bvol = bvol * bsize(KAXIS,lb)
           cvol = bvol * nbinv
           lvol = lvol + bvol
           bsum = 0.
           do k = hg_kli, hg_kui
              do j = hg_jli, hg_jui
                 do i = hg_ili, hg_iui
                    bsum = bsum + unk(isrc,i,j,k,lb)
                 enddo
              enddo
           enddo
           lsum = lsum + bsum * cvol
        endif
     enddo
     
     call mpi_allreduce ( lsum, sum, 1, FLASH_REAL, & 
          MPI_SUM, gr_meshComm, ierr )
     call mpi_allreduce ( lvol, vol, 1, FLASH_REAL, & 
          MPI_SUM, gr_meshComm, ierr )

     gr_hgAvgSource = sum / vol
     
     do lb = 1, lnblocks
        unk(isrc,:,:,:,lb) = unk(isrc,:,:,:,lb) - gr_hgAvgSource
     enddo

  end if



  !=======================================================================

  return
end subroutine gr_hgInitSource
