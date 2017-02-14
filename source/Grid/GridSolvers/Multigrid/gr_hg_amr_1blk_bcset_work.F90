!!****if* source/Grid/GridSolvers/Multigrid/gr_hg_amr_1blk_bcset_work
!!
!! NAME
!!  gr_hg_amr_1blk_bcset_work
!!
!! SYNOPSIS
!! 
!!  call gr_hg_amr_1blk_bcset_work(integer(IN) :: lb, 
!!                           integer(IN) :: idest, 
!!                           integer(IN) :: idiag, 
!!                           integer(IN) :: idir)
!!
!! DESCRIPTION
!!
!!  Set boundary values of the exterior guardcells of a leaf block in the
!!  work array such that the given boundary value is satisfied.  The
!!  handled cases are MG_BND_DIRICHLET, MG_BND_NEUMANN, MG_BND_GIVENVAL
!! 
!!  Periodic boundaries are handled topologically and do not need a case
!!  here.  Allegedly this routine also handles cylindrical and spherical
!!  coordinates.
!!
!!  The implementation contains the same logic as the subroutine gr_hgSetExtBoundary,
!!  which is meant to work with PARAMESH2, except that here, boundary conditions
!!  are applied to only one block. This implementation fits somewhat better into the
!!  PARAMESH3/4 way of doing things.
!!
!!  The implementation is meant for comparisons of PARAMESH2 and PARAMESH3/4 results.
!!  In general, it would not be used when building Multigrid with PARAMESH 3 or 4.
!!
!! ARGUMENTS
!!
!!  lb     - block ID of the block to operate on
!!  idest  - selects the storage space in data_1blk.fh which is to
!!           be used in this call. If the leaf node is having its
!!           guardcells filled then this is set to 1, if its parent
!!           is being filled this is set to 2.
!!  idiag  - ignored variable for diagonal guardcell case
!!  idir   - sweep all directions (0), or x, y, or z only.
!!
!! RESULTS
!!
!!  Leaves some layers of guardcells in the work array for one block filled.
!!
!! NOTES
!!
!!  LBR has checked this (12/27/2007) and it looks OK for the cartesian case.  Meaning
!!   there are no obvious incorrect signs, etc.  But I don't understand how "work" relates to
!!   grid sizes.  Peter suspected there is a problem whenever the refinement boundaries
!!   are different sizes.
!!
!! NOTES
!!
!!  This version should also work with PARAMESH 2, but that is really unnecessary.
!!  When buildign Multigrid with PARAMESH2, gr_hgSetExtBoundary should be used instead.
!!
!! SEE ALSO
!!
!!  gr_hgSetExtBoundary   (this is derived from it)
!!  tot_bnd included with Paramesh2
!!
!!***

#include "Flash.h"

subroutine gr_hg_amr_1blk_bcset_work (lb, idest, &
     idiag, idir)

#ifdef FLASH_GRID_PARAMESH2
  use workspace, ONLY: work
#else
  use workspace, ONLY: work1
#endif
  use tree, ONLY : lrefine,lnblocks,neigh
  use gr_hgData, ONLY: hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui, &
       gr_hgQuadrant, gr_hgGeometry, gr_hgBndTypes, gr_hgCurrentGcReq_extrap

#include "Multigrid.h"
#include "constants.h"

  implicit none

  integer, intent(in) :: lb
  integer, intent(in) :: idest

  integer, intent(IN) :: idiag, idir
  logical :: extrap

  !===============================================================================
!!  extrap - if this .true., then the DIRICHLET fill is assumed NOT
!!           to go to zero at the exterior face.  This is important
!!           as massive oscillations may arise for given-value problems
!!           if the value goes to zero.
  if (idest == 2) then
     extrap = .TRUE.
  else
     extrap = gr_hgCurrentGcReq_extrap
  end if

#ifdef FLASH_GRID_PARAMESH2
#define work1 work
#define idest lb,1
#endif


  ! Directly set any external Dirichlet or Neumann boundaries.

  ! Dirichlet or Given-value boundary conditions ***************************************************************
  if ((gr_hgBndTypes(2*NDIM-1) == MG_BND_DIRICHLET) .or. &
       (gr_hgBndTypes(2*NDIM-1) == MG_BND_GIVENVAL)) then        


        if ((gr_hgGeometry == MG_GEOM_1DCARTESIAN) .or. &
            (gr_hgGeometry == MG_GEOM_2DCARTESIAN) .or. &
            (gr_hgGeometry == MG_GEOM_3DCARTESIAN)) then
           ! Cartesian
           if (neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

              if (extrap) then
                 work1(hg_ili-1,:,:,idest) = 2.*work1(hg_ili,:,:,idest) -    work1(hg_ili+1,:,:,idest)
                 work1(hg_ili-2,:,:,idest) = 3.*work1(hg_ili,:,:,idest) - 2.*work1(hg_ili+1,:,:,idest)
                 work1(hg_ili-3,:,:,idest) = 4.*work1(hg_ili,:,:,idest) - 3.*work1(hg_ili+1,:,:,idest)
                 work1(hg_ili-4,:,:,idest) = 5.*work1(hg_ili,:,:,idest) - 4.*work1(hg_ili+1,:,:,idest)
              else
                 work1(hg_ili-1,:,:,idest) = -   work1(hg_ili,:,:,idest)
                 work1(hg_ili-2,:,:,idest) = -3.*work1(hg_ili,:,:,idest)
                 work1(hg_ili-3,:,:,idest) = -5.*work1(hg_ili,:,:,idest)
                 work1(hg_ili-4,:,:,idest) = -7.*work1(hg_ili,:,:,idest)
              endif

           endif

        else
           ! Non-cartesian  
           if (neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work1(hg_ili-1,:,:,idest) = work1(hg_ili,:,:,idest)

        endif ! end of split between cartesian and non-cartesian

        if (neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

           if (extrap) then
              work1(hg_iui+1,:,:,idest) = 2.*work1(hg_iui,:,:,idest) -    work1(hg_iui-1,:,:,idest)
              work1(hg_iui+2,:,:,idest) = 3.*work1(hg_iui,:,:,idest) - 2.*work1(hg_iui-1,:,:,idest)
              work1(hg_iui+3,:,:,idest) = 4.*work1(hg_iui,:,:,idest) - 3.*work1(hg_iui-1,:,:,idest)
              work1(hg_iui+4,:,:,idest) = 5.*work1(hg_iui,:,:,idest) - 4.*work1(hg_iui-1,:,:,idest)
           else
              work1(hg_iui+1,:,:,idest) = -   work1(hg_iui,:,:,idest)
              work1(hg_iui+2,:,:,idest) = -3.*work1(hg_iui,:,:,idest)
              work1(hg_iui+3,:,:,idest) = -5.*work1(hg_iui,:,:,idest)
              work1(hg_iui+4,:,:,idest) = -7.*work1(hg_iui,:,:,idest)
           endif

           !DEV -- no upper face here for the non-cartesian case.  Seems wrong; now only Cartesian is legal
        endif

        ! Two-dimensional section @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (NDIM >= 2) then

           ! Begin Cartesian or polar  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
           if ((gr_hgGeometry /= MG_GEOM_2DCYLAXISYM) .or. (.not. gr_hgQuadrant)) then

              if (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then
                 if (extrap) then
                    work1(:,hg_jli-1,:,idest) = 2.*work1(:,hg_jli,:,idest) -    work1(:,hg_jli+1,:,idest)
                    work1(:,hg_jli-2,:,idest) = 3.*work1(:,hg_jli,:,idest) - 2.*work1(:,hg_jli+1,:,idest)
                    work1(:,hg_jli-3,:,idest) = 4.*work1(:,hg_jli,:,idest) - 3.*work1(:,hg_jli+1,:,idest)
                    work1(:,hg_jli-4,:,idest) = 5.*work1(:,hg_jli,:,idest) - 4.*work1(:,hg_jli+1,:,idest)
                 else
                    work1(:,hg_jli-1,:,idest) = -   work1(:,hg_jli,:,idest)
                    work1(:,hg_jli-2,:,idest) = -3.*work1(:,hg_jli,:,idest)
                    work1(:,hg_jli-3,:,idest) = -5.*work1(:,hg_jli,:,idest)
                    work1(:,hg_jli-4,:,idest) = -7.*work1(:,hg_jli,:,idest)
                 endif
              endif

              if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
                  (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
                 work1(hg_ili-1,hg_jli-1,:,idest) =  &
                        10.*work1(hg_ili,  hg_jli,  :,idest)   &  
                      - 10.*work1(hg_ili,  hg_jli+1,:,idest) &  
                      +  5.*work1(hg_ili,  hg_jli+2,:,idest) &  
                      -     work1(hg_ili,  hg_jli+3,:,idest) &
                      - 10.*work1(hg_ili+1,hg_jli,  :,idest) &
                      +  5.*work1(hg_ili+1,hg_jli+1,:,idest) & 
                      -     work1(hg_ili+1,hg_jli+2,:,idest) &
                      +  5.*work1(hg_ili+2,hg_jli,  :,idest) &
                      -     work1(hg_ili+2,hg_jli+1,:,idest) &
                      -     work1(hg_ili+3,hg_jli,  :,idest)

                 work1(hg_ili-1,hg_jli-2,:,idest) =  &
                        20.*work1(hg_ili,  hg_jli,  :,idest) &
                      - 30.*work1(hg_ili,  hg_jli+1,:,idest) &
                      + 18.*work1(hg_ili,  hg_jli+2,:,idest) &
                      -  4.*work1(hg_ili,  hg_jli+3,:,idest) &
                      - 15.*work1(hg_ili+1,hg_jli,  :,idest) &
                      + 12.*work1(hg_ili+1,hg_jli+1,:,idest) & 
                      -  3.*work1(hg_ili+1,hg_jli+2,:,idest) &
                      +  6.*work1(hg_ili+2,hg_jli,  :,idest) &
                      -  2.*work1(hg_ili+2,hg_jli+1,:,idest) &
                      -     work1(hg_ili+3,hg_jli,  :,idest)

                 work1(hg_ili-2,hg_jli-1,:,idest) =         &
                        20.*work1(hg_ili,  hg_jli,  :,idest) &
                      - 15.*work1(hg_ili,  hg_jli+1,:,idest) &
                      +  6.*work1(hg_ili,  hg_jli+2,:,idest) &
                      -     work1(hg_ili,  hg_jli+3,:,idest) &
                      - 30.*work1(hg_ili+1,hg_jli,  :,idest) &
                      + 12.*work1(hg_ili+1,hg_jli+1,:,idest) &
                      -  2.*work1(hg_ili+1,hg_jli+2,:,idest) &
                      + 18.*work1(hg_ili+2,hg_jli,  :,idest) &
                      -  3.*work1(hg_ili+2,hg_jli+1,:,idest) &
                      -  4.*work1(hg_ili+3,hg_jli,  :,idest)

                 work1(hg_ili-2,hg_jli-2,:,idest) =  &
                        35.*work1(hg_ili,  hg_jli,  :,idest) &
                      - 42.*work1(hg_ili,  hg_jli+1,:,idest) &
                      + 21.*work1(hg_ili,  hg_jli+2,:,idest) &
                      -  4.*work1(hg_ili,  hg_jli+3,:,idest) &
                      - 42.*work1(hg_ili+1,hg_jli,  :,idest) &
                      + 28.*work1(hg_ili+1,hg_jli+1,:,idest) &
                      -  6.*work1(hg_ili+1,hg_jli+2,:,idest) &
                      + 21.*work1(hg_ili+2,hg_jli,  :,idest) &
                      -  6.*work1(hg_ili+2,hg_jli+1,:,idest) &
                      -  4.*work1(hg_ili+3,hg_jli,  :,idest)
              endif

              if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
                  (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
                 work1(hg_iui+1,hg_jli-1,:,idest) =  &
                        10.*work1(hg_iui,  hg_jli,  :,idest) &
                      - 10.*work1(hg_iui,  hg_jli+1,:,idest) &
                      +  5.*work1(hg_iui,  hg_jli+2,:,idest) &
                      -     work1(hg_iui,  hg_jli+3,:,idest) &
                      - 10.*work1(hg_iui-1,hg_jli,  :,idest) &
                      +  5.*work1(hg_iui-1,hg_jli+1,:,idest) &
                      -     work1(hg_iui-1,hg_jli+2,:,idest) &
                      +  5.*work1(hg_iui-2,hg_jli,  :,idest) &
                      -     work1(hg_iui-2,hg_jli+1,:,idest) &
                      -     work1(hg_iui-3,hg_jli,  :,idest)

                 work1(hg_iui+1,hg_jli-2,:,idest) =  &
                        20.*work1(hg_iui,  hg_jli,  :,idest) &
                      - 30.*work1(hg_iui,  hg_jli+1,:,idest) & 
                      + 18.*work1(hg_iui,  hg_jli+2,:,idest) &
                      -  4.*work1(hg_iui,  hg_jli+3,:,idest) &
                      - 15.*work1(hg_iui-1,hg_jli,  :,idest) &
                      + 12.*work1(hg_iui-1,hg_jli+1,:,idest) &
                      -  3.*work1(hg_iui-1,hg_jli+2,:,idest) &
                      +  6.*work1(hg_iui-2,hg_jli,  :,idest) &
                      -  2.*work1(hg_iui-2,hg_jli+1,:,idest) &
                      -     work1(hg_iui-3,hg_jli,  :,idest)

                 work1(hg_iui+2,hg_jli-1,:,idest) =  &
                        20.*work1(hg_iui,  hg_jli,  :,idest) &
                      - 15.*work1(hg_iui,  hg_jli+1,:,idest) &
                      +  6.*work1(hg_iui,  hg_jli+2,:,idest) &
                      -     work1(hg_iui,  hg_jli+3,:,idest) &
                      - 30.*work1(hg_iui-1,hg_jli,  :,idest) &
                      + 12.*work1(hg_iui-1,hg_jli+1,:,idest) &
                      -  2.*work1(hg_iui-1,hg_jli+2,:,idest) &
                      + 18.*work1(hg_iui-2,hg_jli,  :,idest) &
                      -  3.*work1(hg_iui-2,hg_jli+1,:,idest) &
                      -  4.*work1(hg_iui-3,hg_jli,  :,idest)

                 work1(hg_iui+2,hg_jli-2,:,idest) =  &
                        35.*work1(hg_iui,  hg_jli,  :,idest) &
                      - 42.*work1(hg_iui,  hg_jli+1,:,idest) &
                      + 21.*work1(hg_iui,  hg_jli+2,:,idest) &
                      -  4.*work1(hg_iui,  hg_jli+3,:,idest) &
                      - 42.*work1(hg_iui-1,hg_jli,  :,idest) &
                      + 28.*work1(hg_iui-1,hg_jli+1,:,idest) &
                      -  6.*work1(hg_iui-1,hg_jli+2,:,idest) &
                      + 21.*work1(hg_iui-2,hg_jli,  :,idest) &
                      -  6.*work1(hg_iui-2,hg_jli+1,:,idest) &
                      -  4.*work1(hg_iui-3,hg_jli,  :,idest)
              endif


              ! End 2d cartesian or polar
           else
              ! 2d axisymmetric
              if (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                   work1(:,hg_jli-1,:,idest) = work1(:,hg_jli,:,idest)
           endif
           !  End of split between cartesian / polar / axisymmetric.  From here all 2d cases are done

           if (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

              if (extrap) then
                 work1(:,hg_jui+1,:,idest) = 2.*work1(:,hg_jui,:,idest) -    work1(:,hg_jui-1,:,idest)
                 work1(:,hg_jui+2,:,idest) = 3.*work1(:,hg_jui,:,idest) - 2.*work1(:,hg_jui-1,:,idest)
                 work1(:,hg_jui+3,:,idest) = 4.*work1(:,hg_jui,:,idest) - 3.*work1(:,hg_jui-1,:,idest)
                 work1(:,hg_jui+4,:,idest) = 5.*work1(:,hg_jui,:,idest) - 4.*work1(:,hg_jui-1,:,idest)
              else
                 work1(:,hg_jui+1,:,idest) = -   work1(:,hg_jui,:,idest)
                 work1(:,hg_jui+2,:,idest) = -3.*work1(:,hg_jui,:,idest)
                 work1(:,hg_jui+3,:,idest) = -5.*work1(:,hg_jui,:,idest)
                 work1(:,hg_jui+4,:,idest) = -7.*work1(:,hg_jui,:,idest)
              endif

           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,hg_jui+1,:,idest) =  &
                     10.*work1(hg_ili,  hg_jui,  :,idest) &
                   - 10.*work1(hg_ili,  hg_jui-1,:,idest) &
                   +  5.*work1(hg_ili,  hg_jui-2,:,idest) &
                   -     work1(hg_ili,  hg_jui-3,:,idest) &
                   - 10.*work1(hg_ili+1,hg_jui,  :,idest) &
                   +  5.*work1(hg_ili+1,hg_jui-1,:,idest) &
                   -     work1(hg_ili+1,hg_jui-2,:,idest) &
                   +  5.*work1(hg_ili+2,hg_jui,  :,idest) &
                   -     work1(hg_ili+2,hg_jui-1,:,idest) &
                   -     work1(hg_ili+3,hg_jui,  :,idest)

              work1(hg_ili-1,hg_jui+2,:,idest) =  &
                     20.*work1(hg_ili,  hg_jui,  :,idest) &
                   - 30.*work1(hg_ili,  hg_jui-1,:,idest) &
                   + 18.*work1(hg_ili,  hg_jui-2,:,idest) &
                   -  4.*work1(hg_ili,  hg_jui-3,:,idest) &
                   - 15.*work1(hg_ili+1,hg_jui,  :,idest) &
                   + 12.*work1(hg_ili+1,hg_jui-1,:,idest) &
                   -  3.*work1(hg_ili+1,hg_jui-2,:,idest) &
                   +  6.*work1(hg_ili+2,hg_jui,  :,idest) &
                   -  2.*work1(hg_ili+2,hg_jui-1,:,idest) &
                   -     work1(hg_ili+3,hg_jui,  :,idest)

              work1(hg_ili-2,hg_jui+1,:,idest) =  &
                     20.*work1(hg_ili,  hg_jui,  :,idest) &
                   - 15.*work1(hg_ili,  hg_jui-1,:,idest) &
                   +  6.*work1(hg_ili,  hg_jui-2,:,idest) &
                   -     work1(hg_ili,  hg_jui-3,:,idest) &
                   - 30.*work1(hg_ili+1,hg_jui,  :,idest) &
                   + 12.*work1(hg_ili+1,hg_jui-1,:,idest) &
                   -  2.*work1(hg_ili+1,hg_jui-2,:,idest) &
                   + 18.*work1(hg_ili+2,hg_jui,  :,idest) &
                   -  3.*work1(hg_ili+2,hg_jui-1,:,idest) &
                   -  4.*work1(hg_ili+3,hg_jui,  :,idest)

              work1(hg_ili-2,hg_jui+2,:,idest) =  &
                     35.*work1(hg_ili,  hg_jui,  :,idest) & 
                   - 42.*work1(hg_ili,  hg_jui-1,:,idest) &
                   + 21.*work1(hg_ili,  hg_jui-2,:,idest) &
                   -  4.*work1(hg_ili,  hg_jui-3,:,idest) &
                   - 42.*work1(hg_ili+1,hg_jui,  :,idest) &
                   + 28.*work1(hg_ili+1,hg_jui-1,:,idest) &
                   -  6.*work1(hg_ili+1,hg_jui-2,:,idest) &
                   + 21.*work1(hg_ili+2,hg_jui,  :,idest) &
                   -  6.*work1(hg_ili+2,hg_jui-1,:,idest)&
                   -  4.*work1(hg_ili+3,hg_jui,  :,idest)
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,hg_jui+1,:,idest) =  &
                     10.*work1(hg_iui,  hg_jui,  :,idest) &
                   - 10.*work1(hg_iui,  hg_jui-1,:,idest) &
                   +  5.*work1(hg_iui,  hg_jui-2,:,idest) &
                   -     work1(hg_iui,  hg_jui-3,:,idest) &
                   - 10.*work1(hg_iui-1,hg_jui,  :,idest) &
                   +  5.*work1(hg_iui-1,hg_jui-1,:,idest) &
                   -     work1(hg_iui-1,hg_jui-2,:,idest) &
                   +  5.*work1(hg_iui-2,hg_jui,  :,idest) &
                   -     work1(hg_iui-2,hg_jui-1,:,idest) &
                   -     work1(hg_iui-3,hg_jui,  :,idest)

              work1(hg_iui+1,hg_jui+2,:,idest) =  &
                     20.*work1(hg_iui,  hg_jui,  :,idest) &
                   - 30.*work1(hg_iui,  hg_jui-1,:,idest) &
                   + 18.*work1(hg_iui,  hg_jui-2,:,idest) &
                   -  4.*work1(hg_iui,  hg_jui-3,:,idest) &
                   - 15.*work1(hg_iui-1,hg_jui,  :,idest) &
                   + 12.*work1(hg_iui-1,hg_jui-1,:,idest) &
                   -  3.*work1(hg_iui-1,hg_jui-2,:,idest) &
                   +  6.*work1(hg_iui-2,hg_jui,  :,idest) &
                   -  2.*work1(hg_iui-2,hg_jui-1,:,idest) &
                   -     work1(hg_iui-3,hg_jui,  :,idest)

              work1(hg_iui+2,hg_jui+1,:,idest) =  &
                     20.*work1(hg_iui,  hg_jui,  :,idest) &
                   - 15.*work1(hg_iui,  hg_jui-1,:,idest) &
                   +  6.*work1(hg_iui,  hg_jui-2,:,idest) &
                   -     work1(hg_iui,  hg_jui-3,:,idest) &
                   - 30.*work1(hg_iui-1,hg_jui,  :,idest) &
                   + 12.*work1(hg_iui-1,hg_jui-1,:,idest) &
                   -  2.*work1(hg_iui-1,hg_jui-2,:,idest) &
                   + 18.*work1(hg_iui-2,hg_jui,  :,idest) &
                   -  3.*work1(hg_iui-2,hg_jui-1,:,idest) &
                   -  4.*work1(hg_iui-3,hg_jui,  :,idest)

              work1(hg_iui+2,hg_jui+2,:,idest) =  &
                     35.*work1(hg_iui,  hg_jui,  :,idest) &
                   - 42.*work1(hg_iui,  hg_jui-1,:,idest) &
                   + 21.*work1(hg_iui,  hg_jui-2,:,idest) &
                   -  4.*work1(hg_iui,  hg_jui-3,:,idest) &
                   - 42.*work1(hg_iui-1,hg_jui,  :,idest) &
                   + 28.*work1(hg_iui-1,hg_jui-1,:,idest) &
                   -  6.*work1(hg_iui-1,hg_jui-2,:,idest) &
                   + 21.*work1(hg_iui-2,hg_jui,  :,idest) &
                   -  6.*work1(hg_iui-2,hg_jui-1,:,idest) &
                   -  4.*work1(hg_iui-3,hg_jui,  :,idest)
           endif

        endif
        ! End of two-dimensional section @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ! Start of three-dimensional section @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (NDIM == 3) then

           if (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

              if (extrap) then
                 work1(:,:,hg_kli-1,idest) = 2.*work1(:,:,hg_kli,idest) -    work1(:,:,hg_kli+1,idest)
                 work1(:,:,hg_kli-2,idest) = 3.*work1(:,:,hg_kli,idest) - 2.*work1(:,:,hg_kli+1,idest)
                 work1(:,:,hg_kli-3,idest) = 4.*work1(:,:,hg_kli,idest) - 3.*work1(:,:,hg_kli+1,idest)
                 work1(:,:,hg_kli-4,idest) = 5.*work1(:,:,hg_kli,idest) - 4.*work1(:,:,hg_kli+1,idest)
              else
                 work1(:,:,hg_kli-1,idest) = -   work1(:,:,hg_kli,idest)
                 work1(:,:,hg_kli-2,idest) = -3.*work1(:,:,hg_kli,idest)
                 work1(:,:,hg_kli-3,idest) = -5.*work1(:,:,hg_kli,idest)
                 work1(:,:,hg_kli-4,idest) = -7.*work1(:,:,hg_kli,idest)
              endif

           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,:,hg_kli-1,idest) = &
                   +10*work1(hg_ili+0,:,hg_kli+0,idest)-10*work1(hg_ili+0,:,hg_kli+1,idest) &
                   + 5*work1(hg_ili+0,:,hg_kli+2,idest) -1*work1(hg_ili+0,:,hg_kli+3,idest) &
                   -10*work1(hg_ili+1,:,hg_kli+0,idest) +5*work1(hg_ili+1,:,hg_kli+1,idest) &
                   -1 *work1(hg_ili+1,:,hg_kli+2,idest) +5*work1(hg_ili+2,:,hg_kli+0,idest) &
                   -1 *work1(hg_ili+2,:,hg_kli+1,idest) -1*work1(hg_ili+3,:,hg_kli+0,idest)

              work1(hg_ili-1,:,hg_kli-2,idest) = &
                   +20*work1(hg_ili+0,:,hg_kli+0,idest)-30*work1(hg_ili+0,:,hg_kli+1,idest) &
                   +18*work1(hg_ili+0,:,hg_kli+2,idest) -4*work1(hg_ili+0,:,hg_kli+3,idest) &
                   -15*work1(hg_ili+1,:,hg_kli+0,idest)+12*work1(hg_ili+1,:,hg_kli+1,idest) &
                   -3 *work1(hg_ili+1,:,hg_kli+2,idest) +6*work1(hg_ili+2,:,hg_kli+0,idest) &
                   -2 *work1(hg_ili+2,:,hg_kli+1,idest) -1*work1(hg_ili+3,:,hg_kli+0,idest)

              work1(hg_ili-2,:,hg_kli-1,idest) = &
                   +20*work1(hg_ili+0,:,hg_kli+0,idest)-15*work1(hg_ili+0,:,hg_kli+1,idest) &
                   +6 *work1(hg_ili+0,:,hg_kli+2,idest) -1*work1(hg_ili+0,:,hg_kli+3,idest) &
                   -30*work1(hg_ili+1,:,hg_kli+0,idest)+12*work1(hg_ili+1,:,hg_kli+1,idest) &
                   -2 *work1(hg_ili+1,:,hg_kli+2,idest)+18*work1(hg_ili+2,:,hg_kli+0,idest) &
                   -3 *work1(hg_ili+2,:,hg_kli+1,idest) -4*work1(hg_ili+3,:,hg_kli+0,idest)

              work1(hg_ili-2,:,hg_kli-2,idest) = &
                   +35*work1(hg_ili+0,:,hg_kli+0,idest)-42*work1(hg_ili+0,:,hg_kli+1,idest) &
                   +21*work1(hg_ili+0,:,hg_kli+2,idest) -4*work1(hg_ili+0,:,hg_kli+3,idest) &
                   -42*work1(hg_ili+1,:,hg_kli+0,idest)+28*work1(hg_ili+1,:,hg_kli+1,idest) &
                   -6 *work1(hg_ili+1,:,hg_kli+2,idest)+21*work1(hg_ili+2,:,hg_kli+0,idest) &
                   -6 *work1(hg_ili+2,:,hg_kli+1,idest) -4*work1(hg_ili+3,:,hg_kli+0,idest)
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,:,hg_kli-1,idest) = &
                   +10*work1(hg_iui-0,:,hg_kli+0,idest)-10*work1(hg_iui-0,:,hg_kli+1,idest) &
                   +5 *work1(hg_iui-0,:,hg_kli+2,idest) -1*work1(hg_iui-0,:,hg_kli+3,idest) &
                   -10*work1(hg_iui-1,:,hg_kli+0,idest) +5*work1(hg_iui-1,:,hg_kli+1,idest) &
                   -1 *work1(hg_iui-1,:,hg_kli+2,idest) +5*work1(hg_iui-2,:,hg_kli+0,idest) &
                   -1 *work1(hg_iui-2,:,hg_kli+1,idest) -1*work1(hg_iui-3,:,hg_kli+0,idest)

              work1(hg_iui+1,:,hg_kli-2,idest) = &
                   +20*work1(hg_iui-0,:,hg_kli+0,idest)-30*work1(hg_iui-0,:,hg_kli+1,idest) &
                   +18*work1(hg_iui-0,:,hg_kli+2,idest) -4*work1(hg_iui-0,:,hg_kli+3,idest) &
                   -15*work1(hg_iui-1,:,hg_kli+0,idest)+12*work1(hg_iui-1,:,hg_kli+1,idest) &
                   -3 *work1(hg_iui-1,:,hg_kli+2,idest) +6*work1(hg_iui-2,:,hg_kli+0,idest) &
                   -2 *work1(hg_iui-2,:,hg_kli+1,idest) -1*work1(hg_iui-3,:,hg_kli+0,idest)

              work1(hg_iui+2,:,hg_kli-1,idest) = &
                   +20*work1(hg_iui-0,:,hg_kli+0,idest)-15*work1(hg_iui-0,:,hg_kli+1,idest) &
                   +6 *work1(hg_iui-0,:,hg_kli+2,idest) -1*work1(hg_iui-0,:,hg_kli+3,idest) &
                   -30*work1(hg_iui-1,:,hg_kli+0,idest)+12*work1(hg_iui-1,:,hg_kli+1,idest) &
                   -2 *work1(hg_iui-1,:,hg_kli+2,idest)+18*work1(hg_iui-2,:,hg_kli+0,idest) &
                   -3 *work1(hg_iui-2,:,hg_kli+1,idest) -4*work1(hg_iui-3,:,hg_kli+0,idest)

              work1(hg_iui+2,:,hg_kli-2,idest) = &
                   +35*work1(hg_iui-0,:,hg_kli+0,idest)-42*work1(hg_iui-0,:,hg_kli+1,idest) &
                   +21*work1(hg_iui-0,:,hg_kli+2,idest) -4*work1(hg_iui-0,:,hg_kli+3,idest) &
                   -42*work1(hg_iui-1,:,hg_kli+0,idest)+28*work1(hg_iui-1,:,hg_kli+1,idest) &
                   -6 *work1(hg_iui-1,:,hg_kli+2,idest)+21*work1(hg_iui-2,:,hg_kli+0,idest) &
                   -6 *work1(hg_iui-2,:,hg_kli+1,idest) -4*work1(hg_iui-3,:,hg_kli+0,idest)
           endif

           if ((neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(:,hg_jli-1,hg_kli-1,idest) = &
                   +10*work1(:,hg_jli+0,hg_kli+0,idest)-10*work1(:,hg_jli+0,hg_kli+1,idest) &
                   +5 *work1(:,hg_jli+0,hg_kli+2,idest) -1*work1(:,hg_jli+0,hg_kli+3,idest) &
                   -10*work1(:,hg_jli+1,hg_kli+0,idest) +5*work1(:,hg_jli+1,hg_kli+1,idest) &
                   -1 *work1(:,hg_jli+1,hg_kli+2,idest) +5*work1(:,hg_jli+2,hg_kli+0,idest) &
                   -1 *work1(:,hg_jli+2,hg_kli+1,idest) -1*work1(:,hg_jli+3,hg_kli+0,idest)

              work1(:,hg_jli-1,hg_kli-2,idest) = &
                   +20*work1(:,hg_jli+0,hg_kli+0,idest)-30*work1(:,hg_jli+0,hg_kli+1,idest) &
                   +18*work1(:,hg_jli+0,hg_kli+2,idest) -4*work1(:,hg_jli+0,hg_kli+3,idest) &
                   -15*work1(:,hg_jli+1,hg_kli+0,idest)+12*work1(:,hg_jli+1,hg_kli+1,idest) &
                   -3 *work1(:,hg_jli+1,hg_kli+2,idest) +6*work1(:,hg_jli+2,hg_kli+0,idest) &
                   -2 *work1(:,hg_jli+2,hg_kli+1,idest) -1*work1(:,hg_jli+3,hg_kli+0,idest)

              work1(:,hg_jli-2,hg_kli-1,idest) = &
                   +20*work1(:,hg_jli+0,hg_kli+0,idest)-15*work1(:,hg_jli+0,hg_kli+1,idest) &
                   +6 *work1(:,hg_jli+0,hg_kli+2,idest) -1*work1(:,hg_jli+0,hg_kli+3,idest) &
                   -30*work1(:,hg_jli+1,hg_kli+0,idest)+12*work1(:,hg_jli+1,hg_kli+1,idest) &
                   -2 *work1(:,hg_jli+1,hg_kli+2,idest)+18*work1(:,hg_jli+2,hg_kli+0,idest) &
                   -3 *work1(:,hg_jli+2,hg_kli+1,idest) -4*work1(:,hg_jli+3,hg_kli+0,idest)

              work1(:,hg_jli-2,hg_kli-2,idest) = &
                   +35*work1(:,hg_jli+0,hg_kli+0,idest)-42*work1(:,hg_jli+0,hg_kli+1,idest) &
                   +21*work1(:,hg_jli+0,hg_kli+2,idest) -4*work1(:,hg_jli+0,hg_kli+3,idest) &
                   -42*work1(:,hg_jli+1,hg_kli+0,idest)+28*work1(:,hg_jli+1,hg_kli+1,idest) &
                   -6 *work1(:,hg_jli+1,hg_kli+2,idest)+21*work1(:,hg_jli+2,hg_kli+0,idest) &
                   -6 *work1(:,hg_jli+2,hg_kli+1,idest) -4*work1(:,hg_jli+3,hg_kli+0,idest)
           endif

           if ((neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(:,hg_jui+1,hg_kli-1,idest) = &
                   +10*work1(:,hg_jui-0,hg_kli+0,idest)-10*work1(:,hg_jui-0,hg_kli+1,idest) &
                   + 5*work1(:,hg_jui-0,hg_kli+2,idest) -1*work1(:,hg_jui-0,hg_kli+3,idest) &
                   -10*work1(:,hg_jui-1,hg_kli+0,idest) +5*work1(:,hg_jui-1,hg_kli+1,idest) &
                   -1 *work1(:,hg_jui-1,hg_kli+2,idest) +5*work1(:,hg_jui-2,hg_kli+0,idest) &
                   -1 *work1(:,hg_jui-2,hg_kli+1,idest) -1*work1(:,hg_jui-3,hg_kli+0,idest)

              work1(:,hg_jui+1,hg_kli-2,idest) = &
                   +20*work1(:,hg_jui-0,hg_kli+0,idest)-30*work1(:,hg_jui-0,hg_kli+1,idest) &
                   +18*work1(:,hg_jui-0,hg_kli+2,idest) -4*work1(:,hg_jui-0,hg_kli+3,idest) &
                   -15*work1(:,hg_jui-1,hg_kli+0,idest)+12*work1(:,hg_jui-1,hg_kli+1,idest) &
                   -3 *work1(:,hg_jui-1,hg_kli+2,idest) +6*work1(:,hg_jui-2,hg_kli+0,idest) &
                   -2 *work1(:,hg_jui-2,hg_kli+1,idest) -1*work1(:,hg_jui-3,hg_kli+0,idest)

              work1(:,hg_jui+2,hg_kli-1,idest) = &
                   +20*work1(:,hg_jui-0,hg_kli+0,idest)-15*work1(:,hg_jui-0,hg_kli+1,idest) &
                   +6 *work1(:,hg_jui-0,hg_kli+2,idest) -1*work1(:,hg_jui-0,hg_kli+3,idest) &
                   -30*work1(:,hg_jui-1,hg_kli+0,idest)+12*work1(:,hg_jui-1,hg_kli+1,idest) &
                   -2 *work1(:,hg_jui-1,hg_kli+2,idest)+18*work1(:,hg_jui-2,hg_kli+0,idest) &
                   -3 *work1(:,hg_jui-2,hg_kli+1,idest) -4*work1(:,hg_jui-3,hg_kli+0,idest)

              work1(:,hg_jui+2,hg_kli-2,idest) = &
                   +35*work1(:,hg_jui-0,hg_kli+0,idest)-42*work1(:,hg_jui-0,hg_kli+1,idest) &
                   +21*work1(:,hg_jui-0,hg_kli+2,idest) -4*work1(:,hg_jui-0,hg_kli+3,idest) &
                   -42*work1(:,hg_jui-1,hg_kli+0,idest)+28*work1(:,hg_jui-1,hg_kli+1,idest) &
                   -6 *work1(:,hg_jui-1,hg_kli+2,idest)+21*work1(:,hg_jui-2,hg_kli+0,idest) &
                   -6 *work1(:,hg_jui-2,hg_kli+1,idest) -4*work1(:,hg_jui-3,hg_kli+0,idest)
           endif

           ! Should use the correct 2nd-order expressions for external corners, but for
           ! now we are just setting these cells to zero.
           !  Lower faces of K, high ones are below
           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,hg_jli-1,hg_kli-1,idest) = 0.
              work1(hg_ili-2,hg_jli-1,hg_kli-1,idest) = 0.
              work1(hg_ili-1,hg_jli-2,hg_kli-1,idest) = 0.
              work1(hg_ili-2,hg_jli-2,hg_kli-1,idest) = 0.
              work1(hg_ili-1,hg_jli-1,hg_kli-2,idest) = 0.
              work1(hg_ili-2,hg_jli-1,hg_kli-2,idest) = 0.
              work1(hg_ili-1,hg_jli-2,hg_kli-2,idest) = 0.
              work1(hg_ili-2,hg_jli-2,hg_kli-2,idest) = 0.
           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,hg_jui+1,hg_kli-1,idest) = 0.
              work1(hg_ili-2,hg_jui+1,hg_kli-1,idest) = 0.
              work1(hg_ili-1,hg_jui+2,hg_kli-1,idest) = 0.
              work1(hg_ili-2,hg_jui+2,hg_kli-1,idest) = 0.
              work1(hg_ili-1,hg_jui+1,hg_kli-2,idest) = 0.
              work1(hg_ili-2,hg_jui+1,hg_kli-2,idest) = 0.
              work1(hg_ili-1,hg_jui+2,hg_kli-2,idest) = 0.
              work1(hg_ili-2,hg_jui+2,hg_kli-2,idest) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,hg_jli-1,hg_kli-1,idest) = 0.
              work1(hg_iui+2,hg_jli-1,hg_kli-1,idest) = 0.
              work1(hg_iui+1,hg_jli-2,hg_kli-1,idest) = 0.
              work1(hg_iui+2,hg_jli-2,hg_kli-1,idest) = 0.
              work1(hg_iui+1,hg_jli-1,hg_kli-2,idest) = 0.
              work1(hg_iui+2,hg_jli-1,hg_kli-2,idest) = 0.
              work1(hg_iui+1,hg_jli-2,hg_kli-2,idest) = 0.
              work1(hg_iui+2,hg_jli-2,hg_kli-2,idest) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,hg_jui+1,hg_kli-1,idest) = 0.
              work1(hg_iui+2,hg_jui+1,hg_kli-1,idest) = 0.
              work1(hg_iui+1,hg_jui+2,hg_kli-1,idest) = 0.
              work1(hg_iui+2,hg_jui+2,hg_kli-1,idest) = 0.
              work1(hg_iui+1,hg_jui+1,hg_kli-2,idest) = 0.
              work1(hg_iui+2,hg_jui+1,hg_kli-2,idest) = 0.
              work1(hg_iui+1,hg_jui+2,hg_kli-2,idest) = 0.
              work1(hg_iui+2,hg_jui+2,hg_kli-2,idest) = 0.
           endif


           if (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then
              if (extrap) then
                 work1(:,:,hg_kui+1,idest) = 2.*work1(:,:,hg_kui,idest) -    work1(:,:,hg_kui-1,idest)
                 work1(:,:,hg_kui+2,idest) = 3.*work1(:,:,hg_kui,idest) - 2.*work1(:,:,hg_kui-1,idest)
                 work1(:,:,hg_kui+3,idest) = 4.*work1(:,:,hg_kui,idest) - 3.*work1(:,:,hg_kui-1,idest)
                 work1(:,:,hg_kui+4,idest) = 5.*work1(:,:,hg_kui,idest) - 4.*work1(:,:,hg_kui-1,idest)
              else
                 work1(:,:,hg_kui+1,idest) = -   work1(:,:,hg_kui,idest)
                 work1(:,:,hg_kui+2,idest) = -3.*work1(:,:,hg_kui,idest)
                 work1(:,:,hg_kui+3,idest) = -5.*work1(:,:,hg_kui,idest)
                 work1(:,:,hg_kui+4,idest) = -7.*work1(:,:,hg_kui,idest)
              endif
           endif

           ! DEV could nest these if statements within the previous one, minor savings in execution
           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,:,hg_kui+1,idest) = &
                   +10*work1(hg_ili+0,:,hg_kui-0,idest)-10*work1(hg_ili+0,:,hg_kui-1,idest) &
                   +5 *work1(hg_ili+0,:,hg_kui-2,idest) -1*work1(hg_ili+0,:,hg_kui-3,idest) &
                   -10*work1(hg_ili+1,:,hg_kui-0,idest) +5*work1(hg_ili+1,:,hg_kui-1,idest) &
                   -1 *work1(hg_ili+1,:,hg_kui-2,idest) +5*work1(hg_ili+2,:,hg_kui-0,idest) &
                   -1 *work1(hg_ili+2,:,hg_kui-1,idest) -1*work1(hg_ili+3,:,hg_kui-0,idest)

              work1(hg_ili-1,:,hg_kui+2,idest) = &
                   +20*work1(hg_ili+0,:,hg_kui-0,idest)-30*work1(hg_ili+0,:,hg_kui-1,idest) &
                   +18*work1(hg_ili+0,:,hg_kui-2,idest) -4*work1(hg_ili+0,:,hg_kui-3,idest) &
                   -15*work1(hg_ili+1,:,hg_kui-0,idest)+12*work1(hg_ili+1,:,hg_kui-1,idest) &
                   -3 *work1(hg_ili+1,:,hg_kui-2,idest) +6*work1(hg_ili+2,:,hg_kui-0,idest) &
                   -2 *work1(hg_ili+2,:,hg_kui-1,idest) -1*work1(hg_ili+3,:,hg_kui-0,idest)

              work1(hg_ili-2,:,hg_kui+1,idest) = &
                   +20*work1(hg_ili+0,:,hg_kui-0,idest)-15*work1(hg_ili+0,:,hg_kui-1,idest) &
                   +6 *work1(hg_ili+0,:,hg_kui-2,idest) -1*work1(hg_ili+0,:,hg_kui-3,idest) &
                   -30*work1(hg_ili+1,:,hg_kui-0,idest)+12*work1(hg_ili+1,:,hg_kui-1,idest) &
                   -2 *work1(hg_ili+1,:,hg_kui-2,idest)+18*work1(hg_ili+2,:,hg_kui-0,idest) &
                   -3 *work1(hg_ili+2,:,hg_kui-1,idest) -4*work1(hg_ili+3,:,hg_kui-0,idest)

              work1(hg_ili-2,:,hg_kui+2,idest) = &
                   +35*work1(hg_ili+0,:,hg_kui-0,idest)-42*work1(hg_ili+0,:,hg_kui-1,idest) &
                   +21*work1(hg_ili+0,:,hg_kui-2,idest) -4*work1(hg_ili+0,:,hg_kui-3,idest) &
                   -42*work1(hg_ili+1,:,hg_kui-0,idest)+28*work1(hg_ili+1,:,hg_kui-1,idest) &
                   -6 *work1(hg_ili+1,:,hg_kui-2,idest)+21*work1(hg_ili+2,:,hg_kui-0,idest) &
                   -6 *work1(hg_ili+2,:,hg_kui-1,idest) -4*work1(hg_ili+3,:,hg_kui-0,idest)
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,:,hg_kui+1,idest) = &
                   +10*work1(hg_iui-0,:,hg_kui-0,idest)-10*work1(hg_iui-0,:,hg_kui-1,idest) &
                   +5 *work1(hg_iui-0,:,hg_kui-2,idest) -1*work1(hg_iui-0,:,hg_kui-3,idest) &
                   -10*work1(hg_iui-1,:,hg_kui-0,idest) +5*work1(hg_iui-1,:,hg_kui-1,idest) &
                   -1 *work1(hg_iui-1,:,hg_kui-2,idest) +5*work1(hg_iui-2,:,hg_kui-0,idest) &
                   -1 *work1(hg_iui-2,:,hg_kui-1,idest) -1*work1(hg_iui-3,:,hg_kui-0,idest)

              work1(hg_iui+1,:,hg_kui+2,idest) = &
                   +20*work1(hg_iui-0,:,hg_kui-0,idest)-30*work1(hg_iui-0,:,hg_kui-1,idest) &
                   +18*work1(hg_iui-0,:,hg_kui-2,idest) -4*work1(hg_iui-0,:,hg_kui-3,idest) &
                   -15*work1(hg_iui-1,:,hg_kui-0,idest)+12*work1(hg_iui-1,:,hg_kui-1,idest) &
                   -3 *work1(hg_iui-1,:,hg_kui-2,idest) +6*work1(hg_iui-2,:,hg_kui-0,idest) &
                   -2 *work1(hg_iui-2,:,hg_kui-1,idest) -1*work1(hg_iui-3,:,hg_kui-0,idest)

              work1(hg_iui+2,:,hg_kui+1,idest) = &
                   +20*work1(hg_iui-0,:,hg_kui-0,idest)-15*work1(hg_iui-0,:,hg_kui-1,idest) &
                   +6 *work1(hg_iui-0,:,hg_kui-2,idest) -1*work1(hg_iui-0,:,hg_kui-3,idest) &
                   -30*work1(hg_iui-1,:,hg_kui-0,idest)+12*work1(hg_iui-1,:,hg_kui-1,idest) &
                   -2 *work1(hg_iui-1,:,hg_kui-2,idest)+18*work1(hg_iui-2,:,hg_kui-0,idest) &
                   -3 *work1(hg_iui-2,:,hg_kui-1,idest) -4*work1(hg_iui-3,:,hg_kui-0,idest)

              work1(hg_iui+2,:,hg_kui+2,idest) = &
                   +35*work1(hg_iui-0,:,hg_kui-0,idest)-42*work1(hg_iui-0,:,hg_kui-1,idest) &
                   +21*work1(hg_iui-0,:,hg_kui-2,idest) -4*work1(hg_iui-0,:,hg_kui-3,idest) &
                   -42*work1(hg_iui-1,:,hg_kui-0,idest)+28*work1(hg_iui-1,:,hg_kui-1,idest) &
                   -6 *work1(hg_iui-1,:,hg_kui-2,idest)+21*work1(hg_iui-2,:,hg_kui-0,idest) &
                   -6 *work1(hg_iui-2,:,hg_kui-1,idest) -4*work1(hg_iui-3,:,hg_kui-0,idest)
           endif

           if ((neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(:,hg_jli-1,hg_kui+1,idest) = &
                   +10*work1(:,hg_jli+0,hg_kui-0,idest)-10*work1(:,hg_jli+0,hg_kui-1,idest) &
                   +5 *work1(:,hg_jli+0,hg_kui-2,idest) -1*work1(:,hg_jli+0,hg_kui-3,idest) &
                   -10*work1(:,hg_jli+1,hg_kui-0,idest) +5*work1(:,hg_jli+1,hg_kui-1,idest) &
                   -1 *work1(:,hg_jli+1,hg_kui-2,idest) +5*work1(:,hg_jli+2,hg_kui-0,idest) &
                   -1 *work1(:,hg_jli+2,hg_kui-1,idest) -1*work1(:,hg_jli+3,hg_kui-0,idest)

              work1(:,hg_jli-1,hg_kui+2,idest) = &
                   +20*work1(:,hg_jli+0,hg_kui-0,idest)-30*work1(:,hg_jli+0,hg_kui-1,idest) &
                   +18*work1(:,hg_jli+0,hg_kui-2,idest) -4*work1(:,hg_jli+0,hg_kui-3,idest) &
                   -15*work1(:,hg_jli+1,hg_kui-0,idest)+12*work1(:,hg_jli+1,hg_kui-1,idest) &
                   -3 *work1(:,hg_jli+1,hg_kui-2,idest) +6*work1(:,hg_jli+2,hg_kui-0,idest) &
                   -2 *work1(:,hg_jli+2,hg_kui-1,idest) -1*work1(:,hg_jli+3,hg_kui-0,idest)

              work1(:,hg_jli-2,hg_kui+1,idest) = &
                   +20*work1(:,hg_jli+0,hg_kui-0,idest)-15*work1(:,hg_jli+0,hg_kui-1,idest) &
                   +6*work1(:,hg_jli+0,hg_kui-2,idest) -1*work1(:,hg_jli+0,hg_kui-3,idest) &
                   -30*work1(:,hg_jli+1,hg_kui-0,idest)+12*work1(:,hg_jli+1,hg_kui-1,idest) &
                   -2 *work1(:,hg_jli+1,hg_kui-2,idest)+18*work1(:,hg_jli+2,hg_kui-0,idest) &
                   -3 *work1(:,hg_jli+2,hg_kui-1,idest) -4*work1(:,hg_jli+3,hg_kui-0,idest)

              work1(:,hg_jli-2,hg_kui+2,idest) = &
                   +35*work1(:,hg_jli+0,hg_kui-0,idest)-42*work1(:,hg_jli+0,hg_kui-1,idest) &
                   +21*work1(:,hg_jli+0,hg_kui-2,idest) -4*work1(:,hg_jli+0,hg_kui-3,idest) &
                   -42*work1(:,hg_jli+1,hg_kui-0,idest)+28*work1(:,hg_jli+1,hg_kui-1,idest) &
                   -6 *work1(:,hg_jli+1,hg_kui-2,idest)+21*work1(:,hg_jli+2,hg_kui-0,idest) &
                   -6 *work1(:,hg_jli+2,hg_kui-1,idest) -4*work1(:,hg_jli+3,hg_kui-0,idest)
           endif

           if ((neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(:,hg_jui+1,hg_kui+1,idest) = &
                   +10*work1(:,hg_jui-0,hg_kui-0,idest)-10*work1(:,hg_jui-0,hg_kui-1,idest) &
                   +5 *work1(:,hg_jui-0,hg_kui-2,idest) -1*work1(:,hg_jui-0,hg_kui-3,idest) &
                   -10*work1(:,hg_jui-1,hg_kui-0,idest) +5*work1(:,hg_jui-1,hg_kui-1,idest) &
                   -1 *work1(:,hg_jui-1,hg_kui-2,idest) +5*work1(:,hg_jui-2,hg_kui-0,idest) &
                   -1 *work1(:,hg_jui-2,hg_kui-1,idest) -1*work1(:,hg_jui-3,hg_kui-0,idest)

              work1(:,hg_jui+1,hg_kui+2,idest) = &
                   +20*work1(:,hg_jui-0,hg_kui-0,idest)-30*work1(:,hg_jui-0,hg_kui-1,idest) &
                   +18*work1(:,hg_jui-0,hg_kui-2,idest) -4*work1(:,hg_jui-0,hg_kui-3,idest) &
                   -15*work1(:,hg_jui-1,hg_kui-0,idest)+12*work1(:,hg_jui-1,hg_kui-1,idest) &
                   -3 *work1(:,hg_jui-1,hg_kui-2,idest) +6*work1(:,hg_jui-2,hg_kui-0,idest) &
                   -2 *work1(:,hg_jui-2,hg_kui-1,idest) -1*work1(:,hg_jui-3,hg_kui-0,idest)

              work1(:,hg_jui+2,hg_kui+1,idest) = &
                   +20*work1(:,hg_jui-0,hg_kui-0,idest)-15*work1(:,hg_jui-0,hg_kui-1,idest) &
                   +6 *work1(:,hg_jui-0,hg_kui-2,idest) -1*work1(:,hg_jui-0,hg_kui-3,idest) &
                   -30*work1(:,hg_jui-1,hg_kui-0,idest)+12*work1(:,hg_jui-1,hg_kui-1,idest) &
                   -2 *work1(:,hg_jui-1,hg_kui-2,idest)+18*work1(:,hg_jui-2,hg_kui-0,idest) &
                   -3 *work1(:,hg_jui-2,hg_kui-1,idest) -4*work1(:,hg_jui-3,hg_kui-0,idest)

              work1(:,hg_jui+2,hg_kui+2,idest) = &
                   +35*work1(:,hg_jui-0,hg_kui-0,idest)-42*work1(:,hg_jui-0,hg_kui-1,idest) &
                   +21*work1(:,hg_jui-0,hg_kui-2,idest) -4*work1(:,hg_jui-0,hg_kui-3,idest) &
                   -42*work1(:,hg_jui-1,hg_kui-0,idest)+28*work1(:,hg_jui-1,hg_kui-1,idest) &
                   -6 *work1(:,hg_jui-1,hg_kui-2,idest)+21*work1(:,hg_jui-2,hg_kui-0,idest) &
                   -6 *work1(:,hg_jui-2,hg_kui-1,idest) -4*work1(:,hg_jui-3,hg_kui-0,idest)
           endif

           ! Should use the correct 2nd-order expressions for external corners, but for
           ! now we are just setting these cells to zero.
           !  Upper faces

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,hg_jli-1,hg_kui+1,idest) = 0.
              work1(hg_ili-2,hg_jli-1,hg_kui+1,idest) = 0.
              work1(hg_ili-1,hg_jli-2,hg_kui+1,idest) = 0.
              work1(hg_ili-2,hg_jli-2,hg_kui+1,idest) = 0.
              work1(hg_ili-1,hg_jli-1,hg_kui+2,idest) = 0.
              work1(hg_ili-2,hg_jli-1,hg_kui+2,idest) = 0.
              work1(hg_ili-1,hg_jli-2,hg_kui+2,idest) = 0.
              work1(hg_ili-2,hg_jli-2,hg_kui+2,idest) = 0.
           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_ili-1,hg_jui+1,hg_kui+1,idest) = 0.
              work1(hg_ili-2,hg_jui+1,hg_kui+1,idest) = 0.
              work1(hg_ili-1,hg_jui+2,hg_kui+1,idest) = 0.
              work1(hg_ili-2,hg_jui+2,hg_kui+1,idest) = 0.
              work1(hg_ili-1,hg_jui+1,hg_kui+2,idest) = 0.
              work1(hg_ili-2,hg_jui+1,hg_kui+2,idest) = 0.
              work1(hg_ili-1,hg_jui+2,hg_kui+2,idest) = 0.
              work1(hg_ili-2,hg_jui+2,hg_kui+2,idest) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,hg_jli-1,hg_kui+1,idest) = 0.
              work1(hg_iui+2,hg_jli-1,hg_kui+1,idest) = 0.
              work1(hg_iui+1,hg_jli-2,hg_kui+1,idest) = 0.
              work1(hg_iui+2,hg_jli-2,hg_kui+1,idest) = 0.
              work1(hg_iui+1,hg_jli-1,hg_kui+2,idest) = 0.
              work1(hg_iui+2,hg_jli-1,hg_kui+2,idest) = 0.
              work1(hg_iui+1,hg_jli-2,hg_kui+2,idest) = 0.
              work1(hg_iui+2,hg_jli-2,hg_kui+2,idest) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work1(hg_iui+1,hg_jui+1,hg_kui+1,idest) = 0.
              work1(hg_iui+2,hg_jui+1,hg_kui+1,idest) = 0.
              work1(hg_iui+1,hg_jui+2,hg_kui+1,idest) = 0.
              work1(hg_iui+2,hg_jui+2,hg_kui+1,idest) = 0.
              work1(hg_iui+1,hg_jui+1,hg_kui+2,idest) = 0.
              work1(hg_iui+2,hg_jui+1,hg_kui+2,idest) = 0.
              work1(hg_iui+1,hg_jui+2,hg_kui+2,idest) = 0.
              work1(hg_iui+2,hg_jui+2,hg_kui+2,idest) = 0.
           endif
        endif


     ! Neumann Boundary Conditions ***********************************************************************************
  else if (gr_hgBndTypes(2*NDIM-1) == MG_BND_NEUMANN) then 


        if (neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
             work1(hg_ili-1,:,:,idest) = work1(hg_ili,:,:,idest)
        if (neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
             work1(hg_iui+1,:,:,idest) = work1(hg_iui,:,:,idest)

        if (NDIM >= 2) then
           if (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work1(:,hg_jli-1,:,idest) = work1(:,hg_jli,:,idest)
           if (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work1(:,hg_jui+1,:,idest) = work1(:,hg_jui,:,idest)
        endif

        if (NDIM == 3) then
           if (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work1(:,:,hg_kli-1,idest) = work1(:,:,hg_kli,idest)
           if (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work1(:,:,hg_kui+1,idest) = work1(:,:,hg_kui,idest)
        endif


  endif

  !==============================================================================

  return
end subroutine gr_hg_amr_1blk_bcset_work

