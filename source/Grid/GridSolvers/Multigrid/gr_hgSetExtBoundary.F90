!!****if* source/Grid/GridSolvers/Multigrid/gr_hgSetExtBoundary
!!
!! NAME
!!  gr_hgSetExtBoundary
!!
!! SYNOPSIS
!! 
!!  call gr_hgSetExtBoundary(integer(IN) :: idiag, 
!!                           integer(IN) :: idir, 
!!                           logical(IN) :: extrap)
!!
!! DESCRIPTION
!!
!!  Set boundary values of the exterior guardcells of leaf blocks in the
!!  work array such that the given boundary value is satisfied.  The
!!  handled cases are MG_BND_DIRICHLET, MG_BND_NEUMANN, MG_BND_GIVENVAL
!! 
!!  Periodic boundaries are handled topologically and do not need a case
!!  here.  Allegedly this routine also handles cylindrical and spherical
!!  coordinates.
!!
!! ARGUMENTS
!!
!!  idiag  - ignored variable for diagonal guardcell case
!!  idir   - sweep all directions (0), or x, y, or z only.
!!  extrap - if this .true., then the DIRICHLET fill is assumed NOT
!!           to go to zero at the exterior face.  This is important
!!           as massive oscillations may arise for given-value problems
!!           if the value goes to zero.
!! RESULTS
!!
!!  Leaves four layers of guardcells in the work array filled.
!!
!! NOTES
!!
!! SEE ALSO
!!
!!  tot_bnd in FLASH2 / PARAMESH2 (this is a modification of it)
!!
!!***

subroutine gr_hgSetExtBoundary (idiag, idir, extrap)

  use workspace, ONLY: work
  use tree, ONLY : lrefine,lnblocks,neigh
  use gr_hgdata, ONLY: hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui, &
       gr_hgQuadrant, gr_hgGeometry, gr_hgBndTypes

#include "Multigrid.h"
#include "constants.h"
#include "Flash.h"

  implicit none

  integer, intent(IN) :: idiag, idir
  logical, intent(IN) :: extrap
  integer :: lb

  !===============================================================================


  ! Directly set any external Dirichlet or Neumann boundaries.

  ! Dirichlet or Given-value boundary conditions ***************************************************************
  if ((gr_hgBndTypes(2*NDIM-1) == MG_BND_DIRICHLET) .or. &
       (gr_hgBndTypes(2*NDIM-1) == MG_BND_GIVENVAL)) then        

     do lb = 1, lnblocks

        if ((gr_hgGeometry == MG_GEOM_1DCARTESIAN) .or. &
            (gr_hgGeometry == MG_GEOM_2DCARTESIAN) .or. &
            (gr_hgGeometry == MG_GEOM_3DCARTESIAN)) then
           ! Cartesian
           if (neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

              if (extrap) then
                 work(hg_ili-1,:,:,lb,1) = 2.*work(hg_ili,:,:,lb,1) -    work(hg_ili+1,:,:,lb,1)
                 work(hg_ili-2,:,:,lb,1) = 3.*work(hg_ili,:,:,lb,1) - 2.*work(hg_ili+1,:,:,lb,1)
                 work(hg_ili-3,:,:,lb,1) = 4.*work(hg_ili,:,:,lb,1) - 3.*work(hg_ili+1,:,:,lb,1)
                 work(hg_ili-4,:,:,lb,1) = 5.*work(hg_ili,:,:,lb,1) - 4.*work(hg_ili+1,:,:,lb,1)
              else
                 work(hg_ili-1,:,:,lb,1) = -   work(hg_ili,:,:,lb,1)
                 work(hg_ili-2,:,:,lb,1) = -3.*work(hg_ili,:,:,lb,1)
                 work(hg_ili-3,:,:,lb,1) = -5.*work(hg_ili,:,:,lb,1)
                 work(hg_ili-4,:,:,lb,1) = -7.*work(hg_ili,:,:,lb,1)
              endif

           endif

        else
           ! Non-cartesian  
           if (neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work(hg_ili-1,:,:,lb,1) = work(hg_ili,:,:,lb,1)

        endif ! end of split between cartesian and non-cartesian for ILO_FACE

        !DEV -- The upper 1D face doesn't distinguish between Cartesian and non-Cartesian.  Seems wrong. Now only Cartesian is OK in initialization.
        if (neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

           if (extrap) then
              work(hg_iui+1,:,:,lb,1) = 2.*work(hg_iui,:,:,lb,1) -    work(hg_iui-1,:,:,lb,1)
              work(hg_iui+2,:,:,lb,1) = 3.*work(hg_iui,:,:,lb,1) - 2.*work(hg_iui-1,:,:,lb,1)
              work(hg_iui+3,:,:,lb,1) = 4.*work(hg_iui,:,:,lb,1) - 3.*work(hg_iui-1,:,:,lb,1)
              work(hg_iui+4,:,:,lb,1) = 5.*work(hg_iui,:,:,lb,1) - 4.*work(hg_iui-1,:,:,lb,1)
           else
              work(hg_iui+1,:,:,lb,1) = -   work(hg_iui,:,:,lb,1)
              work(hg_iui+2,:,:,lb,1) = -3.*work(hg_iui,:,:,lb,1)
              work(hg_iui+3,:,:,lb,1) = -5.*work(hg_iui,:,:,lb,1)
              work(hg_iui+4,:,:,lb,1) = -7.*work(hg_iui,:,:,lb,1)
           endif

        endif

        ! Two-dimensional section @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (NDIM >= 2) then

           ! Begin Cartesian or polar  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
           if ((gr_hgGeometry /= MG_GEOM_2DCYLAXISYM) .or. (.not. gr_hgQuadrant)) then

              if (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then
                 if (extrap) then
                    work(:,hg_jli-1,:,lb,1) = 2.*work(:,hg_jli,:,lb,1) -    work(:,hg_jli+1,:,lb,1)
                    work(:,hg_jli-2,:,lb,1) = 3.*work(:,hg_jli,:,lb,1) - 2.*work(:,hg_jli+1,:,lb,1)
                    work(:,hg_jli-3,:,lb,1) = 4.*work(:,hg_jli,:,lb,1) - 3.*work(:,hg_jli+1,:,lb,1)
                    work(:,hg_jli-4,:,lb,1) = 5.*work(:,hg_jli,:,lb,1) - 4.*work(:,hg_jli+1,:,lb,1)
                 else
                    work(:,hg_jli-1,:,lb,1) = -   work(:,hg_jli,:,lb,1)
                    work(:,hg_jli-2,:,lb,1) = -3.*work(:,hg_jli,:,lb,1)
                    work(:,hg_jli-3,:,lb,1) = -5.*work(:,hg_jli,:,lb,1)
                    work(:,hg_jli-4,:,lb,1) = -7.*work(:,hg_jli,:,lb,1)
                 endif
              endif

              if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
                  (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
                 work(hg_ili-1,hg_jli-1,:,lb,1) =  &
                        10.*work(hg_ili,  hg_jli,  :,lb,1)   &  
                      - 10.*work(hg_ili,  hg_jli+1,:,lb,1) &  
                      +  5.*work(hg_ili,  hg_jli+2,:,lb,1) &  
                      -     work(hg_ili,  hg_jli+3,:,lb,1) &
                      - 10.*work(hg_ili+1,hg_jli,  :,lb,1) &
                      +  5.*work(hg_ili+1,hg_jli+1,:,lb,1) & 
                      -     work(hg_ili+1,hg_jli+2,:,lb,1) &
                      +  5.*work(hg_ili+2,hg_jli,  :,lb,1) &
                      -     work(hg_ili+2,hg_jli+1,:,lb,1) &
                      -     work(hg_ili+3,hg_jli,  :,lb,1)

                 work(hg_ili-1,hg_jli-2,:,lb,1) =  &
                        20.*work(hg_ili,  hg_jli,  :,lb,1) &
                      - 30.*work(hg_ili,  hg_jli+1,:,lb,1) &
                      + 18.*work(hg_ili,  hg_jli+2,:,lb,1) &
                      -  4.*work(hg_ili,  hg_jli+3,:,lb,1) &
                      - 15.*work(hg_ili+1,hg_jli,  :,lb,1) &
                      + 12.*work(hg_ili+1,hg_jli+1,:,lb,1) & 
                      -  3.*work(hg_ili+1,hg_jli+2,:,lb,1) &
                      +  6.*work(hg_ili+2,hg_jli,  :,lb,1) &
                      -  2.*work(hg_ili+2,hg_jli+1,:,lb,1) &
                      -     work(hg_ili+3,hg_jli,  :,lb,1)

                 work(hg_ili-2,hg_jli-1,:,lb,1) =         &
                        20.*work(hg_ili,  hg_jli,  :,lb,1) &
                      - 15.*work(hg_ili,  hg_jli+1,:,lb,1) &
                      +  6.*work(hg_ili,  hg_jli+2,:,lb,1) &
                      -     work(hg_ili,  hg_jli+3,:,lb,1) &
                      - 30.*work(hg_ili+1,hg_jli,  :,lb,1) &
                      + 12.*work(hg_ili+1,hg_jli+1,:,lb,1) &
                      -  2.*work(hg_ili+1,hg_jli+2,:,lb,1) &
                      + 18.*work(hg_ili+2,hg_jli,  :,lb,1) &
                      -  3.*work(hg_ili+2,hg_jli+1,:,lb,1) &
                      -  4.*work(hg_ili+3,hg_jli,  :,lb,1)

                 work(hg_ili-2,hg_jli-2,:,lb,1) =  &
                        35.*work(hg_ili,  hg_jli,  :,lb,1) &
                      - 42.*work(hg_ili,  hg_jli+1,:,lb,1) &
                      + 21.*work(hg_ili,  hg_jli+2,:,lb,1) &
                      -  4.*work(hg_ili,  hg_jli+3,:,lb,1) &
                      - 42.*work(hg_ili+1,hg_jli,  :,lb,1) &
                      + 28.*work(hg_ili+1,hg_jli+1,:,lb,1) &
                      -  6.*work(hg_ili+1,hg_jli+2,:,lb,1) &
                      + 21.*work(hg_ili+2,hg_jli,  :,lb,1) &
                      -  6.*work(hg_ili+2,hg_jli+1,:,lb,1) &
                      -  4.*work(hg_ili+3,hg_jli,  :,lb,1)
              endif

              if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
                  (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
                 work(hg_iui+1,hg_jli-1,:,lb,1) =  &
                        10.*work(hg_iui,  hg_jli,  :,lb,1) &
                      - 10.*work(hg_iui,  hg_jli+1,:,lb,1) &
                      +  5.*work(hg_iui,  hg_jli+2,:,lb,1) &
                      -     work(hg_iui,  hg_jli+3,:,lb,1) &
                      - 10.*work(hg_iui-1,hg_jli,  :,lb,1) &
                      +  5.*work(hg_iui-1,hg_jli+1,:,lb,1) &
                      -     work(hg_iui-1,hg_jli+2,:,lb,1) &
                      +  5.*work(hg_iui-2,hg_jli,  :,lb,1) &
                      -     work(hg_iui-2,hg_jli+1,:,lb,1) &
                      -     work(hg_iui-3,hg_jli,  :,lb,1)

                 work(hg_iui+1,hg_jli-2,:,lb,1) =  &
                        20.*work(hg_iui,  hg_jli,  :,lb,1) &
                      - 30.*work(hg_iui,  hg_jli+1,:,lb,1) & 
                      + 18.*work(hg_iui,  hg_jli+2,:,lb,1) &
                      -  4.*work(hg_iui,  hg_jli+3,:,lb,1) &
                      - 15.*work(hg_iui-1,hg_jli,  :,lb,1) &
                      + 12.*work(hg_iui-1,hg_jli+1,:,lb,1) &
                      -  3.*work(hg_iui-1,hg_jli+2,:,lb,1) &
                      +  6.*work(hg_iui-2,hg_jli,  :,lb,1) &
                      -  2.*work(hg_iui-2,hg_jli+1,:,lb,1) &
                      -     work(hg_iui-3,hg_jli,  :,lb,1)

                 work(hg_iui+2,hg_jli-1,:,lb,1) =  &
                        20.*work(hg_iui,  hg_jli,  :,lb,1) &
                      - 15.*work(hg_iui,  hg_jli+1,:,lb,1) &
                      +  6.*work(hg_iui,  hg_jli+2,:,lb,1) &
                      -     work(hg_iui,  hg_jli+3,:,lb,1) &
                      - 30.*work(hg_iui-1,hg_jli,  :,lb,1) &
                      + 12.*work(hg_iui-1,hg_jli+1,:,lb,1) &
                      -  2.*work(hg_iui-1,hg_jli+2,:,lb,1) &
                      + 18.*work(hg_iui-2,hg_jli,  :,lb,1) &
                      -  3.*work(hg_iui-2,hg_jli+1,:,lb,1) &
                      -  4.*work(hg_iui-3,hg_jli,  :,lb,1)

                 work(hg_iui+2,hg_jli-2,:,lb,1) =  &
                        35.*work(hg_iui,  hg_jli,  :,lb,1) &
                      - 42.*work(hg_iui,  hg_jli+1,:,lb,1) &
                      + 21.*work(hg_iui,  hg_jli+2,:,lb,1) &
                      -  4.*work(hg_iui,  hg_jli+3,:,lb,1) &
                      - 42.*work(hg_iui-1,hg_jli,  :,lb,1) &
                      + 28.*work(hg_iui-1,hg_jli+1,:,lb,1) &
                      -  6.*work(hg_iui-1,hg_jli+2,:,lb,1) &
                      + 21.*work(hg_iui-2,hg_jli,  :,lb,1) &
                      -  6.*work(hg_iui-2,hg_jli+1,:,lb,1) &
                      -  4.*work(hg_iui-3,hg_jli,  :,lb,1)
              endif


              ! End 2d cartesian or polar
           else
              ! 2d axisymmetric
              if (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                   work(:,hg_jli-1,:,lb,1) = work(:,hg_jli,:,lb,1)
           endif
           !  End of split between cartesian / polar / axisymmetric.  From here on all 2d cases are done

           if (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

              if (extrap) then
                 work(:,hg_jui+1,:,lb,1) = 2.*work(:,hg_jui,:,lb,1) -    work(:,hg_jui-1,:,lb,1)
                 work(:,hg_jui+2,:,lb,1) = 3.*work(:,hg_jui,:,lb,1) - 2.*work(:,hg_jui-1,:,lb,1)
                 work(:,hg_jui+3,:,lb,1) = 4.*work(:,hg_jui,:,lb,1) - 3.*work(:,hg_jui-1,:,lb,1)
                 work(:,hg_jui+4,:,lb,1) = 5.*work(:,hg_jui,:,lb,1) - 4.*work(:,hg_jui-1,:,lb,1)
              else
                 work(:,hg_jui+1,:,lb,1) = -   work(:,hg_jui,:,lb,1)
                 work(:,hg_jui+2,:,lb,1) = -3.*work(:,hg_jui,:,lb,1)
                 work(:,hg_jui+3,:,lb,1) = -5.*work(:,hg_jui,:,lb,1)
                 work(:,hg_jui+4,:,lb,1) = -7.*work(:,hg_jui,:,lb,1)
              endif

           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,hg_jui+1,:,lb,1) =  &
                     10.*work(hg_ili,  hg_jui,  :,lb,1) &
                   - 10.*work(hg_ili,  hg_jui-1,:,lb,1) &
                   +  5.*work(hg_ili,  hg_jui-2,:,lb,1) &
                   -     work(hg_ili,  hg_jui-3,:,lb,1) &
                   - 10.*work(hg_ili+1,hg_jui,  :,lb,1) &
                   +  5.*work(hg_ili+1,hg_jui-1,:,lb,1) &
                   -     work(hg_ili+1,hg_jui-2,:,lb,1) &
                   +  5.*work(hg_ili+2,hg_jui,  :,lb,1) &
                   -     work(hg_ili+2,hg_jui-1,:,lb,1) &
                   -     work(hg_ili+3,hg_jui,  :,lb,1)

              work(hg_ili-1,hg_jui+2,:,lb,1) =  &
                     20.*work(hg_ili,  hg_jui,  :,lb,1) &
                   - 30.*work(hg_ili,  hg_jui-1,:,lb,1) &
                   + 18.*work(hg_ili,  hg_jui-2,:,lb,1) &
                   -  4.*work(hg_ili,  hg_jui-3,:,lb,1) &
                   - 15.*work(hg_ili+1,hg_jui,  :,lb,1) &
                   + 12.*work(hg_ili+1,hg_jui-1,:,lb,1) &
                   -  3.*work(hg_ili+1,hg_jui-2,:,lb,1) &
                   +  6.*work(hg_ili+2,hg_jui,  :,lb,1) &
                   -  2.*work(hg_ili+2,hg_jui-1,:,lb,1) &
                   -     work(hg_ili+3,hg_jui,  :,lb,1)

              work(hg_ili-2,hg_jui+1,:,lb,1) =  &
                     20.*work(hg_ili,  hg_jui,  :,lb,1) &
                   - 15.*work(hg_ili,  hg_jui-1,:,lb,1) &
                   +  6.*work(hg_ili,  hg_jui-2,:,lb,1) &
                   -     work(hg_ili,  hg_jui-3,:,lb,1) &
                   - 30.*work(hg_ili+1,hg_jui,  :,lb,1) &
                   + 12.*work(hg_ili+1,hg_jui-1,:,lb,1) &
                   -  2.*work(hg_ili+1,hg_jui-2,:,lb,1) &
                   + 18.*work(hg_ili+2,hg_jui,  :,lb,1) &
                   -  3.*work(hg_ili+2,hg_jui-1,:,lb,1) &
                   -  4.*work(hg_ili+3,hg_jui,  :,lb,1)

              work(hg_ili-2,hg_jui+2,:,lb,1) =  &
                     35.*work(hg_ili,  hg_jui,  :,lb,1) & 
                   - 42.*work(hg_ili,  hg_jui-1,:,lb,1) &
                   + 21.*work(hg_ili,  hg_jui-2,:,lb,1) &
                   -  4.*work(hg_ili,  hg_jui-3,:,lb,1) &
                   - 42.*work(hg_ili+1,hg_jui,  :,lb,1) &
                   + 28.*work(hg_ili+1,hg_jui-1,:,lb,1) &
                   -  6.*work(hg_ili+1,hg_jui-2,:,lb,1) &
                   + 21.*work(hg_ili+2,hg_jui,  :,lb,1) &
                   -  6.*work(hg_ili+2,hg_jui-1,:,lb,1)&
                   -  4.*work(hg_ili+3,hg_jui,  :,lb,1)
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,hg_jui+1,:,lb,1) =  &
                     10.*work(hg_iui,  hg_jui,  :,lb,1) &
                   - 10.*work(hg_iui,  hg_jui-1,:,lb,1) &
                   +  5.*work(hg_iui,  hg_jui-2,:,lb,1) &
                   -     work(hg_iui,  hg_jui-3,:,lb,1) &
                   - 10.*work(hg_iui-1,hg_jui,  :,lb,1) &
                   +  5.*work(hg_iui-1,hg_jui-1,:,lb,1) &
                   -     work(hg_iui-1,hg_jui-2,:,lb,1) &
                   +  5.*work(hg_iui-2,hg_jui,  :,lb,1) &
                   -     work(hg_iui-2,hg_jui-1,:,lb,1) &
                   -     work(hg_iui-3,hg_jui,  :,lb,1)

              work(hg_iui+1,hg_jui+2,:,lb,1) =  &
                     20.*work(hg_iui,  hg_jui,  :,lb,1) &
                   - 30.*work(hg_iui,  hg_jui-1,:,lb,1) &
                   + 18.*work(hg_iui,  hg_jui-2,:,lb,1) &
                   -  4.*work(hg_iui,  hg_jui-3,:,lb,1) &
                   - 15.*work(hg_iui-1,hg_jui,  :,lb,1) &
                   + 12.*work(hg_iui-1,hg_jui-1,:,lb,1) &
                   -  3.*work(hg_iui-1,hg_jui-2,:,lb,1) &
                   +  6.*work(hg_iui-2,hg_jui,  :,lb,1) &
                   -  2.*work(hg_iui-2,hg_jui-1,:,lb,1) &
                   -     work(hg_iui-3,hg_jui,  :,lb,1)

              work(hg_iui+2,hg_jui+1,:,lb,1) =  &
                     20.*work(hg_iui,  hg_jui,  :,lb,1) &
                   - 15.*work(hg_iui,  hg_jui-1,:,lb,1) &
                   +  6.*work(hg_iui,  hg_jui-2,:,lb,1) &
                   -     work(hg_iui,  hg_jui-3,:,lb,1) &
                   - 30.*work(hg_iui-1,hg_jui,  :,lb,1) &
                   + 12.*work(hg_iui-1,hg_jui-1,:,lb,1) &
                   -  2.*work(hg_iui-1,hg_jui-2,:,lb,1) &
                   + 18.*work(hg_iui-2,hg_jui,  :,lb,1) &
                   -  3.*work(hg_iui-2,hg_jui-1,:,lb,1) &
                   -  4.*work(hg_iui-3,hg_jui,  :,lb,1)

              work(hg_iui+2,hg_jui+2,:,lb,1) =  &
                     35.*work(hg_iui,  hg_jui,  :,lb,1) &
                   - 42.*work(hg_iui,  hg_jui-1,:,lb,1) &
                   + 21.*work(hg_iui,  hg_jui-2,:,lb,1) &
                   -  4.*work(hg_iui,  hg_jui-3,:,lb,1) &
                   - 42.*work(hg_iui-1,hg_jui,  :,lb,1) &
                   + 28.*work(hg_iui-1,hg_jui-1,:,lb,1) &
                   -  6.*work(hg_iui-1,hg_jui-2,:,lb,1) &
                   + 21.*work(hg_iui-2,hg_jui,  :,lb,1) &
                   -  6.*work(hg_iui-2,hg_jui-1,:,lb,1) &
                   -  4.*work(hg_iui-3,hg_jui,  :,lb,1)
           endif

        endif
        ! End of two-dimensional section @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ! Start of three-dimensional section @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (NDIM == 3) then

           if (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then

              if (extrap) then
                 work(:,:,hg_kli-1,lb,1) = 2.*work(:,:,hg_kli,lb,1) -    work(:,:,hg_kli+1,lb,1)
                 work(:,:,hg_kli-2,lb,1) = 3.*work(:,:,hg_kli,lb,1) - 2.*work(:,:,hg_kli+1,lb,1)
                 work(:,:,hg_kli-3,lb,1) = 4.*work(:,:,hg_kli,lb,1) - 3.*work(:,:,hg_kli+1,lb,1)
                 work(:,:,hg_kli-4,lb,1) = 5.*work(:,:,hg_kli,lb,1) - 4.*work(:,:,hg_kli+1,lb,1)
              else
                 work(:,:,hg_kli-1,lb,1) = -   work(:,:,hg_kli,lb,1)
                 work(:,:,hg_kli-2,lb,1) = -3.*work(:,:,hg_kli,lb,1)
                 work(:,:,hg_kli-3,lb,1) = -5.*work(:,:,hg_kli,lb,1)
                 work(:,:,hg_kli-4,lb,1) = -7.*work(:,:,hg_kli,lb,1)
              endif

           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,:,hg_kli-1,lb,1) = &
                   +10*work(hg_ili+0,:,hg_kli+0,lb,1)-10*work(hg_ili+0,:,hg_kli+1,lb,1) &
                   + 5*work(hg_ili+0,:,hg_kli+2,lb,1) -1*work(hg_ili+0,:,hg_kli+3,lb,1) &
                   -10*work(hg_ili+1,:,hg_kli+0,lb,1) +5*work(hg_ili+1,:,hg_kli+1,lb,1) &
                   -1 *work(hg_ili+1,:,hg_kli+2,lb,1) +5*work(hg_ili+2,:,hg_kli+0,lb,1) &
                   -1 *work(hg_ili+2,:,hg_kli+1,lb,1) -1*work(hg_ili+3,:,hg_kli+0,lb,1)

              work(hg_ili-1,:,hg_kli-2,lb,1) = &
                   +20*work(hg_ili+0,:,hg_kli+0,lb,1)-30*work(hg_ili+0,:,hg_kli+1,lb,1) &
                   +18*work(hg_ili+0,:,hg_kli+2,lb,1) -4*work(hg_ili+0,:,hg_kli+3,lb,1) &
                   -15*work(hg_ili+1,:,hg_kli+0,lb,1)+12*work(hg_ili+1,:,hg_kli+1,lb,1) &
                   -3 *work(hg_ili+1,:,hg_kli+2,lb,1) +6*work(hg_ili+2,:,hg_kli+0,lb,1) &
                   -2 *work(hg_ili+2,:,hg_kli+1,lb,1) -1*work(hg_ili+3,:,hg_kli+0,lb,1)

              work(hg_ili-2,:,hg_kli-1,lb,1) = &
                   +20*work(hg_ili+0,:,hg_kli+0,lb,1)-15*work(hg_ili+0,:,hg_kli+1,lb,1) &
                   +6 *work(hg_ili+0,:,hg_kli+2,lb,1) -1*work(hg_ili+0,:,hg_kli+3,lb,1) &
                   -30*work(hg_ili+1,:,hg_kli+0,lb,1)+12*work(hg_ili+1,:,hg_kli+1,lb,1) &
                   -2 *work(hg_ili+1,:,hg_kli+2,lb,1)+18*work(hg_ili+2,:,hg_kli+0,lb,1) &
                   -3 *work(hg_ili+2,:,hg_kli+1,lb,1) -4*work(hg_ili+3,:,hg_kli+0,lb,1)

              work(hg_ili-2,:,hg_kli-2,lb,1) = &
                   +35*work(hg_ili+0,:,hg_kli+0,lb,1)-42*work(hg_ili+0,:,hg_kli+1,lb,1) &
                   +21*work(hg_ili+0,:,hg_kli+2,lb,1) -4*work(hg_ili+0,:,hg_kli+3,lb,1) &
                   -42*work(hg_ili+1,:,hg_kli+0,lb,1)+28*work(hg_ili+1,:,hg_kli+1,lb,1) &
                   -6 *work(hg_ili+1,:,hg_kli+2,lb,1)+21*work(hg_ili+2,:,hg_kli+0,lb,1) &
                   -6 *work(hg_ili+2,:,hg_kli+1,lb,1) -4*work(hg_ili+3,:,hg_kli+0,lb,1)
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,:,hg_kli-1,lb,1) = &
                   +10*work(hg_iui-0,:,hg_kli+0,lb,1)-10*work(hg_iui-0,:,hg_kli+1,lb,1) &
                   +5 *work(hg_iui-0,:,hg_kli+2,lb,1) -1*work(hg_iui-0,:,hg_kli+3,lb,1) &
                   -10*work(hg_iui-1,:,hg_kli+0,lb,1) +5*work(hg_iui-1,:,hg_kli+1,lb,1) &
                   -1 *work(hg_iui-1,:,hg_kli+2,lb,1) +5*work(hg_iui-2,:,hg_kli+0,lb,1) &
                   -1 *work(hg_iui-2,:,hg_kli+1,lb,1) -1*work(hg_iui-3,:,hg_kli+0,lb,1)

              work(hg_iui+1,:,hg_kli-2,lb,1) = &
                   +20*work(hg_iui-0,:,hg_kli+0,lb,1)-30*work(hg_iui-0,:,hg_kli+1,lb,1) &
                   +18*work(hg_iui-0,:,hg_kli+2,lb,1) -4*work(hg_iui-0,:,hg_kli+3,lb,1) &
                   -15*work(hg_iui-1,:,hg_kli+0,lb,1)+12*work(hg_iui-1,:,hg_kli+1,lb,1) &
                   -3 *work(hg_iui-1,:,hg_kli+2,lb,1) +6*work(hg_iui-2,:,hg_kli+0,lb,1) &
                   -2 *work(hg_iui-2,:,hg_kli+1,lb,1) -1*work(hg_iui-3,:,hg_kli+0,lb,1)

              work(hg_iui+2,:,hg_kli-1,lb,1) = &
                   +20*work(hg_iui-0,:,hg_kli+0,lb,1)-15*work(hg_iui-0,:,hg_kli+1,lb,1) &
                   +6 *work(hg_iui-0,:,hg_kli+2,lb,1) -1*work(hg_iui-0,:,hg_kli+3,lb,1) &
                   -30*work(hg_iui-1,:,hg_kli+0,lb,1)+12*work(hg_iui-1,:,hg_kli+1,lb,1) &
                   -2 *work(hg_iui-1,:,hg_kli+2,lb,1)+18*work(hg_iui-2,:,hg_kli+0,lb,1) &
                   -3 *work(hg_iui-2,:,hg_kli+1,lb,1) -4*work(hg_iui-3,:,hg_kli+0,lb,1)

              work(hg_iui+2,:,hg_kli-2,lb,1) = &
                   +35*work(hg_iui-0,:,hg_kli+0,lb,1)-42*work(hg_iui-0,:,hg_kli+1,lb,1) &
                   +21*work(hg_iui-0,:,hg_kli+2,lb,1) -4*work(hg_iui-0,:,hg_kli+3,lb,1) &
                   -42*work(hg_iui-1,:,hg_kli+0,lb,1)+28*work(hg_iui-1,:,hg_kli+1,lb,1) &
                   -6 *work(hg_iui-1,:,hg_kli+2,lb,1)+21*work(hg_iui-2,:,hg_kli+0,lb,1) &
                   -6 *work(hg_iui-2,:,hg_kli+1,lb,1) -4*work(hg_iui-3,:,hg_kli+0,lb,1)
           endif

           if ((neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(:,hg_jli-1,hg_kli-1,lb,1) = &
                   +10*work(:,hg_jli+0,hg_kli+0,lb,1)-10*work(:,hg_jli+0,hg_kli+1,lb,1) &
                   +5 *work(:,hg_jli+0,hg_kli+2,lb,1) -1*work(:,hg_jli+0,hg_kli+3,lb,1) &
                   -10*work(:,hg_jli+1,hg_kli+0,lb,1) +5*work(:,hg_jli+1,hg_kli+1,lb,1) &
                   -1 *work(:,hg_jli+1,hg_kli+2,lb,1) +5*work(:,hg_jli+2,hg_kli+0,lb,1) &
                   -1 *work(:,hg_jli+2,hg_kli+1,lb,1) -1*work(:,hg_jli+3,hg_kli+0,lb,1)

              work(:,hg_jli-1,hg_kli-2,lb,1) = &
                   +20*work(:,hg_jli+0,hg_kli+0,lb,1)-30*work(:,hg_jli+0,hg_kli+1,lb,1) &
                   +18*work(:,hg_jli+0,hg_kli+2,lb,1) -4*work(:,hg_jli+0,hg_kli+3,lb,1) &
                   -15*work(:,hg_jli+1,hg_kli+0,lb,1)+12*work(:,hg_jli+1,hg_kli+1,lb,1) &
                   -3 *work(:,hg_jli+1,hg_kli+2,lb,1) +6*work(:,hg_jli+2,hg_kli+0,lb,1) &
                   -2 *work(:,hg_jli+2,hg_kli+1,lb,1) -1*work(:,hg_jli+3,hg_kli+0,lb,1)

              work(:,hg_jli-2,hg_kli-1,lb,1) = &
                   +20*work(:,hg_jli+0,hg_kli+0,lb,1)-15*work(:,hg_jli+0,hg_kli+1,lb,1) &
                   +6 *work(:,hg_jli+0,hg_kli+2,lb,1) -1*work(:,hg_jli+0,hg_kli+3,lb,1) &
                   -30*work(:,hg_jli+1,hg_kli+0,lb,1)+12*work(:,hg_jli+1,hg_kli+1,lb,1) &
                   -2 *work(:,hg_jli+1,hg_kli+2,lb,1)+18*work(:,hg_jli+2,hg_kli+0,lb,1) &
                   -3 *work(:,hg_jli+2,hg_kli+1,lb,1) -4*work(:,hg_jli+3,hg_kli+0,lb,1)

              work(:,hg_jli-2,hg_kli-2,lb,1) = &
                   +35*work(:,hg_jli+0,hg_kli+0,lb,1)-42*work(:,hg_jli+0,hg_kli+1,lb,1) &
                   +21*work(:,hg_jli+0,hg_kli+2,lb,1) -4*work(:,hg_jli+0,hg_kli+3,lb,1) &
                   -42*work(:,hg_jli+1,hg_kli+0,lb,1)+28*work(:,hg_jli+1,hg_kli+1,lb,1) &
                   -6 *work(:,hg_jli+1,hg_kli+2,lb,1)+21*work(:,hg_jli+2,hg_kli+0,lb,1) &
                   -6 *work(:,hg_jli+2,hg_kli+1,lb,1) -4*work(:,hg_jli+3,hg_kli+0,lb,1)
           endif

           if ((neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(:,hg_jui+1,hg_kli-1,lb,1) = &
                   +10*work(:,hg_jui-0,hg_kli+0,lb,1)-10*work(:,hg_jui-0,hg_kli+1,lb,1) &
                   + 5*work(:,hg_jui-0,hg_kli+2,lb,1) -1*work(:,hg_jui-0,hg_kli+3,lb,1) &
                   -10*work(:,hg_jui-1,hg_kli+0,lb,1) +5*work(:,hg_jui-1,hg_kli+1,lb,1) &
                   -1 *work(:,hg_jui-1,hg_kli+2,lb,1) +5*work(:,hg_jui-2,hg_kli+0,lb,1) &
                   -1 *work(:,hg_jui-2,hg_kli+1,lb,1) -1*work(:,hg_jui-3,hg_kli+0,lb,1)

              work(:,hg_jui+1,hg_kli-2,lb,1) = &
                   +20*work(:,hg_jui-0,hg_kli+0,lb,1)-30*work(:,hg_jui-0,hg_kli+1,lb,1) &
                   +18*work(:,hg_jui-0,hg_kli+2,lb,1) -4*work(:,hg_jui-0,hg_kli+3,lb,1) &
                   -15*work(:,hg_jui-1,hg_kli+0,lb,1)+12*work(:,hg_jui-1,hg_kli+1,lb,1) &
                   -3 *work(:,hg_jui-1,hg_kli+2,lb,1) +6*work(:,hg_jui-2,hg_kli+0,lb,1) &
                   -2 *work(:,hg_jui-2,hg_kli+1,lb,1) -1*work(:,hg_jui-3,hg_kli+0,lb,1)

              work(:,hg_jui+2,hg_kli-1,lb,1) = &
                   +20*work(:,hg_jui-0,hg_kli+0,lb,1)-15*work(:,hg_jui-0,hg_kli+1,lb,1) &
                   +6 *work(:,hg_jui-0,hg_kli+2,lb,1) -1*work(:,hg_jui-0,hg_kli+3,lb,1) &
                   -30*work(:,hg_jui-1,hg_kli+0,lb,1)+12*work(:,hg_jui-1,hg_kli+1,lb,1) &
                   -2 *work(:,hg_jui-1,hg_kli+2,lb,1)+18*work(:,hg_jui-2,hg_kli+0,lb,1) &
                   -3 *work(:,hg_jui-2,hg_kli+1,lb,1) -4*work(:,hg_jui-3,hg_kli+0,lb,1)

              work(:,hg_jui+2,hg_kli-2,lb,1) = &
                   +35*work(:,hg_jui-0,hg_kli+0,lb,1)-42*work(:,hg_jui-0,hg_kli+1,lb,1) &
                   +21*work(:,hg_jui-0,hg_kli+2,lb,1) -4*work(:,hg_jui-0,hg_kli+3,lb,1) &
                   -42*work(:,hg_jui-1,hg_kli+0,lb,1)+28*work(:,hg_jui-1,hg_kli+1,lb,1) &
                   -6 *work(:,hg_jui-1,hg_kli+2,lb,1)+21*work(:,hg_jui-2,hg_kli+0,lb,1) &
                   -6 *work(:,hg_jui-2,hg_kli+1,lb,1) -4*work(:,hg_jui-3,hg_kli+0,lb,1)
           endif

           ! Should use the correct 2nd-order expressions for external corners, but for
           ! now we are just setting these cells to zero.
           !  Lower faces of K, high ones are below
           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,hg_jli-1,hg_kli-1,lb,1) = 0.
              work(hg_ili-2,hg_jli-1,hg_kli-1,lb,1) = 0.
              work(hg_ili-1,hg_jli-2,hg_kli-1,lb,1) = 0.
              work(hg_ili-2,hg_jli-2,hg_kli-1,lb,1) = 0.
              work(hg_ili-1,hg_jli-1,hg_kli-2,lb,1) = 0.
              work(hg_ili-2,hg_jli-1,hg_kli-2,lb,1) = 0.
              work(hg_ili-1,hg_jli-2,hg_kli-2,lb,1) = 0.
              work(hg_ili-2,hg_jli-2,hg_kli-2,lb,1) = 0.
           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,hg_jui+1,hg_kli-1,lb,1) = 0.
              work(hg_ili-2,hg_jui+1,hg_kli-1,lb,1) = 0.
              work(hg_ili-1,hg_jui+2,hg_kli-1,lb,1) = 0.
              work(hg_ili-2,hg_jui+2,hg_kli-1,lb,1) = 0.
              work(hg_ili-1,hg_jui+1,hg_kli-2,lb,1) = 0.
              work(hg_ili-2,hg_jui+1,hg_kli-2,lb,1) = 0.
              work(hg_ili-1,hg_jui+2,hg_kli-2,lb,1) = 0.
              work(hg_ili-2,hg_jui+2,hg_kli-2,lb,1) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,hg_jli-1,hg_kli-1,lb,1) = 0.
              work(hg_iui+2,hg_jli-1,hg_kli-1,lb,1) = 0.
              work(hg_iui+1,hg_jli-2,hg_kli-1,lb,1) = 0.
              work(hg_iui+2,hg_jli-2,hg_kli-1,lb,1) = 0.
              work(hg_iui+1,hg_jli-1,hg_kli-2,lb,1) = 0.
              work(hg_iui+2,hg_jli-1,hg_kli-2,lb,1) = 0.
              work(hg_iui+1,hg_jli-2,hg_kli-2,lb,1) = 0.
              work(hg_iui+2,hg_jli-2,hg_kli-2,lb,1) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,hg_jui+1,hg_kli-1,lb,1) = 0.
              work(hg_iui+2,hg_jui+1,hg_kli-1,lb,1) = 0.
              work(hg_iui+1,hg_jui+2,hg_kli-1,lb,1) = 0.
              work(hg_iui+2,hg_jui+2,hg_kli-1,lb,1) = 0.
              work(hg_iui+1,hg_jui+1,hg_kli-2,lb,1) = 0.
              work(hg_iui+2,hg_jui+1,hg_kli-2,lb,1) = 0.
              work(hg_iui+1,hg_jui+2,hg_kli-2,lb,1) = 0.
              work(hg_iui+2,hg_jui+2,hg_kli-2,lb,1) = 0.
           endif


           if (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) then
              if (extrap) then
                 work(:,:,hg_kui+1,lb,1) = 2.*work(:,:,hg_kui,lb,1) -    work(:,:,hg_kui-1,lb,1)
                 work(:,:,hg_kui+2,lb,1) = 3.*work(:,:,hg_kui,lb,1) - 2.*work(:,:,hg_kui-1,lb,1)
                 work(:,:,hg_kui+3,lb,1) = 4.*work(:,:,hg_kui,lb,1) - 3.*work(:,:,hg_kui-1,lb,1)
                 work(:,:,hg_kui+4,lb,1) = 5.*work(:,:,hg_kui,lb,1) - 4.*work(:,:,hg_kui-1,lb,1)
              else
                 work(:,:,hg_kui+1,lb,1) = -   work(:,:,hg_kui,lb,1)
                 work(:,:,hg_kui+2,lb,1) = -3.*work(:,:,hg_kui,lb,1)
                 work(:,:,hg_kui+3,lb,1) = -5.*work(:,:,hg_kui,lb,1)
                 work(:,:,hg_kui+4,lb,1) = -7.*work(:,:,hg_kui,lb,1)
              endif
           endif

           ! DEV could nest these if statements within the previous one, minor savings in execution
           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,:,hg_kui+1,lb,1) = &
                   +10*work(hg_ili+0,:,hg_kui-0,lb,1)-10*work(hg_ili+0,:,hg_kui-1,lb,1) &
                   +5 *work(hg_ili+0,:,hg_kui-2,lb,1) -1*work(hg_ili+0,:,hg_kui-3,lb,1) &
                   -10*work(hg_ili+1,:,hg_kui-0,lb,1) +5*work(hg_ili+1,:,hg_kui-1,lb,1) &
                   -1 *work(hg_ili+1,:,hg_kui-2,lb,1) +5*work(hg_ili+2,:,hg_kui-0,lb,1) &
                   -1 *work(hg_ili+2,:,hg_kui-1,lb,1) -1*work(hg_ili+3,:,hg_kui-0,lb,1)

              work(hg_ili-1,:,hg_kui+2,lb,1) = &
                   +20*work(hg_ili+0,:,hg_kui-0,lb,1)-30*work(hg_ili+0,:,hg_kui-1,lb,1) &
                   +18*work(hg_ili+0,:,hg_kui-2,lb,1) -4*work(hg_ili+0,:,hg_kui-3,lb,1) &
                   -15*work(hg_ili+1,:,hg_kui-0,lb,1)+12*work(hg_ili+1,:,hg_kui-1,lb,1) &
                   -3 *work(hg_ili+1,:,hg_kui-2,lb,1) +6*work(hg_ili+2,:,hg_kui-0,lb,1) &
                   -2 *work(hg_ili+2,:,hg_kui-1,lb,1) -1*work(hg_ili+3,:,hg_kui-0,lb,1)

              work(hg_ili-2,:,hg_kui+1,lb,1) = &
                   +20*work(hg_ili+0,:,hg_kui-0,lb,1)-15*work(hg_ili+0,:,hg_kui-1,lb,1) &
                   +6 *work(hg_ili+0,:,hg_kui-2,lb,1) -1*work(hg_ili+0,:,hg_kui-3,lb,1) &
                   -30*work(hg_ili+1,:,hg_kui-0,lb,1)+12*work(hg_ili+1,:,hg_kui-1,lb,1) &
                   -2 *work(hg_ili+1,:,hg_kui-2,lb,1)+18*work(hg_ili+2,:,hg_kui-0,lb,1) &
                   -3 *work(hg_ili+2,:,hg_kui-1,lb,1) -4*work(hg_ili+3,:,hg_kui-0,lb,1)

              work(hg_ili-2,:,hg_kui+2,lb,1) = &
                   +35*work(hg_ili+0,:,hg_kui-0,lb,1)-42*work(hg_ili+0,:,hg_kui-1,lb,1) &
                   +21*work(hg_ili+0,:,hg_kui-2,lb,1) -4*work(hg_ili+0,:,hg_kui-3,lb,1) &
                   -42*work(hg_ili+1,:,hg_kui-0,lb,1)+28*work(hg_ili+1,:,hg_kui-1,lb,1) &
                   -6 *work(hg_ili+1,:,hg_kui-2,lb,1)+21*work(hg_ili+2,:,hg_kui-0,lb,1) &
                   -6 *work(hg_ili+2,:,hg_kui-1,lb,1) -4*work(hg_ili+3,:,hg_kui-0,lb,1)
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,:,hg_kui+1,lb,1) = &
                   +10*work(hg_iui-0,:,hg_kui-0,lb,1)-10*work(hg_iui-0,:,hg_kui-1,lb,1) &
                   +5 *work(hg_iui-0,:,hg_kui-2,lb,1) -1*work(hg_iui-0,:,hg_kui-3,lb,1) &
                   -10*work(hg_iui-1,:,hg_kui-0,lb,1) +5*work(hg_iui-1,:,hg_kui-1,lb,1) &
                   -1 *work(hg_iui-1,:,hg_kui-2,lb,1) +5*work(hg_iui-2,:,hg_kui-0,lb,1) &
                   -1 *work(hg_iui-2,:,hg_kui-1,lb,1) -1*work(hg_iui-3,:,hg_kui-0,lb,1)

              work(hg_iui+1,:,hg_kui+2,lb,1) = &
                   +20*work(hg_iui-0,:,hg_kui-0,lb,1)-30*work(hg_iui-0,:,hg_kui-1,lb,1) &
                   +18*work(hg_iui-0,:,hg_kui-2,lb,1) -4*work(hg_iui-0,:,hg_kui-3,lb,1) &
                   -15*work(hg_iui-1,:,hg_kui-0,lb,1)+12*work(hg_iui-1,:,hg_kui-1,lb,1) &
                   -3 *work(hg_iui-1,:,hg_kui-2,lb,1) +6*work(hg_iui-2,:,hg_kui-0,lb,1) &
                   -2 *work(hg_iui-2,:,hg_kui-1,lb,1) -1*work(hg_iui-3,:,hg_kui-0,lb,1)

              work(hg_iui+2,:,hg_kui+1,lb,1) = &
                   +20*work(hg_iui-0,:,hg_kui-0,lb,1)-15*work(hg_iui-0,:,hg_kui-1,lb,1) &
                   +6 *work(hg_iui-0,:,hg_kui-2,lb,1) -1*work(hg_iui-0,:,hg_kui-3,lb,1) &
                   -30*work(hg_iui-1,:,hg_kui-0,lb,1)+12*work(hg_iui-1,:,hg_kui-1,lb,1) &
                   -2 *work(hg_iui-1,:,hg_kui-2,lb,1)+18*work(hg_iui-2,:,hg_kui-0,lb,1) &
                   -3 *work(hg_iui-2,:,hg_kui-1,lb,1) -4*work(hg_iui-3,:,hg_kui-0,lb,1)

              work(hg_iui+2,:,hg_kui+2,lb,1) = &
                   +35*work(hg_iui-0,:,hg_kui-0,lb,1)-42*work(hg_iui-0,:,hg_kui-1,lb,1) &
                   +21*work(hg_iui-0,:,hg_kui-2,lb,1) -4*work(hg_iui-0,:,hg_kui-3,lb,1) &
                   -42*work(hg_iui-1,:,hg_kui-0,lb,1)+28*work(hg_iui-1,:,hg_kui-1,lb,1) &
                   -6 *work(hg_iui-1,:,hg_kui-2,lb,1)+21*work(hg_iui-2,:,hg_kui-0,lb,1) &
                   -6 *work(hg_iui-2,:,hg_kui-1,lb,1) -4*work(hg_iui-3,:,hg_kui-0,lb,1)
           endif

           if ((neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(:,hg_jli-1,hg_kui+1,lb,1) = &
                   +10*work(:,hg_jli+0,hg_kui-0,lb,1)-10*work(:,hg_jli+0,hg_kui-1,lb,1) &
                   +5 *work(:,hg_jli+0,hg_kui-2,lb,1) -1*work(:,hg_jli+0,hg_kui-3,lb,1) &
                   -10*work(:,hg_jli+1,hg_kui-0,lb,1) +5*work(:,hg_jli+1,hg_kui-1,lb,1) &
                   -1 *work(:,hg_jli+1,hg_kui-2,lb,1) +5*work(:,hg_jli+2,hg_kui-0,lb,1) &
                   -1 *work(:,hg_jli+2,hg_kui-1,lb,1) -1*work(:,hg_jli+3,hg_kui-0,lb,1)

              work(:,hg_jli-1,hg_kui+2,lb,1) = &
                   +20*work(:,hg_jli+0,hg_kui-0,lb,1)-30*work(:,hg_jli+0,hg_kui-1,lb,1) &
                   +18*work(:,hg_jli+0,hg_kui-2,lb,1) -4*work(:,hg_jli+0,hg_kui-3,lb,1) &
                   -15*work(:,hg_jli+1,hg_kui-0,lb,1)+12*work(:,hg_jli+1,hg_kui-1,lb,1) &
                   -3 *work(:,hg_jli+1,hg_kui-2,lb,1) +6*work(:,hg_jli+2,hg_kui-0,lb,1) &
                   -2 *work(:,hg_jli+2,hg_kui-1,lb,1) -1*work(:,hg_jli+3,hg_kui-0,lb,1)

              work(:,hg_jli-2,hg_kui+1,lb,1) = &
                   +20*work(:,hg_jli+0,hg_kui-0,lb,1)-15*work(:,hg_jli+0,hg_kui-1,lb,1) &
                   +6*work(:,hg_jli+0,hg_kui-2,lb,1) -1*work(:,hg_jli+0,hg_kui-3,lb,1) &
                   -30*work(:,hg_jli+1,hg_kui-0,lb,1)+12*work(:,hg_jli+1,hg_kui-1,lb,1) &
                   -2 *work(:,hg_jli+1,hg_kui-2,lb,1)+18*work(:,hg_jli+2,hg_kui-0,lb,1) &
                   -3 *work(:,hg_jli+2,hg_kui-1,lb,1) -4*work(:,hg_jli+3,hg_kui-0,lb,1)

              work(:,hg_jli-2,hg_kui+2,lb,1) = &
                   +35*work(:,hg_jli+0,hg_kui-0,lb,1)-42*work(:,hg_jli+0,hg_kui-1,lb,1) &
                   +21*work(:,hg_jli+0,hg_kui-2,lb,1) -4*work(:,hg_jli+0,hg_kui-3,lb,1) &
                   -42*work(:,hg_jli+1,hg_kui-0,lb,1)+28*work(:,hg_jli+1,hg_kui-1,lb,1) &
                   -6 *work(:,hg_jli+1,hg_kui-2,lb,1)+21*work(:,hg_jli+2,hg_kui-0,lb,1) &
                   -6 *work(:,hg_jli+2,hg_kui-1,lb,1) -4*work(:,hg_jli+3,hg_kui-0,lb,1)
           endif

           if ((neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(:,hg_jui+1,hg_kui+1,lb,1) = &
                   +10*work(:,hg_jui-0,hg_kui-0,lb,1)-10*work(:,hg_jui-0,hg_kui-1,lb,1) &
                   +5 *work(:,hg_jui-0,hg_kui-2,lb,1) -1*work(:,hg_jui-0,hg_kui-3,lb,1) &
                   -10*work(:,hg_jui-1,hg_kui-0,lb,1) +5*work(:,hg_jui-1,hg_kui-1,lb,1) &
                   -1 *work(:,hg_jui-1,hg_kui-2,lb,1) +5*work(:,hg_jui-2,hg_kui-0,lb,1) &
                   -1 *work(:,hg_jui-2,hg_kui-1,lb,1) -1*work(:,hg_jui-3,hg_kui-0,lb,1)

              work(:,hg_jui+1,hg_kui+2,lb,1) = &
                   +20*work(:,hg_jui-0,hg_kui-0,lb,1)-30*work(:,hg_jui-0,hg_kui-1,lb,1) &
                   +18*work(:,hg_jui-0,hg_kui-2,lb,1) -4*work(:,hg_jui-0,hg_kui-3,lb,1) &
                   -15*work(:,hg_jui-1,hg_kui-0,lb,1)+12*work(:,hg_jui-1,hg_kui-1,lb,1) &
                   -3 *work(:,hg_jui-1,hg_kui-2,lb,1) +6*work(:,hg_jui-2,hg_kui-0,lb,1) &
                   -2 *work(:,hg_jui-2,hg_kui-1,lb,1) -1*work(:,hg_jui-3,hg_kui-0,lb,1)

              work(:,hg_jui+2,hg_kui+1,lb,1) = &
                   +20*work(:,hg_jui-0,hg_kui-0,lb,1)-15*work(:,hg_jui-0,hg_kui-1,lb,1) &
                   +6 *work(:,hg_jui-0,hg_kui-2,lb,1) -1*work(:,hg_jui-0,hg_kui-3,lb,1) &
                   -30*work(:,hg_jui-1,hg_kui-0,lb,1)+12*work(:,hg_jui-1,hg_kui-1,lb,1) &
                   -2 *work(:,hg_jui-1,hg_kui-2,lb,1)+18*work(:,hg_jui-2,hg_kui-0,lb,1) &
                   -3 *work(:,hg_jui-2,hg_kui-1,lb,1) -4*work(:,hg_jui-3,hg_kui-0,lb,1)

              work(:,hg_jui+2,hg_kui+2,lb,1) = &
                   +35*work(:,hg_jui-0,hg_kui-0,lb,1)-42*work(:,hg_jui-0,hg_kui-1,lb,1) &
                   +21*work(:,hg_jui-0,hg_kui-2,lb,1) -4*work(:,hg_jui-0,hg_kui-3,lb,1) &
                   -42*work(:,hg_jui-1,hg_kui-0,lb,1)+28*work(:,hg_jui-1,hg_kui-1,lb,1) &
                   -6 *work(:,hg_jui-1,hg_kui-2,lb,1)+21*work(:,hg_jui-2,hg_kui-0,lb,1) &
                   -6 *work(:,hg_jui-2,hg_kui-1,lb,1) -4*work(:,hg_jui-3,hg_kui-0,lb,1)
           endif

           ! Should use the correct 2nd-order expressions for external corners, but for
           ! now we are just setting these cells to zero.
           !  Upper faces

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,hg_jli-1,hg_kui+1,lb,1) = 0.
              work(hg_ili-2,hg_jli-1,hg_kui+1,lb,1) = 0.
              work(hg_ili-1,hg_jli-2,hg_kui+1,lb,1) = 0.
              work(hg_ili-2,hg_jli-2,hg_kui+1,lb,1) = 0.
              work(hg_ili-1,hg_jli-1,hg_kui+2,lb,1) = 0.
              work(hg_ili-2,hg_jli-1,hg_kui+2,lb,1) = 0.
              work(hg_ili-1,hg_jli-2,hg_kui+2,lb,1) = 0.
              work(hg_ili-2,hg_jli-2,hg_kui+2,lb,1) = 0.
           endif

           if ((neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_ili-1,hg_jui+1,hg_kui+1,lb,1) = 0.
              work(hg_ili-2,hg_jui+1,hg_kui+1,lb,1) = 0.
              work(hg_ili-1,hg_jui+2,hg_kui+1,lb,1) = 0.
              work(hg_ili-2,hg_jui+2,hg_kui+1,lb,1) = 0.
              work(hg_ili-1,hg_jui+1,hg_kui+2,lb,1) = 0.
              work(hg_ili-2,hg_jui+1,hg_kui+2,lb,1) = 0.
              work(hg_ili-1,hg_jui+2,hg_kui+2,lb,1) = 0.
              work(hg_ili-2,hg_jui+2,hg_kui+2,lb,1) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,hg_jli-1,hg_kui+1,lb,1) = 0.
              work(hg_iui+2,hg_jli-1,hg_kui+1,lb,1) = 0.
              work(hg_iui+1,hg_jli-2,hg_kui+1,lb,1) = 0.
              work(hg_iui+2,hg_jli-2,hg_kui+1,lb,1) = 0.
              work(hg_iui+1,hg_jli-1,hg_kui+2,lb,1) = 0.
              work(hg_iui+2,hg_jli-1,hg_kui+2,lb,1) = 0.
              work(hg_iui+1,hg_jli-2,hg_kui+2,lb,1) = 0.
              work(hg_iui+2,hg_jli-2,hg_kui+2,lb,1) = 0.
           endif

           if ((neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) .and. &
               (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY)) then
              work(hg_iui+1,hg_jui+1,hg_kui+1,lb,1) = 0.
              work(hg_iui+2,hg_jui+1,hg_kui+1,lb,1) = 0.
              work(hg_iui+1,hg_jui+2,hg_kui+1,lb,1) = 0.
              work(hg_iui+2,hg_jui+2,hg_kui+1,lb,1) = 0.
              work(hg_iui+1,hg_jui+1,hg_kui+2,lb,1) = 0.
              work(hg_iui+2,hg_jui+1,hg_kui+2,lb,1) = 0.
              work(hg_iui+1,hg_jui+2,hg_kui+2,lb,1) = 0.
              work(hg_iui+2,hg_jui+2,hg_kui+2,lb,1) = 0.
           endif
        endif
     enddo


     ! Neumann Boundary Conditions ***********************************************************************************
  else if (gr_hgBndTypes(2*NDIM-1) == MG_BND_NEUMANN) then 

     do lb = 1, lnblocks

        if (neigh(1,ILO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
             work(hg_ili-1,:,:,lb,1) = work(hg_ili,:,:,lb,1)
        if (neigh(1,IHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
             work(hg_iui+1,:,:,lb,1) = work(hg_iui,:,:,lb,1)

        if (NDIM >= 2) then
           if (neigh(1,JLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work(:,hg_jli-1,:,lb,1) = work(:,hg_jli,:,lb,1)
           if (neigh(1,JHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work(:,hg_jui+1,:,lb,1) = work(:,hg_jui,:,lb,1)
        endif

        if (NDIM == 3) then
           if (neigh(1,KLO_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work(:,:,hg_kli-1,lb,1) = work(:,:,hg_kli,lb,1)
           if (neigh(1,KHI_FACE,lb) <= PARAMESH_PHYSICAL_BOUNDARY) &
                work(:,:,hg_kui+1,lb,1) = work(:,:,hg_kui,lb,1)
        endif

     enddo

  endif

  !==============================================================================

  return
end subroutine gr_hgSetExtBoundary

