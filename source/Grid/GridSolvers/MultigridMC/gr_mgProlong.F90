!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgProlong
!!
!! NAME
!!
!!  gr_mgProlong
!!
!! SYNOPSIS
!!
!!  call gr_mgProlong(integer(in) :: level,
!!                    integer(in) :: ifrom,
!!                    integer(in) :: ito,
!!                    integer(in) :: add)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   ifrom : 
!!
!!   ito : 
!!
!!   add : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_prolong()

!  Description: Prolongate a multigrid variable from one mesh level to the
!               next finer level.  All blocks on the "from" level are assumed
!               to contain valid data.  On output, all blocks on the "to"
!               level should contain valid data.

!  Parameters:  level       Level to prolongate from.
!               ifrom       Index of the variable to prolongate from.
!               ito         Index of the variable to prolongate into.
!               add         If /= 0, add to contents of ito; otherwise replace.


      subroutine gr_mgProlong (level, ifrom, ito, add)

!===============================================================================

      use gr_mgData

      use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

      use tree, only : nodetype,newchild,lrefine,child

      use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                    Grid_getBlkPtr,       &
                                    Grid_releaseBlkPtr

      use workspace, ONLY: work

      implicit none
#include "constants.h"
#include "Flash.h"
      integer, intent(in) :: level, ifrom, ito, add

      integer :: lb

      real    :: normf, normt

      integer, parameter :: mchild2 = 8
      integer, dimension(mchild2) :: child2

      integer, save :: numPEs2, myPE2
      integer :: lnblocks2

      integer, pointer, dimension(:), save      :: nodetype2
      logical, pointer, dimension(:), save      :: newchild2
      real, pointer, dimension(:,:,:,:), save :: unkt
      logical, save :: first_call = .true.

      integer, parameter :: k2d = K2D
      integer, parameter :: k3d = K3D

      integer :: iw

!===============================================================================

      if (first_call) then
         nodetype2 => nodetype 
         newchild2 => newchild 
         myPE2   = gr_meshMe 
         numPEs2 = gr_meshNumProcs
      end if
      
!               Get data from the database.

      call Grid_getLocalNumBlks(lnblocks2)

!               The PARAMESH library only supplies prolongation routines to work
!               on the "work" array.  So we copy the variable to be prolongated
!               into work(), interpolate that array, and copy the result into
!               the appropriate level of the output array.  However, we have to
!               trick the mesh package into thinking all of the "to"-level
!               blocks are new children so that the supplied prolongation
!               routine can do its work.

!               Alter the node type and new child arrays.


      iw = interp_work - 1
      do lb = 1, lnblocks2
          ! Point to blocks center vars:
          call Grid_getBlkPtr(lb,unkt,CENTER)
          work(ili-iw*1:iui+iw*1,jli-iw*k2d:jui+iw*k2d,kli-iw*k3d:kui+iw*k3d,lb,1) = & 
     &    unkt(ifrom,ili-iw*1:iui+iw*1,jli-iw*k2d:jui+iw*k2d,kli-iw*k3d:kui+iw*k3d)
          ! Point to blocks center vars:
          call Grid_releaseBlkPtr(lb,unkt,CENTER)
      enddo

      call amr_mg_prolong (NumPes2, MyPe2, level+1)


!               Now transfer the prolongated work data to the output level
!               and variable.  Add or replace as requested.

      if (add == 0) then
        do lb = 1, lnblocks2
          if (lrefine(lb) == level+1) then

          ! Point to blocks center vars:
          call Grid_getBlkPtr(lb,unkt,CENTER)
 
          unkt(ito,ili:iui,jli:jui,kli:kui) = & 
     &        work(ili:iui,jli:jui,kli:kui,lb,1)

          ! Point to blocks center vars:
          call Grid_releaseBlkPtr(lb,unkt,CENTER)

          endif
        enddo
      else
        do lb = 1, lnblocks2
          if (lrefine(lb) == level+1) then

          ! Point to blocks center vars:
          call Grid_getBlkPtr(lb,unkt,CENTER)
 
            unkt(ito,ili:iui,jli:jui,kli:kui) = & 
     &        unkt(ito,ili:iui,jli:jui,kli:kui) + & 
     &        work(ili:iui,jli:jui,kli:kui,lb,1)

          ! Point to blocks center vars:
          call Grid_releaseBlkPtr(lb,unkt,CENTER)

          endif
        enddo
      endif

!===============================================================================

      return
      end
