!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleSetInnerZoneGrid
!!
!! NAME
!!
!!  gr_mpoleSetInnerZoneGrid
!!
!! SYNOPSIS
!!
!!  call gr_mpoleSetInnerZoneGrid (integer (in)    :: nRlocal,
!!                                 integer (in)    :: nRinnerZone,
!!                                 integer (in)    :: nPinnerZone,
!!                                 real    (inout) :: RinnerZone (1:nRinnerZone))
!!
!! DESCRIPTION
!!
!!  This routine sets up the inner zone radial grid from all local inner zone
!!  radii on each processor. When exiting this routine, all processors will have
!!  a copy of the inner zone grid and its associated data.
!!
!! ARGUMENTS
!!
!!  nRlocal     : number of inner zone radii on the local processor
!!  nRinnerZone : total number of inner zone radii
!!  nPinnerZone : total number of processors containing the inner zone
!!  RinnerZone  : the collection of all inner zone radii
!!
!!***

subroutine gr_mpoleSetInnerZoneGrid (nRlocal,     &
                                     nRinnerZone, &
                                     nPinnerZone, &
                                     RinnerZone   )

  use Driver_interface,  ONLY : Driver_abortFlash

  use Grid_data,         ONLY : gr_meshMe,  &
                                gr_meshComm

  use gr_mpoleInterface, ONLY : gr_mpoleHeapsort

  use gr_mpoleData,      ONLY : gr_mpoleDrInnerZone,            &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxR,                   &
                                gr_mpoleInnerZoneMaxR,          &
                                gr_mpoleInnerZoneQmax,          &
                                gr_mpoleInnerZoneDrRadii,       &
                                gr_mpoleInnerZoneQlower,        &
                                gr_mpoleInnerZoneQupper,        &
                                gr_mpoleInnerZoneSize,          &
                                gr_mpoleInnerZoneResolution,    &
                                gr_mpoleInnerZoneResolutionInv, &
                                gr_mpoleOuterZoneExists

  implicit none

#include "Flash.h"
#include "constants.h"

  include "Flash_mpi.h"

  integer, intent (in)    :: nRlocal
  integer, intent (in)    :: nRinnerZone
  integer, intent (in)    :: nPinnerZone
  real,    intent (inout) :: RinnerZone (1:nRinnerZone)

  logical :: invokeSend

  integer :: DrUnit
  integer :: error
  integer :: messageTag
  integer :: n,nMpiRecvCalls
  integer :: nRgrid,nRreceived
  integer :: Q
  integer :: usedSpace,freeSpace

  integer :: status (MPI_STATUS_SIZE)

  real    :: prunedR, removeR
  real    :: rmax, rmaxPrev
!
!
!     ...Those processors which contain inner zone radii will send them to the
!        master processor, which collects them into the big array. After all inner
!        zone radii were collected on the master processor, he orders them and
!        sets up the inner zone radial grid. The inner zone radial grid is then
!        broadcast to all processors.
!
!        These are the individual steps taken:
!
!              1) Local processors -> send local inner zone radii
!              2) Master processor -> collect all local inner zone radii
!                                     into global inner zone radii array
!              3) Master processor -> express the true value global inner
!                                     zone radii in inner zone atomic
!                                     distance units
!              4) Master processor -> order global inner zone radii into
!                                     increasing order
!              5) Master processor -> determine # of inner zone grid radii
!              6) Master processor -> allocate inner zone grid
!              7) Master processor -> calculate inner zone grid radii
!              8) Master processor -> deallocate global inner zone radii array
!              9) All   processors -> broadcast # of inner zone grid radii
!             10) Local processors -> deallocate local inner zone radii array
!             11) Local processors -> allocate inner zone grid
!             12) All   processors -> broadcast inner zone grid
!
!
  messageTag = 1

  invokeSend = (nRlocal > 0) .and. (gr_meshMe /= MASTER_PE)

  if (invokeSend) then

      call MPI_Send  (RinnerZone,  &
                      nRlocal,     &
                      FLASH_REAL,  &
                      MASTER_PE,   &
                      messageTag,  &
                      gr_meshComm, &
                      error        )
  end if
!
!
!     ...Master's actions.
!
!
  if (gr_meshMe == MASTER_PE) then

      usedSpace = nRlocal
      freeSpace = nRinnerZone - nRlocal

      if (nRlocal > 0) then
          nMpiRecvCalls = nPinnerZone - 1
      else
          nMpiRecvCalls = nPinnerZone
      end if

      do n = 1,nMpiRecvCalls

         call MPI_Recv (RinnerZone (usedSpace+1), &
                        freeSpace,                &
                        FLASH_REAL,               &
                        MPI_ANY_SOURCE,           &
                        MPI_ANY_TAG,              &
                        gr_meshComm,              &
                        status,                   &
                        error                     )

         call MPI_get_Count (status,     &
                             FLASH_REAL, &
                             nRreceived, &
                             error       )

         usedSpace = usedSpace + nRreceived
         freeSpace = freeSpace - nRreceived
      end do 
!
!
!     ...At this point the master has the total number of inner zone radii.
!        Convert the radii to inner zone atomic distance units and sort.
!        Determine the number of different inner grid points depending on
!        the resolution given.
!
!        The fortran intrinsic function 'ceiling' cannot be used to honor
!        the inner zone resolution (i.e. the inner zone bin width). If the
!        inner zone resolution is low (of the order of 0.1 for example) the
!        use of the ceiling function is fine, but for high resolution
!        (1.e-8 for example) the ceiling function (being an integer) steps
!        out of the integer representability bounds. Therefore, it is much
!        better to use the modular function.
!
!
      RinnerZone = gr_mpoleDrInnerZoneInv * RinnerZone

      call gr_mpoleHeapsort (nRinnerZone,RinnerZone)

      removeR  = mod (RinnerZone (1) , gr_mpoleInnerZoneResolution)     ! unwanted digits to be removed
      prunedR  = RinnerZone (1) - removeR                               ! leaves the significant digits
      rmaxPrev = prunedR + gr_mpoleInnerZoneResolution                  ! puts upper bound on radius

      nRgrid = 1
      do n = 2,nRinnerZone

         removeR  = mod (RinnerZone (n) , gr_mpoleInnerZoneResolution)
         prunedR  = RinnerZone (n) - removeR
         rmax     = prunedR + gr_mpoleInnerZoneResolution

         if (rmax > rmaxPrev) then
             nRgrid = nRgrid + 1
             rmaxPrev = rmax
         end if
      end do

      if (nRgrid == 0) then
          call Driver_abortFlash ('[gr_mpoleRad3Dcartesian] ERROR: no inner zone grid radii found')
      end if
!
!
!     ...If an outer zone exists, we have found the the number of radial grid points.
!        If no outer zone exists, we eventually have to add an extra radial bin, if the
!        largest domain radius is larger than the largest inner zone radius. This is a
!        necessary precaution when cell face positions are addressed, which lay outside
!        the current determined (cell center based!) inner zone range. Collect all the
!        inner radial grid points into the inner zone radii array.
!
!
      allocate (gr_mpoleInnerZoneDrRadii (1:nRgrid+1))  ! add extra inner radial bin just in case

      removeR = mod (RinnerZone (1) , gr_mpoleInnerZoneResolution)
      prunedR = RinnerZone (1) - removeR
      rmax    = prunedR + gr_mpoleInnerZoneResolution

      gr_mpoleInnerZoneDrRadii (1) = rmax

      nRgrid = 1
      do n = 2,nRinnerZone

         removeR  = mod (RinnerZone (n) , gr_mpoleInnerZoneResolution)
         prunedR  = RinnerZone (n) - removeR
         rmax     = prunedR + gr_mpoleInnerZoneResolution

         if (rmax > gr_mpoleInnerZoneDrRadii (nRgrid)) then
             nRgrid = nRgrid + 1
             gr_mpoleInnerZoneDrRadii (nRgrid) = rmax
         end if
      end do

      if (.not. gr_mpoleOuterZoneExists) then

           rmax = RinnerZone (nRinnerZone) * gr_mpoleDrInnerZone

           if (rmax < gr_mpoleMaxR) then
               nRgrid = nRgrid + 1
               gr_mpoleInnerZoneDrRadii (nRgrid) = gr_mpoleMaxR * gr_mpoleDrInnerZoneInv
           end if

      end if

      gr_mpoleInnerZoneQmax = nRgrid

  end if
!
!
!     ...Send the inner grid to all other processors.
!
!
  call MPI_Bcast (gr_mpoleInnerZoneQmax, &
                  1,                     &
                  FLASH_INTEGER,         &
                  MASTER_PE,             &
                  gr_meshComm,           &
                  error                  )

  if (gr_meshMe /= MASTER_PE) then
      allocate   (gr_mpoleInnerZoneDrRadii (1:gr_mpoleInnerZoneQmax))
  end if

  call MPI_Bcast (gr_mpoleInnerZoneDrRadii, &
                  gr_mpoleInnerZoneQmax,    &
                  FLASH_REAL,               &
                  MASTER_PE,                &
                  gr_meshComm,              &
                  error                     )
!
!
!     ...Create the inner zone radii locator arrays to help search through
!        the inner zone grid. This is done on all processors.
!
!
  allocate (gr_mpoleInnerZoneQlower (1:gr_mpoleInnerZoneSize))
  allocate (gr_mpoleInnerZoneQupper (1:gr_mpoleInnerZoneSize))

  gr_mpoleInnerZoneQlower (:) = gr_mpoleInnerZoneQmax
  gr_mpoleInnerZoneQupper (:) = gr_mpoleInnerZoneQmax

  do Q = gr_mpoleInnerZoneQmax,1,-1
     DrUnit = int (gr_mpoleInnerZoneDrRadii (Q))
     do n = min(gr_mpoleInnerZoneSize,DrUnit+1),1,-1
        gr_mpoleInnerZoneQlower (n) = Q
     end do
  end do

  do Q = gr_mpoleInnerZoneQmax,1,-1
     DrUnit = int (gr_mpoleInnerZoneDrRadii (Q))
     do n = DrUnit,1,-1
        gr_mpoleInnerZoneQupper (n) = Q
     end do
  end do
!
!
!       ...Ready!
!
!
  return
end subroutine gr_mpoleSetInnerZoneGrid
