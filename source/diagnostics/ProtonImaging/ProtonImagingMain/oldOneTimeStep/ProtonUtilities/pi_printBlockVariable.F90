!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_printBlockVariable
!!
!! NAME
!!
!!  pi_printBlockVariable
!!
!! SYNOPSIS
!!
!!  call pi_printBlockVariable (integer, intent (in) :: blockID,
!!                              integer, intent (in) :: variable,
!!                              integer, intent (in) :: fileUnit)
!!
!! DESCRIPTION
!!
!!  Prints a certain cell variable content of a block on the current processor to the file
!!  associated with the unit number 'fileUnit'. To help locating the block, its bounding
!!  Box is also printed. The whole block is printed, excluding the guard cells.
!!
!! ARGUMENTS
!!
!!  blockID  : the block number ID on the current processor
!!  variable : the integer handle for the variable to be printed
!!  fileUnit : the unit number for the printout file
!!
!! NOTES
!!
!!  Only certain variables (those relevant for proton imaging) can be printed right now.
!!  Any attempt to print a different variable will result in an informal message, but the
!!  calculation will continue. The variables that can be printed are:
!!
!!               1) magnetic x-component (handle MAGX_VAR)
!!               2) magnetic y-component (handle MAGY_VAR)
!!               3) magnetic z-component (handle MAGZ_VAR)
!!               4) electric x-component (handle ELEX_VAR)
!!               5) electric y-component (handle ELEY_VAR)
!!               6) electric z-component (handle ELEZ_VAR)
!!               7) opaque solid object  (handle BDRY_VAR)
!!
!!***

subroutine pi_printBlockVariable (blockID, variable, fileUnit)

  use Driver_interface,       ONLY : Driver_abortFlash

  use Grid_interface,         ONLY : Grid_getBlkBC,          &
                                     Grid_getBlkBoundBox,    &
                                     Grid_getBlkIndexLimits, &
                                     Grid_getBlkPtr,         &
                                     Grid_getCellCoords,     &
                                     Grid_getDeltas,         &
                                     Grid_releaseBlkPtr

  use Logfile_interface,      ONLY : Logfile_stampMessage

  use pi_interface,           ONLY : pi_printMatrix

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
  integer, intent (in) :: variable
  integer, intent (in) :: fileUnit

  character (len = MAX_STRING_LENGTH) :: titleIndex
  character (len = MAX_STRING_LENGTH) :: titleTotal
  character (len = 45               ) :: titleVariable

  integer :: ib, ie, jb, je
  integer :: imaxBlock, jmaxBlock, kmaxBlock
  integer :: iminBlock, jminBlock, kminBlock
  integer :: k

  real    :: xminBlock, yminBlock, zminBlock
  real    :: xmaxBlock, ymaxBlock, zmaxBlock

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)
  real    :: bndBox      (LOW:HIGH,1:MDIM)

  real, allocatable :: matrix   (:,:)
  real, pointer     :: solnData (:,:,:,:)
!
!
!     ...Check and store the variable to be printed into text form. Return with
!        message to the logfile, if the variable cannot be printed.
!
!
  if (     variable == MAGX_VAR) then
      titleVariable = " Block variable: Magnetic x-component (MAGX) "
  else if (variable == MAGY_VAR) then
      titleVariable = " Block variable: Magnetic y-component (MAGY) "
  else if (variable == MAGZ_VAR) then
      titleVariable = " Block variable: Magnetic z-component (MAGZ) "
  else if (variable == ELEX_VAR) then
      titleVariable = " Block variable: Electric x-component (ELEX) "
  else if (variable == ELEY_VAR) then
      titleVariable = " Block variable: Electric y-component (ELEY) "
  else if (variable == ELEZ_VAR) then
      titleVariable = " Block variable: Electric z-component (ELEZ) "
  else if (variable == BDRY_VAR) then
      titleVariable = " Block variable: Opaque solid object  (BDRY) "
  else
      call Logfile_stampMessage ("pi_printBlockVariable: Block variable not recognized!")
      return
  end if
!
!
!     ...Get some block info.
!
!
  if (blockID < 1) then
      call Driver_abortFlash ("pi_printBlockVariable: Invalid block ID < 1 ")
  end if

  call Grid_getBlkBoundBox    (blockID, bndBox)
  call Grid_getBlkPtr         (blockID, solnData, CENTER)
  call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC, CENTER)

  iminBlock = blkLimits   (LOW ,IAXIS)
  imaxBlock = blkLimits   (HIGH,IAXIS)
  jminBlock = blkLimits   (LOW ,JAXIS)
  jmaxBlock = blkLimits   (HIGH,JAXIS)
  kminBlock = blkLimits   (LOW ,KAXIS)
  kmaxBlock = blkLimits   (HIGH,KAXIS)

  xminBlock = bndBox      (LOW ,IAXIS)
  xmaxBlock = bndBox      (HIGH,IAXIS)
  yminBlock = bndBox      (LOW ,JAXIS)
  ymaxBlock = bndBox      (HIGH,JAXIS)
  zminBlock = bndBox      (LOW ,KAXIS)
  zmaxBlock = bndBox      (HIGH,KAXIS)
!
!
!     ...Print block location info.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Block bounding box "
  write (fileUnit,'(/)')
  write (fileUnit,'(a,3es14.6)') " xmin, ymin, zmin = ",xminBlock,yminBlock,zminBlock
  write (fileUnit,'(a,3es14.6)') " xmax, ymax, zmax = ",xmaxBlock,ymaxBlock,zmaxBlock
  write (fileUnit,'(/)')
!
!
!     ...Print out the data in matrix form for the i,j indices, one k index at a time.
!
!
  ib = iminBlock
  ie = imaxBlock
  jb = jminBlock
  je = jmaxBlock

  allocate (matrix (ib:ie,jb:je))

  do k = kminBlock, kmaxBlock

     write (titleIndex,'(a,i0)') " Plane (i,j) cut for Index k = ",k
     titleTotal = titleVariable // titleIndex

     matrix (ib:ie,jb:je) = solnData (variable,ib:ie,jb:je,k)

     call pi_printMatrix (fileUnit,       &
                          titleTotal,     &
                          ib, ie, jb, je, &
                          ib, ie, jb, je, &
                          matrix)

  end do
!
!
!     ...Final chores.
!
!
  deallocate (matrix)

  call Grid_releaseBlkPtr (blockID, solnData, CENTER)
!
!
!     ...Ready!
!
!
  return
end subroutine pi_printBlockVariable
