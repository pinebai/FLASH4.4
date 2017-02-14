!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printBlockVariable
!!
!! NAME
!!
!!  ed_printBlockVariable
!!
!! SYNOPSIS
!!
!!  call ed_printBlockVariable (integer, intent (in) :: blockID,
!!                              integer, intent (in) :: variable,
!!                              integer, intent (in) :: fileUnit)
!!
!! DESCRIPTION
!!
!!  Prints a certain cell variable content of a block on the current processor to the file
!!  associated with the unit number 'fileUnit'. To help locating the block, its bounding
!!  Box is also printed. The whole block is printed, including the guard cells.
!!
!! ARGUMENTS
!!
!!  blockID  : the block number ID on the current processor
!!  variable : the integer handle for the variable to be printed
!!  fileUnit : the unit number for the printout file
!!
!! NOTES
!!
!!  Only certain variables (those relevant for the laser deposition) can be printed
!!  right now. Any attempt to print a different variable will result in an informal
!!  message, but the calculation will continue. The variables that can be printed
!!  are:
!!
!!               1) energy deposition    (handle ed_depoVar, e.g., DEPO_VAR)
!!               2) mass density         (handle DENS_VAR)
!!               3) electron temperature (handle TELE_VAR)
!!
!!***

subroutine ed_printBlockVariable (blockID, variable, fileUnit)

  use Driver_interface,       ONLY : Driver_abortFlash

  use Grid_interface,         ONLY : Grid_getBlkBC,          &
                                     Grid_getBlkBoundBox,    &
                                     Grid_getBlkIndexLimits, &
                                     Grid_getBlkPtr,         &
                                     Grid_getCellCoords,     &
                                     Grid_getDeltas,         &
                                     Grid_releaseBlkPtr

  use Logfile_interface,      ONLY : Logfile_stampMessage

  use ed_interface,           ONLY : ed_printMatrix

  use EnergyDeposition_data,  ONLY : ed_depoVar,             &
                                     ed_depoVarName
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
  integer :: imaxGuard, jmaxGuard, kmaxGuard
  integer :: iminGuard, jminGuard, kminGuard
  integer :: k
  integer :: nGuardCellsI, nGuardCellsJ, nGuardCellsK

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
  if (     variable == ed_depoVar) then
      titleVariable = " Block variable:    Energy Deposition ("//ed_depoVarName//")."
  else if (variable == DENS_VAR) then
      titleVariable = " Block variable:         Mass Density (DENS)."
  else if (variable == TELE_VAR) then
      titleVariable = " Block variable: Electron Temperature (TELE)."
  else
      call Logfile_stampMessage ("ed_printBlockVariable: Variable not recognized!")
      return
  end if
!
!
!     ...Get some block info.
!
!
  if (blockID < 1) then
      call Driver_abortFlash ("ed_printBlockVariable: Invalid block ID < 1 ")
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

  iminGuard = blkLimitsGC (LOW ,IAXIS)
  imaxGuard = blkLimitsGC (HIGH,IAXIS)
  jminGuard = blkLimitsGC (LOW ,JAXIS)
  jmaxGuard = blkLimitsGC (HIGH,JAXIS)
  kminGuard = blkLimitsGC (LOW ,KAXIS)
  kmaxGuard = blkLimitsGC (HIGH,KAXIS)

  xminBlock = bndBox      (LOW ,IAXIS)
  xmaxBlock = bndBox      (HIGH,IAXIS)
  yminBlock = bndBox      (LOW ,JAXIS)
  ymaxBlock = bndBox      (HIGH,JAXIS)
  zminBlock = bndBox      (LOW ,KAXIS)
  zmaxBlock = bndBox      (HIGH,KAXIS)

  nGuardCellsI = iminBlock - iminGuard
  nGuardCellsJ = jminBlock - jminGuard
  nGuardCellsK = kminBlock - kminGuard
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
!     ...Print number of guard cells in each direction.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a,i2)') " Number of guard cells in i-direction = ",nGuardCellsI
  write (fileUnit,'(a,i2)') " Number of guard cells in j-direction = ",nGuardCellsJ
  write (fileUnit,'(a,i2)') " Number of guard cells in k-direction = ",nGuardCellsK
  write (fileUnit,'(/)')
!
!
!     ...Print out the data in matrix form for the i,j indices, one k index
!        at a time.
!
!
  ib = iminGuard
  ie = imaxGuard
  jb = jminGuard
  je = jmaxGuard

  allocate (matrix (ib:ie,jb:je))

  do k = kminGuard, kmaxGuard

     write (titleIndex,'(a,i0)') " Plane (i,j) cut for Index k = ",k
     titleTotal = titleVariable // titleIndex

     matrix (ib:ie,jb:je) = solnData (variable,ib:ie,jb:je,k)

     call ed_printMatrix (fileUnit,       &
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
end subroutine ed_printBlockVariable
