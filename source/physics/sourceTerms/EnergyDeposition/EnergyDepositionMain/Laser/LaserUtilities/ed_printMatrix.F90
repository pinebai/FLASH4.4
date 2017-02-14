!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printMatrix
!!
!! NAME
!!
!!  ed_printMatrix
!!
!! SYNOPSIS
!!
!!  call ed_printMatrix (integer,           intent (in) :: fileUnit,
!!                       character (len=*), intent (in) :: title,
!!                       integer,           intent (in) :: rowMinMatrix,
!!                       integer,           intent (in) :: rowMaxMatrix,
!!                       integer,           intent (in) :: colMinMatrix,
!!                       integer,           intent (in) :: colMaxMatrix,
!!                       integer,           intent (in) :: rowMinPrint,
!!                       integer,           intent (in) :: rowMaxPrint,
!!                       integer,           intent (in) :: colMinPrint,
!!                       integer,           intent (in) :: colMaxPrint,
!!                       real,              intent (in) :: matrix (:,:))
!!
!! DESCRIPTION
!!
!!  General matrix printing routine. The matrix is printed, together with the descriptive
!!  title, to the file associated with the unit number 'fileUnit'. The dimensions of the
!!  matrix need not to be equal to the actual matrix section printed, however the size
!!  of the matrix will be checked against the requested index printing range.
!!
!! ARGUMENTS
!!
!!  fileUnit     : the unit number for the printout file
!!  title        : the printout title
!!  rowMinMatrix : the minimum row index of the matrix
!!  rowMaxMatrix : the maximum row index of the matrix
!!  colMinMatrix : the minimum column index of the matrix
!!  colMaxMatrix : the maximum column index of the matrix
!!  rowMinMatrix : the minimum row index of the matrix section to be printed
!!  rowMaxPrint  : the maximum row index of the matrix section to be printed
!!  colMinPrint  : the minimum column index of the matrix section to be printed
!!  colMaxPrint  : the maximum column index of the matrix section to be printed
!!  matrixPrint  : the matrix to be printed
!!
!! NOTES
!!
!!  The columns of the matrix are printed such that only 10 columns are printed
!!  simultaneously on each line.
!!
!!***

subroutine ed_printMatrix (fileUnit,                   &
                           title,                      &
                           rowMinMatrix, rowMaxMatrix, &
                           colMinMatrix, colMaxMatrix, &
                           rowMinPrint , rowMaxPrint,  &
                           colMinPrint , colMaxPrint,  &
                           matrix                      )

  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

  integer,           intent (in) :: fileUnit
  character (len=*), intent (in) :: title
  integer,           intent (in) :: rowMinMatrix, rowMaxMatrix
  integer,           intent (in) :: colMinMatrix, colMaxMatrix
  integer,           intent (in) :: rowMinPrint , rowMaxPrint
  integer,           intent (in) :: colMinPrint , colMaxPrint
  real,              intent (in) :: matrix (rowMinMatrix : rowMaxMatrix , colMinMatrix : colMaxMatrix)

  integer :: block
  integer :: col, colBeg, colEnd
  integer :: nColBlocks, nColTotal
  integer :: row

  integer, parameter :: nColPerBlock = 10
!
!
!     ...Print out the title.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') title
  write (fileUnit,'(/)')
!
!
!     ...Check, if printing index ranges are ok.
!
!
  if (rowMaxPrint > rowMaxMatrix .or. rowMinPrint < rowMinMatrix) then
      call Driver_abortFlash ("ed_printMatrix: Matrix row printing range out of bounds")
  end if

  if (colMaxPrint > colMaxMatrix .or. colMinPrint < colMinMatrix) then
      call Driver_abortFlash ("ed_printMatrix: Matrix column printing range out of bounds")
  end if

  if (rowMinPrint > rowMaxPrint .or. colMinPrint > colMaxPrint) then
      call Driver_abortFlash ("ed_printMatrix: No matrix printing range!")
  end if
!
!
!     ...Everything ok. Start the printing.
!
!
  nColTotal  = colMaxPrint - colMinPrint + 1
  nColBlocks = nColTotal / nColPerBlock
!
!
!     ...Print all full column blocks.
!
!
  colEnd = colMinPrint - 1

  do block = 1,nColBlocks

     colBeg = colEnd + 1
     colEnd = colEnd + nColPerBlock

     write (fileUnit,'(/,1x,10i14)') (col , col = colBeg , colEnd)
     write (fileUnit,'(/)')

     do row = rowMinPrint ,rowMaxPrint
        write (fileUnit,'(i6,10es14.6)') row, (matrix (row,col), col = colBeg , colEnd)
     end do

  end do
!
!
!     ...Print remaining partial column block (if any).
!
!
  colBeg = colEnd + 1
  
  if (colBeg > colMaxPrint) return

  write (fileUnit,'(/,1x,10i14)') (col , col = colBeg , colMaxPrint)
  write (fileUnit,'(/)')

  do row = rowMinPrint ,rowMaxPrint
     write (fileUnit,'(i6,10es14.6)') row, (matrix (row,col), col = colBeg , colMaxPrint)
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_printMatrix
