!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_printMatrix
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
!!  of the matrix will be checked with the index printing range.
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

  implicit none

  integer,           intent (in) :: fileUnit
  character (len=*), intent (in) :: title
  integer,           intent (in) :: rowMinMatrix, rowMaxMatrix
  integer,           intent (in) :: colMinMatrix, colMaxMatrix
  integer,           intent (in) :: rowMinPrint , rowMaxPrint
  integer,           intent (in) :: colMinPrint , colMaxPrint
  real,              intent (in) :: matrix (rowMinMatrix : rowMaxMatrix , colMinMatrix : colMaxMatrix)

  return
end subroutine ed_printMatrix
