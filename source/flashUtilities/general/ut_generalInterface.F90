!!****ih* source/flashUtilities/general/ut_generalInterface
!!
!! NAME
!!
!!  ut_generalInterface
!!
!! SYNOPSIS
!!
!!  use ut_generalInterface
!!
!! DESCRIPTION
!!
!!  Interface module for some general utilities.
!!
!!***

! Modification history:
!     Created   March 2016  KW

module ut_generalInterface

  interface
     integer function ut_getFreeFileUnit ()
     end function ut_getFreeFileUnit
  end interface

  interface
     subroutine ut_printMatrix (fileUnit,                   &
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
       real,              intent (in) :: matrix (rowMinMatrix:rowMaxMatrix , colMinMatrix:colMaxMatrix)
     end subroutine ut_printMatrix
  end interface


end module ut_generalInterface
