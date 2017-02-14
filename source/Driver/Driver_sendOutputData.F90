!!****f* source/Driver/Driver_sendOutputData
!!
!! NAME
!!  Driver_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Driver_sendOutputData()
!!  
!! DESCRIPTION 
!! This routine sends the scalar variables owned by the Driver unit
!! like time and dt to the IO unit, to be written to a checkpoint file.
!!
!!
!!***


!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"


subroutine Driver_sendOutputData()


implicit none
end subroutine Driver_sendOutputData

