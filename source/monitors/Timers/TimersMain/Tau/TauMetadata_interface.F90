!!****if* source/monitors/Timers/TimersMain/Tau/TauMetadata_interface
!!
!! NAME
!!  write_tau_metadata_str
!!
!! SYNOPSIS
!!
!!  write_tau_metadata_str (string(IN)  :: name,
!!                          string(IN)  :: value)
!!  write_tau_metadata_int (string(IN)  :: name,
!!                          integer(IN) :: value)
!!  write_tau_metadata_real(string(IN)  :: name,
!!                          real(IN)    :: value)
!!  write_tau_metadata_log (string(IN)  :: name,
!!                          logical(IN) :: value)
!!
!! DESCRIPTION
!!
!!  Write metadata fields to a Tau profile.
!!
!! ARGUMENTS
!!
!!  name  -- the name of the Tau metadata field to write/create
!!  value -- the value to associate with the specified Tau metadata field
!!
!!***

module TauMetadata_interface

implicit none

contains

!==============================================================================

subroutine write_tau_metadata_str(name, value)

  character(len=*), intent(IN) :: name, value

  call TAU_METADATA(name, value)

end subroutine write_tau_metadata_str

subroutine write_tau_metadata_int(name, value)

  character(len=*), intent(IN) :: name
  integer, intent(IN)          :: value
  character(len=20)            :: valstr

  write(valstr,'(I20)') value
  call TAU_METADATA(name, trim(valstr))

end subroutine write_tau_metadata_int

subroutine write_tau_metadata_real(name, value)

  character(len=*), intent(IN) :: name
  real, intent(IN)             :: value
  character(len=20)            :: valstr

  write(valstr,'(ES20.13)') value
  call TAU_METADATA(name, trim(valstr))

end subroutine write_tau_metadata_real

subroutine write_tau_metadata_log(name, value)

  character(len=*), intent(IN) :: name
  logical, intent(IN)          :: value
  character(len=20)            :: valstr

  write(valstr,'(L20)') value
  call TAU_METADATA(name, trim(valstr))

end subroutine write_tau_metadata_log

!==============================================================================

end module TauMetadata_interface
