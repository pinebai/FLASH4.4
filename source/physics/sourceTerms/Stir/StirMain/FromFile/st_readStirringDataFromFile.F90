!!****if* source/physics/sourceTerms/Stir/StirMain/FromFile/st_readStirringDataFromFile
!!
!! NAME
!!  st_readStirringDataFromFile
!!
!! SYNOPSIS
!!  call st_readStirringDataFromFile()
!!
!! DESCRIPTION
!!  reads the stirring data necessary to construct
!!  the physical force field from file
!!
!! ARGUMENTS
!!  infile     : file name of file containing the forcing sequence
!!  time       : current simulation time
!!  timeinfile : time in forcing file for which the modes/aka/akb are to be updated
!!
!! AUTHOR: Christoph Federrath, 2008
!!
!!***

subroutine st_readStirringDataFromFile(infile, time, timeinfile)

  use Stir_data, ONLY: st_globalMe, &
       st_mode, st_decay, st_aka, st_akb, st_energy, st_nmodes, &
       st_dtUpdateAccel, &
       st_solweight, st_solweightnorm, st_spectform, &
       st_ampl, st_stirmin, st_stirmax, &
       st_OUphases

  implicit none

#include "constants.h"

  character (len=80), intent(in) :: infile
  real, intent(in)               :: time
  real, intent(out)              :: timeinfile

  logical, parameter  :: Debug = .false.
  logical             :: opened_successful = .false.
  integer             :: dr_myPE, nsteps, step, stepinfile, desired_step, iostat
  real                :: end_time


  opened_successful = .false.
  do while (.not. opened_successful)
    open (unit=42, file=infile, iostat=iostat, status='OLD', action='READ', &
          access='SEQUENTIAL', form='UNFORMATTED')
    ! header contains number of times and number of modes, end time, autocorrelation time, ...
    if (iostat.eq.0) then
       if (Debug) write (*,'(A)') 'reading header...'
       read (unit=42) nsteps, st_nmodes, end_time, st_decay, st_energy, st_solweight, &
                      st_solweightnorm, st_stirmin, st_stirmax, st_spectform
       if (Debug) write (*,'(A)') '...finished reading header'
       opened_successful = .true.
    else
       write (*,'(A,I6,2A)') &
            '[',st_globalMe,'] st_readStirringDataFromFile: could not open file for read. filename: ', trim(infile)
       write (*,'(A,I6,A)') '[',st_globalMe,'] Trying again...'
       call sleep(1) ! wait a second...
    endif
  enddo

  ! these are in the global contex
  st_dtUpdateAccel = end_time/nsteps
  desired_step = floor(time/st_dtUpdateAccel)

  do step = 0, nsteps
  
      if (Debug) write (*,'(A,I6)') 'step = ', step
      read (unit=42) stepinfile, timeinfile, &
                     st_mode    (:, 1:  st_nmodes), &
                     st_aka     (:, 1:  st_nmodes), &
                     st_akb     (:, 1:  st_nmodes), &
                     st_ampl    (   1:  st_nmodes), &
                     st_OUphases(   1:6*st_nmodes)

      if (step.ne.stepinfile) write(*,'(A,I6)') 'st_readStirringDataFromFile: something wrong! step = ', step
      if (desired_step.eq.step) then
         if (st_globalMe == MASTER_PE) write(*,'(A,I5,2(A,ES12.5))') 'Read new forcing pattern: #', step, &
                                                                 ' time=', time, ' time_infile=', timeinfile
         close (unit=42)
         exit ! the loop
      endif
      
  enddo

  return

end subroutine st_readStirringDataFromFile
