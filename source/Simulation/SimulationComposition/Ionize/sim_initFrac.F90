!!****if* source/Simulation/SimulationComposition/Ionize/sim_initFrac
!!
!! NAME
!!
!!  sim_initFrac
!!
!!
!! SYNOPSIS
!!  call sim_initFrac()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the fractions of various ions.
!!
!!
!!  ARGUMENTS : There are no arguments in this subroutine
!!
!!***


subroutine sim_initFrac()

!!           WRITTEN BY: S. ORLANDO, November 2001

!
  use Simulation_speciesData, ONLY : SPEC_NUM, SPEC_NIMAX, sim_specSelected,&
       sim_specEpRatio,sim_specAbundance, sim_specNumElect, sim_specXfrac, &
       sim_ELEC, sim_HYD
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getProperty

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

!
!  N_e = 1.213 N_p for typical cosmic abundance 
!    (Cox, A.N. 2000, ``Allen's Astrophysical Quantities'', 4th ed)
!
  real :: sum_mas, ck_sum
  integer :: id, nel, nion, inn, i
!
  real, dimension(NSPECIES) :: mabar
!
  real, parameter :: epscheck = 1.0e-4
!
!
  character*255    error_message


  !mabar is appropriate name.  This is the average atomic mass
  !over all ions of a particular element.  For instance, Helium and 
  !its two ions will be assigned the same average atomic mass.
  do i=1, NSPECIES 
     call Multispecies_getProperty(SPECIES_BEGIN + (i-1), A, mabar(i))
  end do



  sum_mas = mabar(sim_HYD)+sim_specEpRatio*mabar(sim_ELEC)
!
  id = 3
  do nel = 1,SPEC_NUM
     if (sim_specSelected(nel)) then
        sum_mas = sum_mas+sim_specAbundance(nel)*mabar(id)
        nion = sim_specNumElect(nel)+1
        id = id+nion
     endif
  enddo
!
  sim_specXfrac(sim_HYD) = mabar(sim_HYD)/sum_mas
  sim_specXfrac(sim_ELEC) = sim_specEpRatio*mabar(sim_ELEC)/sum_mas
!
  ck_sum = sim_specXfrac(sim_HYD)+sim_specXfrac(sim_ELEC)
  !
  id = 3
  do nel = 1,SPEC_NUM
     if (sim_specSelected(nel)) then
        nion = sim_specNumElect(nel)+1
        do inn = 1,nion
           sim_specXfrac(id) = sim_specAbundance(nel)*mabar(id)/sum_mas
           if (inn.eq.1) ck_sum = ck_sum+sim_specXfrac(id)
           id = id+1
        enddo
     endif
  enddo
  !
  if (abs(ck_sum-1.0).gt.epscheck) then
     write(error_message,'(''Error in routine sim_initFrac: wrong sum of fluid mass fraction. ck_sum = '',f10.5)') ck_sum
     write(*,*) error_message
     call Driver_abortFlash(error_message)
  endif
  
  !
  !
  return
end subroutine sim_initFrac
