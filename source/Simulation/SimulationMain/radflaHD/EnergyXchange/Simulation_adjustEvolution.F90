!!****if* source/Simulation/SimulationMain/radflaHD/EnergyXchange/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

#include "constants.h"
#include "Flash.h"

  use Simulation_data

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped

  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real,    intent(in) :: dt
  real,    intent(in) :: stime

  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: axis(MDIM)
  integer :: i
  integer :: j
  integer :: k
  integer :: lb

  real, pointer :: blkPtr(:,:,:,:)
  real    :: rho
  real    :: tradActual

  if(sim_globalMe /= MASTER_PE) return

  ! Get the temperatures in any cell:
  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
!     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!
!              ! The radiation energy density is held constant, so we reset it
!              ! here
!              axis(IAXIS) = i
!              axis(JAXIS) = j
!              axis(KAXIS) = k
!              call RadTrans_mgdEFromT(blklst(lb), axis, sim_trad, tradActual)
!              blkPtr(TRAD_VAR,i,j,k) = tradActual
!           enddo
!        end do
!     end do
!     call Eos_wrapped(MODE_DENS_TEMP_GATHER, blkLimits, blklst(lb))
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! Write temperatures (and whatever else) to a file. If
              ! you add additional quantities here, you might want to
              ! also change the file header in Simulation_init.F90:
              rho = blkPtr(DENS_VAR,i,j,k)
              write(sim_fileUnitT, '(i10,7(1pe15.6))') &
                   nstep, stime, dt, &
                   blkPtr(TION_VAR,i,j,k), &
                   blkPtr(TELE_VAR,i,j,k), &
                   blkPtr(TRAD_VAR,i,j,k), &
                   blkPtr(CVEL_VAR,i,j,k), &
                   blkPtr(CVIO_VAR,i,j,k)
              write(sim_fileUnitE, '(i10,8(1pe15.6))') &
                   nstep, stime, dt, &
                   blkPtr(EION_VAR,i,j,k)*rho, &
                   blkPtr(EELE_VAR,i,j,k)*rho, &
                   blkPtr(EINT_VAR,i,j,k)*rho, &
                   blkPtr(ERAD_VAR,i,j,k)*rho, &
                   blkPtr(CVEL_VAR,i,j,k), &
                   blkPtr(CVIO_VAR,i,j,k)
              go to 100
              exit
           enddo
        end do
     end do
     100 continue
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     exit
  end do
  
end subroutine Simulation_adjustEvolution
