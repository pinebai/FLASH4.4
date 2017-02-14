!!****if* source/physics/RadTrans/RadTransMain/MGD/ExpRelax/RadTrans
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This subroutine performs the radiatiative transfer calculation
!!      for this step using multigroup diffusion theory.
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!   dt     : The time step
!!   pass   : reverses solve direction
!!
!!***
!#define DEBUG_RADTRANS
subroutine RadTrans(nblk, blklst, dt, pass)

  use rt_data, ONLY: rt_useMGD
  use rt_data, ONLY: rt_mgdDomainBC
  use rt_data, ONLY: rt_mgdBounds
  use rt_data, ONLY: rt_mgdFlMode
  use rt_data, ONLY: rt_mgdFlCoef
  use rt_data, ONLY: rt_mgdNumGroups
  use rt_data, ONLY: rt_mgdBcVals
  use rt_data, ONLY: rt_mgdthetaImplct

  use rt_expData, ONLY: rt_expRelaxMaxIter

  use RadTrans_data, ONLY: rt_speedlt
  use RadTrans_data, ONLY: rt_radconst
  use RadTrans_data, ONLY: rt_meshMe
  use RadTrans_data, ONLY: rt_meshCopyCount, rt_acrossMe, rt_globalComm

  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr
  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_fillGuardCells
  use Grid_interface, ONLY: GRID_PDE_BND_PERIODIC
  use Grid_interface, ONLY: GRID_PDE_BND_NEUMANN
  use Grid_interface, ONLY: GRID_PDE_BND_DIRICHLET

  use Diffuse_interface, ONLY: Diffuse_solveScalar
  use Diffuse_interface, ONLY: Diffuse_fluxLimiter

  use Eos_interface, ONLY: Eos_wrapped
  use Eos_interface, ONLY: Eos

  use Opacity_interface, ONLY: Opacity

  use Roots_interface, ONLY: Roots_x4Polynomial

  use Timers_interface, ONLY: Timers_start
  use Timers_interface, ONLY: Timers_stop 

  implicit none

#include "Eos.h"
#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt
  integer, intent(in), optional :: pass

  ! Local variables:
  integer :: lb, gloc, gvar
  integer :: i, ii
  integer :: j
  integer :: k
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: bcTypes(6)
  integer :: nr
  integer :: ierr
  integer :: ni
  integer :: num_iter

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real :: bcValues(2,6)
  real :: absorb_opac
  real :: emit_opac
  real :: trans_opac
  real :: alpha, logAlpha
  real :: deeledurad
  real :: soln(NUNK_VARS)
  real :: cvele
  real :: utot
  real :: eelin
  real :: urlin
  real :: dens
  real :: eeold
  real :: teold
  real :: urnew
  real :: ucheck
  real :: phiold
  real :: coef0
  real :: coef1
  real :: coef2
  real :: coef3
  real :: r1
  real :: r2
  real :: max_change
  real :: global_max_change

  real :: root (1:4,1:2)

  real, pointer :: blkPtr(:,:,:,:)
  real    :: f

  real, parameter :: tolerance = 1.0e-2

  logical, parameter :: recalibrateGroupsAtStart = .FALSE.
  logical, parameter :: recalibrateGroupsAtEnd   = .TRUE.

  !=========================================================================
  if(.not. rt_useMGD) return

  call Timers_start("RadTrans") 

  if(rt_mgdNumGroups /= 1) then
     call Driver_abortFlash("[RadTrans] ExpRelax radiation only works with 1 group")
  end if

!*** NEW STUFF START ***
  if (recalibrateGroupsAtStart .OR. recalibrateGroupsAtEnd) then
  ! The specific radiation energy (ERAD_VAR) may have changed due to
  ! hydrodynamic work. The group specific internal energies must be
  ! updated accordingly. This involves adding up the energy in each
  ! group and comparing it to ERAD_VAR. I will use MGDC_VAR as a
  ! temporary storage location for the sum of the group energies.

  ! Zero out MGDC_VAR:
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              blkPtr(MGDC_VAR, i, j, k) = 0.0
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

  ! Compute the total radiation energy from the group energies:
  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     gvar = MGDR_NONREP_LOC2UNK(gloc)

     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 blkPtr(MGDC_VAR,i,j,k) = &
                      blkPtr(MGDC_VAR,i,j,k) + blkPtr(gvar,i,j,k)
              enddo
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
  end do

  call RadTrans_sumEnergy(MGDC_VAR, nblk, blklst)
  end if                          !recalibrateGroupsAtStart .OR. recalibrateGroupsAtEnd

  if (recalibrateGroupsAtStart) then
  ! Scale the group energies according to ERAD_VAR:
  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     gvar = MGDR_NONREP_LOC2UNK(gloc)

     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 if (blkPtr(MGDC_VAR,i,j,k) == 0.0) then
                    f = 0.0
                 else
                    f = blkPtr(ERAD_VAR,i,j,k) / blkPtr(MGDC_VAR,i,j,k)
                 end if
                 blkPtr(gvar,i,j,k) = blkPtr(gvar,i,j,k) * f

              enddo
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
  end do
  end if                          !recalibrateGroupsAtStart
!*** NEW STUFF END ***

  ! Convert from specific radiation energy to radiation energy density
  ! and store the old enegy density
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              blkPtr(ERAD_VAR,i,j,k) = blkPtr(ERAD_VAR,i,j,k) * blkPtr(DENS_VAR,i,j,k)
              blkPtr(UROL_VAR,i,j,k) = blkPtr(ERAD_VAR,i,j,k)
              blkPtr(EEOL_VAR,i,j,k) = blkPtr(EELE_VAR,i,j,k)
              blkPtr(TEOL_VAR,i,j,k) = blkPtr(TELE_VAR,i,j,k)
              blkPtr(UTOT_VAR,i,j,k) = blkPtr(ERAD_VAR,i,j,k) + &
                   blkPtr(EELE_VAR,i,j,k) * blkPtr(DENS_VAR,i,j,k)
              call eos_cvele(blkPtr(:,i,j,k), blkPtr(CVEL_VAR,i,j,k))
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

  ! Setup boundary conditions:
  bcTypes(:) = rt_mgdDomainBC(1,:)
  bcValues(1,:) = rt_mgdBcVals(1,:)
  bcValues(2,:) = -1.0
  
  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW .or. bcTypes == REFLECTING)
     bcTypes = GRID_PDE_BND_NEUMANN
  elsewhere (bcTypes == VACUUM)
     bcTypes = VACUUM
  end where

  num_iter = 0
  do ni = 1, rt_expRelaxMaxIter

#ifdef DEBUG_RADTRANS
     if (ni > 1) then
        if (rt_meshMe == 0) print*,'RadTrans doing iteration #',ni,' ...'
     end if
#endif

     ! Set the various terms that go into the diffusion calculation
     ! and compute the emission term:
     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        ! Loop over cells + guard cells in block:
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 ! print "(i8,3(1pe19.10))", i, &
                 !      blkPtr(EELE_VAR,i,j,k), &
                 !      blkPtr(ERAD_VAR,i,j,k)/blkPtr(DENS_VAR,i,j,k), &
                 !      blkPtr(EELE_VAR,i,j,k) + blkPtr(ERAD_VAR,i,j,k)/blkPtr(DENS_VAR,i,j,k)

                 ! At this point, the only trusted variable is:
                 ! 
                 ! blkPtr(ERAD_VAR) -> URAD^(n-1)_+
                 ! blkPtr(EELE_VAR) -> EELE^(n-1)_+
                 ! blkPtr(TELE_VAR) -> TELE^(n-1)_+
                 !
                 ! 
                 ! And we also trust:
                 !
                 ! blkPtr(UROL_VAR) -> URAD_-
                 ! blkPtr(EEOL_VAR) -> EELE_-
                 ! blkPtr(TEOL_VAR) -> TELE_-
                 !
                 ! Which are the old values, at the beginning of the
                 ! time step

                 ! Set the leading coefficient is 1 for radiation
                 ! diffusion:
                 blkPtr(DFCF_VAR, i, j, k) = 1.0

                 ! *******************************************
                 ! *   SET DIFFUSION COEF AND FLUX LIMITER   *
                 ! *******************************************

                 ! Set the value of the flux limit:
                 blkPtr(FLLM_VAR,i,j,k) = rt_mgdFlCoef * rt_speedlt * blkPtr(ERAD_VAR,i,j,k)

                 ! Get the opacity:
                 call Opacity(blkPtr(:,i,j,k), 1, absorb_opac, emit_opac, trans_opac)

                 ! Set the diffusion coefficient:
                 blkPtr(COND_VAR, i, j, k) = rt_speedlt / (3.0 * trans_opac)


                 ! **************************************************
                 ! *   UPDATE: blkPtr(EELE) -> EELE[URAD^(n-1)_+]   *
                 ! **************************************************

                 ! URAD has been updated. Now, we have to updated tele
                 ! based on the new value of erad. So in each
                 ! iteration of this loop, the only thing that changes
                 ! is urad...

                 ! Compute alpha and dphidurad...
                 soln(:) = blkPtr(:,i,j,k)
                 soln(TELE_VAR) = soln(TEOL_VAR)
                 call compute_alpha(soln, alpha, logAlpha, deeledurad)

                 ! Compute a new temperature:
                 call advance_temperature_carefully(&
                      alpha, &
                      -logAlpha, &
                      soln(ERAD_VAR), &
                      soln(TEOL_VAR), &
                      soln(TELE_VAR) )

                 ! Now compute a new electron specific internal energy:
                 call eos_eele(soln, soln(EELE_VAR))

                 ! Alternatively, assuming a const specifc heat, we
                 ! could updated using temperature:
                 ! soln(EELE_VAR) = soln(EEOL_VAR) + soln(CVEL_VAR) * (soln(TELE_VAR) - soln(TEOL_VAR))

                 ! Set:
                 ! 
                 ! blkPtr(EELE) -> EELE[URAD^(n-1)_+]
                 ! blkPtr(TELE) -> TELE[URAD^(n-1)_+]
                 blkPtr(TELE_VAR,i,j,k) = soln(TELE_VAR)
                 blkPtr(EELE_VAR,i,j,k) = soln(EELE_VAR)

                 ! Store deraddurad[URAD^(n-1)_+] and alpha[URAD^(n-1)_+]
                 !
                 ! It is needed to evaluate EELE^(n)_{+,lin} after the
                 ! diffusion solve
                 blkPtr(DEDU_VAR,i,j,k) = deeledurad
                 blkPtr(ALPH_VAR,i,j,k) = alpha
                 blkPtr(LGAL_VAR,i,j,k) = logAlpha

                 ! **************************************
                 ! *   COMPUTE REMAINING COEFFICIENTS   *
                 ! **************************************

                 ! Set the linear term:
                 blkPtr(ABSR_VAR, i, j, k) = blkPtr(DENS_VAR,i,j,k) / dt * deeledurad                 

                 ! Set the constant term:
                 blkPtr(EMIS_VAR,i,j,k) = blkPtr(DENS_VAR,i,j,k) / dt * ( &
                      blkPtr(EEOL_VAR,i,j,k) - blkPtr(EELE_VAR,i,j,k) + &
                      deeledurad*blkPtr(ERAD_VAR,i,j,k))

                 ! **********************************
                 ! *   SET INITIAL VALUE FOR ERAD   *
                 ! **********************************

                 blkPtr(EELE_VAR,i,j,k) = blkPtr(EELE_VAR,i,j,k) - &
                      deeledurad * blkPtr(ERAD_VAR,i,j,k)

                 ! ERAD_VAR is now initialized to:
                 ! blkPtr(ERAD_VAR) -> ERAD_-
                 blkPtr(ERAD_VAR,i,j,k) = blkPtr(UROL_VAR,i,j,k)

              enddo
           enddo
        enddo
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do

     ! Diffuse_fluxLimiter needs to know the value of ERAD_VAR in guard cells:
     call Grid_fillGuardCells(CENTER, ALLDIR)

     ! Apply the flux limiter:
     call Diffuse_fluxLimiter(COND_VAR, ERAD_VAR, FLLM_VAR, rt_mgdFlMode, nblk, blklst)

     ! Solve the diffusion equation for this group:
     call Diffuse_solveScalar( &
          ERAD_VAR, COND_VAR, DFCF_VAR, &
          bcTypes, bcValues, dt, 1.0, 1.0, rt_mgdthetaImplct, &
          pass, nblk,blklst, ABSR_VAR, EMIS_VAR)
     
     max_change = 0.0

     ! After this diffusion solve, ERAD_VAR has been updated. It
     ! should now have the value:
     !
     ! ERAD^(n)_{+,lin}
     !
     ! There are just two steps left:
     !
     ! 1) Compute E^(n)_+
     ! 2) Compute e^(n)_+
     do lb = 1, nblk
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)

        ! Loop over cells + guard cells in block:
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 ! ****************************************************
                 ! *   UPDATE: blkPtr(EELE_VAR) -> EELE^(n)_{+,lin}   *
                 ! ****************************************************

                 blkPtr(EELE_VAR,i,j,k) = blkPtr(EELE_VAR,i,j,k) 

                 blkPtr(EELE_VAR,i,j,k) = blkPtr(EELE_VAR,i,j,k) + &
                      blkPtr(DEDU_VAR,i,j,k) * blkPtr(ERAD_VAR,i,j,k)

                 ! blkPtr(ETOT_VAR,i,j,k) = blkPtr(EELE_VAR,i,j,k) + blkPtr(ERAD

                 ! ***************************************************

                 ! At this point, the variables are:
                 ! 
                 ! blkPtr(ERAD_VAR) -> URAD^(n)_{+,lin}
                 ! blkPtr(EELE_VAR) -> EELE^(n)_{+,lin}
                 ! blkPtr(TELE_VAR) -> TELE[URAD^(n-1)_+]
                 !
                 ! The final step involves advancing to:
                 !
                 ! blkPtr(ERAD_VAR) -> URAD^(n)_+
                 ! blkPtr(EELE_VAR) -> EELE^(n)_+
                 ! blkPtr(TELE_VAR) -> TELE^(n)_+
                 !
                 ! The total energy in each cell has been set. Since
                 ! EELE is a function of URAD in this method, there is
                 ! only one value of URAD which is consistent with
                 ! EELE and gives the correct total energy for the
                 ! cell. In the constant specific heat approximation,
                 ! this leads to a quartic algebraic equation...

                 ! print "(i8,3(1pe19.10))", i, &
                 !      blkPtr(EELE_VAR,i,j,k), &
                 !      blkPtr(ERAD_VAR,i,j,k)/blkPtr(DENS_VAR,i,j,k), &
                 !      blkPtr(EELE_VAR,i,j,k) + blkPtr(ERAD_VAR,i,j,k)/blkPtr(DENS_VAR,i,j,k)

                 ! Get the cell specific heat:
                 soln(:) = blkPtr(:,i,j,k)
                 soln(TELE_VAR) = soln(TEOL_VAR)
                 cvele = soln(CVEL_VAR)

                 ! print *, "cvele = ", cvele, rt_radconst


                 teold = blkPtr(TEOL_VAR,i,j,k)
                 eeold = blkPtr(EEOL_VAR,i,j,k)
                 urnew = blkPtr(ERAD_VAR,i,j,k)
                 dens = blkPtr(DENS_VAR,i,j,k)
                 alpha = blkPtr(ALPH_VAR,i,j,k)                 
                 logAlpha = blkPtr(LGAL_VAR,i,j,k)
                 phiold = rt_radconst * teold**4

                 eelin = blkPtr(EELE_VAR,i,j,k)
                 urlin = blkPtr(ERAD_VAR,i,j,k)
                 utot = dens*eelin + urlin

                 if (.TRUE. .OR.(alpha < 1.0 .AND. -logAlpha > 4096.0 * epsilon(logAlpha))) then
                 ! Set quartic equation coefficients:
                    if (-logAlpha .LE. 4096.0 * epsilon(logAlpha)) then
                       ! Trust logAlpha more than alpha, and use the approximation
                       ! alpha = exp(logAlpha) ~= 1+logAlpha,
                       ! if alpha is very close to 1. - KW
                       coef0 = (-logAlpha)/rt_radconst * &
                      (dens*eeold - dens*cvele*teold - utot) &
                             - alpha*phiold / rt_radconst
                    coef1 = dens*cvele*(-logAlpha)/rt_radconst
                 else if (alpha .NE. 1.0) then
                    coef0 = (1-alpha)/rt_radconst * &
                      (dens*eeold - dens*cvele*teold - utot - alpha/(1-alpha)*phiold)
                    coef1 = dens*cvele*(1-alpha)/rt_radconst
                 else 
                    coef0 =  - alpha*phiold /rt_radconst
                    coef1 =  0.0
                 end if
                 coef2 = 0.0
                 coef3 = 0.0

                 call Roots_x4Polynomial (coef3, coef2, coef1, coef0,    nr, root)

                 if(nr == 0) then
                    print*,'Roots:',(cmplx(root(ii,Re),root(ii,Im)),ii=1,4),&
                           'Coefs:',coef1,coef0
                    call Driver_abortFlash("[RadTrans] No real roots!")
                 end if

                 r1 = root (1,Re)
                 r2 = root (2,Re)

                 if(r2 >= 0) then
                    print*,'Roots:',(cmplx(root(ii,Re),root(ii,Im)),ii=1,4),&
                           'Coefs:',coef1,coef0
                    call Driver_abortFlash("[RadTrans] Too many real roots!")
                 end if

                end if

                 blkPtr(TELE_VAR,i,j,k) = r1
                 blkPtr(EELE_VAR,i,j,k) = eeold + cvele * (blkPtr(TELE_VAR,i,j,k) - teold)
                 if (alpha .NE. 1.0 .OR. logAlpha .NE. 0.0) then
                    if (-logAlpha .LE. 4096.0 * epsilon(logAlpha)) then
                       blkPtr(ERAD_VAR,i,j,k) = rt_radconst / (-logAlpha) * &
                            (blkPtr(TELE_VAR,i,j,k)**4 - alpha*teold**4)
                    else
                       blkPtr(ERAD_VAR,i,j,k) = rt_radconst / (1-alpha) * &
                            (blkPtr(TELE_VAR,i,j,k)**4 - alpha*teold**4)
                    end if
                 end if

                 ! Compute energy check:
                 ucheck = blkPtr(EELE_VAR,i,j,k) * blkPtr(DENS_VAR,i,j,k) + &
                      blkPtr(ERAD_VAR,i,j,k) 

                 if( abs((ucheck - utot)/utot) > 1.0e-6 ) then
                    ! TODO: This is very susceptible to roundoff...fix it.
                    write(*,*) ' ucheck,utot,abs(ucheck - utot)/utot = ', ucheck,utot,abs(ucheck - utot)/utot
                    call Driver_abortFlash("[RadTrans] Big error")
                 end if

                 ! print "(i8,4(1pe19.10))", i, &
                 !      blkPtr(EELE_VAR,i,j,k) * blkPtr(DENS_VAR,i,j,k), &
                 !      urnew, &
                 !      utot, ucheck

                 utot = blkPtr(UTOT_VAR,i,j,k)
                 blkPtr(UTOT_VAR,i,j,k) = blkPtr(ERAD_VAR,i,j,k) + &
                      blkPtr(EELE_VAR,i,j,k) * blkPtr(DENS_VAR,i,j,k)

                 max_change = max( max_change, &
                      abs((utot-blkPtr(UTOT_VAR,i,j,k))/(0.5*(utot+blkPtr(UTOT_VAR,i,j,k)))) )
              end do
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
     
     ! Call mpi_allreduce to find the global maximum change in the
     ! total energy:
     num_iter = num_iter + 1
     ! I am not sure if rt_globalComm is the right one to use here...
     call mpi_allreduce(max_change, global_max_change, 1, FLASH_REAL, & 
          MPI_MAX, rt_globalComm, ierr )
     if(ierr /= MPI_SUCCESS) call Driver_abortFlash("[RadTrans] Error in allreduce")
     if(global_max_change <= tolerance) exit
  end do

  if (rt_meshMe == 0) then
     if (rt_expRelaxMaxIter == 1 .OR. num_iter .NE. 1) then
        print '(a,i8)', "[RadTrans] Number of Iterations: ", num_iter
     end if
  end if

  ! ***************************
  ! *   SOLUTION CALCULATED   *
  ! ***************************

  ! Convert from radiation energy density back to specific radiation energy:
  do lb = 1, nblk
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
!*** NEW FINAL STUFF START ***
            if (recalibrateGroupsAtEnd) then
             if (recalibrateGroupsAtStart) then
  ! Reusing old urad values saved before the main loop of this routine

  ! Scale the group energies according to ERAD_VAR:
              do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
                 gvar = MGDR_NONREP_LOC2UNK(gloc)

                 if (blkPtr(UROL_VAR,i,j,k) == 0.0) then
                    f = 0.0
                 else
                    f = blkPtr(ERAD_VAR,i,j,k) / blkPtr(UROL_VAR,i,j,k)
                 end if
                 blkPtr(gvar,i,j,k) = blkPtr(gvar,i,j,k) * f
              end do
             end if
            end if
!*** NEW FINAL STUFF END ***
              blkPtr(ERAD_VAR,i,j,k) = blkPtr(ERAD_VAR,i,j,k)/blkPtr(DENS_VAR,i,j,k)
!*** NEW FINAL STUFF START ***
            if (recalibrateGroupsAtEnd) then
             if (.NOT. recalibrateGroupsAtStart) then
  ! Reusing summed MGDC_VAR values computed on entry to this routine

  ! Scale the group energies according to ERAD_VAR:
              do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
                 gvar = MGDR_NONREP_LOC2UNK(gloc)

                 if (blkPtr(MGDC_VAR,i,j,k) == 0.0) then
                    f = 0.0
                 else
                    f = blkPtr(ERAD_VAR,i,j,k) / blkPtr(MGDC_VAR,i,j,k)
                 end if
                 blkPtr(gvar,i,j,k) = blkPtr(gvar,i,j,k) * f
              end do
             end if
            end if
!*** NEW FINAL STUFF END ***
              
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)

     call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blklst(lb))     
  end do
  
  ! At this point we know the old temperature/energy and the new
  ! temperature/energy. This is a good time to compute the radtrans
  ! time step using whatever criteria we want...
  call compute_dt

  call Timers_stop("RadTrans") 
  return

contains

  ! *************************************
  ! *                                   *
  ! *     SUBROUTINE: COMPUTE_ALPHA     *
  ! *                                   *
  ! *************************************
  subroutine compute_alpha(soln_org, alpha_global, logAlphaGlobal, deeledurad)
    implicit none

    real, parameter :: dt_growth = 1.10
    real, parameter :: tau_tolerance = 0.10
    real, parameter :: dt_red_factor = 0.10

    integer, parameter :: buffer_size = 10000

    real, intent(in ) :: soln_org(NUNK_VARS)
    real, intent(out) :: alpha_global, logAlphaGlobal
    real, intent(out) :: deeledurad

    real :: alpha, logAlpha
    real :: tau
    real :: tau_old
    real :: urad
    real :: new_dt
    real :: new_time
    real :: tele
    real :: tele_new
    real :: soln(NUNK_VARS)
    real :: cur
    real :: cvele

    real, save :: tau_buffer(buffer_size)
    real, save :: alpha_buffer(buffer_size)
    real, save :: logAlphaBuffer(buffer_size)

    integer :: count
    integer :: i, j

    logical :: done

    ! These quantities do not change throughout the integration:
    urad = soln_org(ERAD_VAR)

    ! These quantities do change:
    soln(:) = soln_org(:)
    tele = soln(TELE_VAR)
       
    done = .false.

    ! *********************************
    ! *   COMPUTE INITIAL TIME STEP   *
    ! *********************************

    ! We must compute the initial time step size for the ODE
    ! integration. If the temperature does not change much (which will
    ! be the case for many cells) then we can take a single
    ! iteration...

    ! Compute the original value of tau:
    call compute_tau(soln, tau_old)

    ! Start with actual dt size:
    new_dt = dt
    tau = tau_old
    do while (.true.)

       ! Using the old value of tau, advance the temperature...
       call advance_temperature_carefully(exp(-new_dt/tau_old), new_dt/tau_old, max(urad,0.0), tele, tele_new)

       ! Compute tau using the new temperature
       soln(TELE_VAR) = max(tele_new,1.0d-80)
       call compute_tau(soln, tau)

       ! Check to see how much tau changed...
       if(abs(tau - tau_old) / (0.5*(tau + tau_old)) <= tau_tolerance) exit

       ! Tau changed a lot, reduce the dt and try again...
       new_dt = new_dt * dt_red_factor
    end do

    ! We now have a good initial time step in new_dt...

    ! *************************
    ! *   STEP THROUGH TIME   *
    ! *************************

    ! Reset the temperature:
    soln(TELE_VAR) = soln_org(TELE_VAR)
    
    ! Step forward in time...
    new_time = 0.0
    alpha_global = 1.0; logAlphaGlobal = 0.0
    count = 0
    do while(.not. done) 

       if (new_dt + new_time >= dt) then
          new_dt = dt - new_time
          done = .true.
       end if

       ! Compute tau again:
       soln(TELE_VAR) = tele
       call compute_tau(soln, tau)

       ! Advance the temperature:
       alpha = exp(-new_dt/tau)
       alpha_global = alpha_global * alpha
       logAlpha = -new_dt/tau
       logAlphaGlobal = logAlphaGlobal + logAlpha
       call advance_temperature_carefully(alpha, new_dt/tau, urad, tele, tele_new)
       tele = tele_new

       ! print "(8(1pe15.6))", dt, new_time, new_time + new_dt, new_dt, tau, urad, tele/11604.5221, alpha_global

       ! Update time and check to see if we are done:
       new_time = new_time + new_dt

       ! Compute a new dt:
       if(abs(tau_old-tau) == 0.0) then
          new_dt = new_dt*dt_growth
       else
          new_dt = min( new_dt*dt_growth, (tau_tolerance*tau/abs(tau_old-tau)) * new_dt)
       end if

       tau_old = tau

       ! Store the tau and alpha values in buffers so that they can be
       ! used to compute dphidurad:
       count = count + 1
       
       if(count > buffer_size) then
          call Driver_abortFlash("[RadTrans] buffer_size too small")
       end if

       tau_buffer(count) = tau
       alpha_buffer(count) = alpha
       logAlphaBuffer(count) = logAlpha
    end do


    ! *******************************
    ! *   COMPUTE D(EELE)/D(URAD)   *
    ! *******************************

    soln(TELE_VAR) = tele
    call compute_tau(soln, tau)
    cvele = soln(CVEL_VAR)
    
    deeledurad = 0.0
    ! print *, count
    do i = 1, count
       cur = tau_buffer(i) * (1.0 - alpha_buffer(i))

       do j = i+1, count
          cur = cur * alpha_buffer(j)
       end do

       deeledurad = deeledurad + cur 
    end do

    deeledurad = deeledurad / tau
    
    deeledurad = deeledurad * cvele / (4*rt_radconst*tele**3)

    ! print *, 1 - alpha_global, dphidurad, dphidurad / (1-alpha) * tau, tau
    
  end subroutine compute_alpha

  subroutine eos_cvele(soln, cvele)
    implicit none

    real, intent(in ) :: soln(NUNK_VARS)
    real, intent(out) :: cvele

    logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
    real, dimension(EOS_NUM) :: eos_arr
    real :: massfrac(NSPECIES)
    
    ! Set inputs:
    eos_arr(EOS_DENS) = soln(DENS_VAR)
    eos_arr(EOS_TEMP) = soln(TEMP_VAR)
    eos_arr(EOS_TEMPION) = soln(TION_VAR)
    eos_arr(EOS_TEMPELE) = soln(TELE_VAR)
    eos_arr(EOS_TEMPRAD) = soln(TRAD_VAR)

    ! Call Eos:
    mask = .false.
    mask(EOS_CVELE) = .true.
    mask(EOS_CVION) = .true.
    mask(EOS_DET)   = .true.
    massfrac = soln(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1)
    call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
    cvele = eos_arr(EOS_CVELE)

  end subroutine eos_cvele


  subroutine eos_eele(soln, eele)
    implicit none

    real, intent(in ) :: soln(NUNK_VARS)
    real, intent(out) :: eele

    logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
    real, dimension(EOS_NUM) :: eos_arr
    real :: massfrac(NSPECIES)
    
    ! Set inputs:
    eos_arr(EOS_DENS) = soln(DENS_VAR)
    eos_arr(EOS_TEMP) = soln(TEMP_VAR)
    eos_arr(EOS_TEMPION) = soln(TION_VAR)
    eos_arr(EOS_TEMPELE) = soln(TELE_VAR)
    eos_arr(EOS_TEMPRAD) = soln(TRAD_VAR)

    ! Call Eos:
    mask = .false.
    mask(EOS_EINTELE) = .true.
    massfrac = soln(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1)
    call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
    eele = eos_arr(EOS_EINTELE)

  end subroutine eos_eele


  subroutine advance_temperature(alpha, urad, tele_old, tele_new)
    implicit none

    real, intent(in ) :: alpha
    real, intent(in ) :: urad
    real, intent(in ) :: tele_old
    real, intent(out) :: tele_new

    real :: phi
   
    ! Compute a new phi and electron temperature...
    phi  = rt_radconst * tele_old**4
    if (alpha < 1.0) &
        phi = urad + alpha * (phi - urad)
    tele_new = (phi / rt_radconst)**0.25

  end subroutine advance_temperature

  subroutine advance_temperature_carefully(alpha, x, urad, tele_old, tele_new)
    implicit none

    real, intent(in ) :: alpha
    real, intent(in ) :: x
    real, intent(in ) :: urad
    real, intent(in ) :: tele_old
    real, intent(out) :: tele_new

    real :: phi
   
    ! Compute a new phi and electron temperature...
    phi  = rt_radconst * tele_old**4
    if (alpha < 1.0 .AND. x > 4096.0 * epsilon(x)) then
       phi = urad + alpha * (phi - urad)
    else
       ! Trust x more than alpha, and use the approximation
       ! alpha = exp(-x) ~= 1-x,
       ! if alpha is very close to 1. - KW
       phi = alpha*phi + urad*x
    end if
    if (.NOT.(phi .GT. 0.0)) then
       print*,'advance_temperature_carefully: about to return nonpositive phi=',phi,'!'
       print*,'alpha, x, urad, tele_old=',alpha, x, urad, tele_old
    end if
    tele_new = (phi / rt_radconst)**0.25

  end subroutine advance_temperature_carefully


  subroutine compute_tau(soln, tau)
    implicit none

    real, intent(in ) :: soln(NUNK_VARS)
    real, intent(out) :: tau

    real :: cvele
    real :: absorb_opac
    real :: emit_opac
    real :: trans_opac

    ! Compute the EOS:
    cvele = soln(CVEL_VAR)
    
    ! Compute a new opacity:
    call Opacity(soln, 1, absorb_opac, emit_opac, trans_opac)

    ! compute tau...
    tau = soln(DENS_VAR) * cvele / &
         (4*rt_speedlt*rt_radconst*absorb_opac*soln(TELE_VAR)**3)
   
  end subroutine compute_tau


  ! **********************************
  ! *                                *
  ! *     SUBROUTINE: COMPUTE_DT     *
  ! *                                *
  ! **********************************
  subroutine compute_dt
    use rt_data, ONLY: rt_precomputedDt
    use rt_data, ONLY: rt_precomputedMinLoc
    use rt_data, ONLY: rt_tempChangeRelTol
    implicit none

    real :: new_dt
    real :: new_val
    real :: old_val
    real :: dens

    rt_precomputedDt = 1.0e+100

    do lb = 1, nblk
       call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
       call Grid_getBlkPtr(blklst(lb), blkPtr)
       
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                dens = blkPtr(DENS_VAR,i,j,k)

                ! old_val = blkPtr(UROL_VAR,i,j,k) + dens * blkPtr(EEOL_VAR,i,j,k)
                ! new_val = dens * (blkPtr(ERAD_VAR,i,j,k) + blkPtr(EELE_VAR,i,j,k))

                old_val = blkPtr(TEOL_VAR,i,j,k)
                new_val = blkPtr(TELE_VAR,i,j,k)

                ! old_val = (blkPtr(UROL_VAR,i,j,k)/rt_radconst)**0.25
                ! new_val = blkPtr(TRAD_VAR,i,j,k)
                
                if (old_val .eq. new_val) then
                   new_dt = 1.0e+100
                else
                   new_dt = dt*rt_tempChangeRelTol*old_val / abs(old_val - new_val)
                endif

                if(new_dt < rt_precomputedDt) then                   
                   rt_precomputedDt = new_dt
                   rt_precomputedMinLoc(1) = i
                   rt_precomputedMinLoc(2) = j
                   rt_precomputedMinLoc(3) = k
                   rt_precomputedMinLoc(4) = blklst(lb)
                   rt_precomputedMinLoc(5) = rt_meshMe
                end if
             enddo
          end do
       end do
       call Grid_releaseBlkPtr(blklst(lb), blkPtr)
    end do
  end subroutine compute_dt

end subroutine RadTrans
