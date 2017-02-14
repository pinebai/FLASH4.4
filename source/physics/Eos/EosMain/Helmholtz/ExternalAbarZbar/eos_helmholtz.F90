!!****if* source/physics/Eos/EosMain/Helmholtz/ExternalAbarZbar/eos_helmholtz
!!
!! NAME
!!
!! eos_helmholtz
!!
!! SYNOPSIS
!!
!! call eos_helmholtz(integer(IN) :: mode,
!!                integer(IN) :: vecLen,
!!                real(INOUT) :: eosData(vecLen*EOS_NUM),
!!      optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!   Driver for the Helmholtz and Nadyozhin equations of state.
!!   See the NOTES section for important information about this implementation.
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  Data is packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  energy and temperature are generated as output.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!
!!  ARGUMENTS
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  massFrac : If present, contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!             Not used in this implementation.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             The Helmholtz EOS kernel calculation ignores the mask setting and calculates
!!             all derivatives, whether needed or not.  This routine does not return
!!             derivatives if the mask is requested, but the calculation is not speeded up
!!             by setting the mask.
!!
!!
!! PARAMETERS
!!
!!  eos_tol    Controls the accuracy of the Newton Rhapson iterations for MODE_DENS_EI and 
!!             MODE_DENS_PRES.
!!
!! EXAMPLE
!!
!! --- A single-point at a time example, does not calculate derivatives ---
!!
!!  use Eos_interface, ONLY:  Eos
!!  use eos_localInterface, ONLY : eos_externalComputeAbarZbar
!!  use Grid_interface ! ....
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for DENS_VAR,TEMP_VAR, etc.
!!  #include "Eos.h"         ! for EOS_NUM, EOS_TEMP, etc.
!!
!!  real  :: temp_zone, rho_zone, ptot, eint, entr, gamma
!!  real, dimension(EOS_NUM)  :: eosData
!!  real, dimension(SPECIES_BEGIN:UNK_VARS_END,1) :: xn ! NSPECIES mass fractions, followed by scalars
!!  integer, dimension(2,MDIM)                 :: blockRange,blockExtent
!!
!!
!!  !  Initialize whatever species and mass scalar data are used as inputs
!!  !  by eos_externalComputeAbarZbar and the <Unit>_computeAbarZbar
!!  !  implementations that it calls.
!!  xn(SPECIES_BEGIN:SPECIES_END,1) = 1.0e-12
!!  xn(C12_SPEC,1) = 1.0
!!  xn(MASS_SCALARS_BEGIN:UNK_VARS_END,1) = 0.0
!!  xn(YE_MSCALAR,1) = 0.5
!!
!!  !.... initialize temp_zone, rho_zone, etc.
!!
!!  call Grid_getBlkIndexLimits(blockId,blockRange,blockExtent)
!!  do k = blockRange(LOW,KAXIS), blockRange(HIGH,KAXIS)
!!     do j = blockRange(LOW,JAXIS),blockRange(HIGH,JAXIS)
!!        do i = blockRange(LOW,IAXIS),blockRange(HIGH,IAXIS)
!!
!!           eosData(EOS_TEMP) = temp_zone
!!           eosData(EOS_DENS) = rho_zone
!!
!!  !.... Now eos_externalComputeAbarZbar is called to
!!  !     compute eosData(EOS_ABAR) and eosData(EOS_ZBAR).
!!           call eos_externalComputeAbarZbar(xn, &
!!                                            eosData(EOS_ABAR:EOS_ABAR), &
!!                                            eosData(EOS_ZBAR:EOS_ZBAR))
!!           call Eos(MODE_DENS_TEMP,1,eosData)
!!
!!           ptot = eosData(EOS_PRES)
!!           eint = eosData(EOS_EINT)  
!!           entr = eosData(EOS_ENTR)
!!           gamma = eosData(EOS_GAMC)
!!           
!!           call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
!!           call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
!!           call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
!!           call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
!!               if you want ENER_VAR, calculate it from eint and kinetic energy
!!           call Grid_putPointData(blockId,CENTER,GAMC_VAR,EXTERIOR,iPosition,gamma)
!!           call Grid_putPointData(blockId,CENTER,GAME_VAR,EXTERIOR,iPosition,(ptot/(etot*sim_rhoAmbient) + 1.0))
!!
!!         enddo  ! end of k loop
!!     enddo     ! end of j loop
!!  enddo        ! end of i loop
!!
!! ------------------ Row at a time example, with derivates (based on Eos_unitTest) --------
!!
!!  use Eos_interface, ONLY:  Eos
!!  use Grid_interface ! ....
!!  #include "constants.h"   ! for MODE_DENS_TEMP, LOW,HIGH,IAXIS,JAXIS,KAXIS
!!  #include "Flash.h"       ! for DENS_VAR, TEMP_VAR, etc.
!!  #include "Eos.h"         ! for EOS_NUM, EOS_DENS, EOS_TEMP, etc.
!!  integer veclen, isize, jsize, ksize, i,j,k, e
!!  real, dimension(:), allocatable :: eosData
!!  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
!!  real, allocatable, dimension(:,:,:,:) :: derivedVariables
!!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!!
!!   ! in the Eos_unitTest, this loops over all blocks.... here is a modified snippet from inside
!!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!
!!    !  Allocate the necessary arrays for an entire block of data
!!    isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
!!    jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
!!    ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)
!!    vecLen=isize
!!    allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
!!    allocate(eosData(vecLen*EOS_NUM))
!!    mask = .true.
!!
!!    ! indices into the first location for these variables
!!    pres = (EOS_PRES-1)*vecLen
!!    dens = (EOS_DENS-1)*vecLen
!!    temp = (EOS_TEMP-1)*vecLen
!!    abar = (EOS_ABAR-1)*vecLen
!!    zbar = (EOS_ZBAR-1)*vecLen
!!
!!
!!    call Grid_getBlkPtr(blockID,solnData)
!!    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
!!
!!           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
!!           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
!!           ! Eos Helmholtz needs a good initial estimate of temperature no matter what the mode
!!           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
!!
!!           ! ... also initialize the Abar and Zbar elements of eosData:
!!
!!           call eos_externalComputeAbarZbar(solnData(SPECIES_BEGIN:UNK_VARS_END,ib:ie,j,k), &
!!                                            eosData(abar+1:abar+vecLen), &
!!                                            eosData(zbar+1:zbar+vecLen))
!!           call Eos(MODE_DENS_PRES,vecLen,eosData,mask=mask)
!!
!!           do e=EOS_VARS+1,EOS_NUM
!!              m = (e-1)*vecLen
!!              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
!!           end do
!!        end do
!!     end do
!!
!! NOTES
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, and MODE_DENS_PRES are defined in constants.h.
!!
!!  This routine is private to the Eos unit and should be called directly only
!!  from routines that are part of the Eos unit.
!!  All routines calling this routine directly must include a 
!!     use eos_localInterface
!!  statement, preferable with "ONLY" attribute, e.g.,
!!     use eos_localInterface, ONLY:  eos_helmholtz
!!
!!  Code outside of the Eos unit should call this Helmholtz implementation only
!!  indirectly, for example, by invoking the public Eos routine.
!!  Code calling the Eos routine routine must include a 
!!     use Eos_interface 
!!  statement, preferable with "ONLY" attribute, e.g.,
!!     use Eos_interface, ONLY:  Eos
!!
!!  The Helmholtz equation of state calculations are iterative for any mode other
!!  than MODE_DENS_TEMP.  Therefore, the intial estimates for temperature and density
!!  must be pretty good upon entering Eos with any other MODE_....or the calculations will
!!  not converge.
!!
!!  This algorithm uses a data table helm_table.dat which contains the coefficients for
!!  one of the interpolating algorithms.  Upon first entry to the Eos, a binary version of this
!!  table (helm_table.bdat) is created for speed of access.  This binary file should NOT be
!!  carried across machine architectures or even compilers.
!!
!!  When USE_EOS_YE is defined, this routine is replaced by the one in 
!!  physics/Eos/EosMain/Helmholtz/Ye
!!
!!  When operating in MODE_DENS_EI, the INPUT energy is updated.  This change of an input parameter
!!     can be overridden by setting the runtime parameter eos_forceConstantInput to true.
!!     Noted below, see comments prefaced with ConstantInput.
!!  Similarly, when operating in MODE_DENS_PRES, the INPUT pressure is updated.  Physicists need
!!     to be aware of this.  Similarly can be overridden with the runtime parameter/
!!
!!  The accuracy can be adjusted with the parameter eos_tol.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!
!!*** 


subroutine eos_helmholtz(mode,vecLen,eosData,massFrac,mask)

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY:  Logfile_stampMessage

  use eos_helmInterface, ONLY : eos_helm

  use eos_helmData, ONLY: eos_tol, eos_maxNewton, &
       eos_forceConstantInput, eos_larget, eos_fluffDens
  use Eos_data, ONLY: eos_smallt, eos_meshMe

  use eos_vecData, ONLY:  tempRow, denRow, etotRow, abarRow, zbarRow, &
       gamcRow, ptotRow, deaRow, dezRow, stotRow, dsdRow, dstRow, &
       detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real, INTENT(inout), dimension(vecLen*EOS_NUM) :: eosData
  real, optional,INTENT(in), dimension(vecLen*NSPECIES) :: massFrac !UNUSED
  ! must correspond to dimensions of Eos_wrapped
  logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask

  integer :: i, k
  integer :: vecBegin,vecEnd
  integer :: pres, temp, dens, gamc, eint
  integer :: abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, dea, dez, pel, ne, eta, c_v, c_p
  real    :: abarInv, zbarFrac
  real    :: told

  ! declare some local storage for the results of the Newton iteration
  real,dimension(vecLen)::  ewantRow, tnew, error,pwantRow, bnd_lo, bnd_hi
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
  real,dimension(vecLen)::  psaveRow, esaveRow
  logical  :: firstCall = .true.

  !      Fill the pipe with the initial temperature, density, and composition.
  !      The composition is parametrized by abar and zbar, which are computed
  !      from the mass fractions xn().

  ! warn user if mask is used. Helmholtz EOS doesn't use it
  if (present(mask)) then
     if (firstCall) then
        call Logfile_stampMessage("[EOS Helmholtz] WARNING!  Mask setting does not speed up Eos Helmholtz calls")
        write(*,*)"[EOS Helmholtz] WARNING!  Mask setting does not speed up Eos Helmholtz calls"
        firstCall = .false.
     end if
  end if

  vecBegin = 1
  vecEnd = vecLen

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   ! in flash2 eos_helm, this is etot
  gamc = (EOS_GAMC-1)*vecLen   ! in flash2 eos_helm, this is gamc
  abar = (EOS_ABAR-1)*vecLen   
  zbar = (EOS_ZBAR-1)*vecLen   
  entr = (EOS_ENTR-1)*vecLen

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(vecLen)
#endif

  do k = 1, vecLen

     tempRow(k)    = eosData(temp+k)
     denRow(k)     = eosData(dens+k)
     ewantRow(k)   = eosData(eint+k)   ! store desired internal energy for mode=2 case
     pwantRow(k)   = eosData(pres+k)   ! store desired pressure for mode=3 case

     abarRow(k) = eosData(abar+k)
     zbarRow(k) = eosData(zbar+k)
     ! Note in Eos.F90, we assume the user knows what he's doing.  Eos_wrapped does not.

  enddo

  !! Save the input parameters for later reference
  esaveRow = ewantRow
  psaveRow = pwantRow

  !! set wide bracket bounds
  bnd_lo(:) = eos_smallt
  bnd_hi(:) = eos_larget


  !==============================================================================

  !      MODE_DENS_TEMP  temperature and density given

  !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
  !      with the thermodynamic quantities returned by the EOS.

  if (mode==MODE_DENS_TEMP) then

     call eos_helm(1,vecLen)

     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     !==============================================================================
     !      MODE_DENS_EI  internal energy and density given

  else if (mode==MODE_DENS_EI) then

     ! Initialize the errors
     error(:) = 0.0e0

     ! Do the first eos call with all the zones in the pipe
     !  NOTE that eos_helm can ONLY operate in the equivalent of
     !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
     !  So if you send in a crappy temperature here, you'll get a crappy starting
     !  position and the iteration won't converge.
     !  Initial temperature here is what is stored in the grid, even though we 
     !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)

     ! catch really terrible guesses, should really only be from solution pollution
     do k = vecBegin, vecEnd
        if (tempRow(k) > eos_larget) tempRow(k) = eos_larget
        if (tempRow(k) < eos_smallt) tempRow(k) = eos_smallt
     enddo

     call eos_helm(vecBegin,vecEnd)
     !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow


     !  Create initial condition
     do k = vecBegin, vecEnd
        !  ewantRow is our desired EI input

        tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k))  & 
             &           / detRow(k)

        ! if we jump out of the brackets, check if there is a solution
        if ( tnew(k) < bnd_lo(k) ) then
           told = tempRow(k)
           tempRow(k) = bnd_lo(k)
           call eos_helm(k,k)
           ! assume eint is a monotonically increasing function of T
           if ( ewantRow(k) <= etotRow(k) ) then
              ! this will force exit from the loop via error=0
              tnew(k) = tempRow(k)
           else
              ! bisection
              tempRow(k) = told
              bnd_hi(k) = told
              tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
           endif
        else if ( tnew(k) > bnd_hi(k) ) then
           told = tempRow(k)
           tempRow(k) = bnd_hi(k)
           call eos_helm(k,k)
           ! assume eint is a monotonically increasing function of T
           if ( ewantRow(k) >= etotRow(k) ) then
              ! this will force exit from the loop via error=0
              tnew(k) = tempRow(k)
           else
              ! bisection
              tempRow(k) = told
              bnd_lo(k) = told
              tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
           endif
        else
           ! using derivative as an indicator of whether the old guess was
           ! an upper lower bound, update bracket range
           ! (assumes derivative is always well-behaved)
           if ( tnew(k) >= tempRow(k) ) then
               bnd_lo(k) = tempRow(k)
           else
               bnd_hi(k) = tempRow(k)
           endif
        endif

        ! Don't allow the temperature to change by more than an order of magnitude 
        ! in a single iteration
        ! (this is probably not needed anymore with bracketing, but might reduce
        !  the number of iterations in some bad cases (?) )
!        if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
!             &           10.e0*tempRow(k)
!        if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
!             &           0.1e0*tempRow(k)

        ! Compute the error
        error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

        ! Store the new temperature
        ! but not if the guess was already good enough
        ! this prevents values from shifting if the EOS is called more than once
        ! without any changes in between (e.g. in guard cells)
        if ( .not. (error(k)< eos_tol) ) then
           tempRow(k) = tnew(k)
        endif

     enddo

     ! Loop over the zones individually now
     do k = vecBegin, vecEnd
        do i = 2, eos_maxNewton
           if (error(k)< eos_tol) goto 70

           call eos_helm(k,k)

           tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k))  & 
                &              / detRow(k)
           ! if we jump out of the brackets, check if there is a solution
           if ( tnew(k) < bnd_lo(k) ) then
              told = tempRow(k)
              tempRow(k) = bnd_lo(k)
              call eos_helm(k,k)
              ! assume eint is a monotonically increasing function of T
              if ( ewantRow(k) <= etotRow(k) ) then
                 ! this will force exit from the loop via error=0
                 tnew(k) = tempRow(k)
              else
                 ! bisection
                 tempRow(k) = told
                 bnd_hi(k) = told
                 tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
              endif
           else if ( tnew(k) > bnd_hi(k) ) then
              told = tempRow(k)
              tempRow(k) = bnd_hi(k)
              call eos_helm(k,k)
              ! assume eint is a monotonically increasing function of T
              if ( ewantRow(k) >= etotRow(k) ) then
                 ! this will force exit from the loop via error=0
                 tnew(k) = tempRow(k)
              else
                 ! bisection
                 tempRow(k) = told
                 bnd_lo(k) = told
                 tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
              endif
           else
              ! using derivative as an indicator of whether the old guess was
              ! an upper lower bound, update bracket range
              ! (assumes derivative is always well-behaved)
              if ( tnew(k) >= tempRow(k) ) then
                  bnd_lo(k) = tempRow(k)
              else
                  bnd_hi(k) = tempRow(k)
              endif
           endif


           ! Don't allow the temperature to change by more than an order of magnitude 
           ! in a single iteration
           ! (this is probably not needed anymore with bracketing, but might reduce
           !  the number of iterations in some bad cases (?) )
           if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
                &              10.e0*tempRow(k)
           if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
                &              0.1e0*tempRow(k)

           ! Compute the error
           error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

           ! Store the new temperature
           tempRow(k) = tnew(k)

        end do  ! end of Newton iterations loop.  Failure drops below, success goes to 70

        ! Land here if too many iterations are needed

        print *, ' '
        print *, 'Newton-Raphson failed in subroutine eos_helmholtz'
        print *, '(e and rho as input):'
        print *, ' '
        print *, 'too many iterations', eos_maxNewton
        print *, ' '
        print *, ' temp = ', tempRow(k)
        print *, ' dens = ', denRow(k)
        print *, ' abar = ', abarRow(k)
        print *, ' zbar = ', zbarRow(k)
        print *, ' pres = ', ptotRow(k)
        print *, ' ewant= ', ewantRow(k)

        call Driver_abortFlash('[Eos] Error: too many iterations in Helmholtz Eos')


        ! Land here if the Newton iteration converged
        !  jumps out of the iterations, but then continues to the next vector location

70      continue           
!   debugging checks for reporting hitting boundaries disabled
!        if ( tempRow(k) == eos_smallt) then
!              print *, 'eos: at dens =', denRow(k), ' shooting for e = ', ewantRow(k), &
!                    'needed faking temperature at smallt= ', tempRow(k)
!        endif
!        if (bnd_hi(k) == eos_larget) then
!           if ( abs(tempRow(k) - eos_larget)/eos_larget < 2.0*eos_tol ) then
!              print *, 'eos: at dens =', denRow(k), ' shooting for e = ', ewantRow(k), &
!                    'needed faking temperature at larget= ', tempRow(k)
!           endif
!        endif

     end do

     ! Crank through the entire eos one last time

     call eos_helm(vecBegin,vecEnd)

     ! Fill the FLASH arrays with the results.  

     !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma and entropy)
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     !  Update the energy to be the true energy, instead of the energy we were trying to meet
     !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
     do k = 1, vecLen
        if (eos_forceConstantInput) then
           eosData(eint+k) = esaveRow(k)
        else
           ! check if we hit the bounds, reset energy in fluff only
           if (tempRow(k) == eos_smallt .and. denRow(k) < eos_fluffDens) then
!              print *, 'replacing energy'
              eosData(eint+k) = etotRow(k)
           else
              eosData(eint+k) = esaveRow(k)
           endif
        endif
     enddo


     !==============================================================================

     !      MODE_DENS_PRES  pressure and density given

  else if (mode==MODE_DENS_PRES) then

     ! Initialize the errors
     error(:) = 0.0e0

     ! Do the first eos call with all the zones in the pipe
     call eos_helm(vecBegin,vecEnd)

     do k = vecBegin, vecEnd

        tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k))  & 
             &           / dptRow(k)

        ! Don't allow the temperature to change by more than an order of magnitude 
        ! in a single iteration
        if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
             &           10.e0*tempRow(k)
        if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
             &           0.1e0*tempRow(k)

        ! Compute the error
        error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

        ! Store the new temperature
        tempRow(k) = tnew(k)

        ! Check if we are freezing, if so set the temperature to smallt, and adjust 
        ! the error so we don't wait for this one
        if (tempRow(k) .LT. eos_smallt) then
           tempRow(k) = eos_smallt
           error(k)    = 0.1*eos_tol
        endif

     enddo

     ! Loop over the zones individually now
     do k = vecBegin, vecEnd

        do i = 2, eos_maxNewton

           if (error(k) .LT. eos_tol) goto 170

           ! do eos only over this single item
           call eos_helm(k,k)

           tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k))  & 
                &              / dptRow(k)

           ! Don't allow the temperature to change by more than an order of magnitude 
           ! in a single iteration
           if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
                &              10.e0*tempRow(k)
           if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
                &              0.1e0*tempRow(k)

           ! Compute the error
           error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

           ! Store the new temperature
           tempRow(k) = tnew(k)

           ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
           ! the error so we don't wait for this one
           if (tempRow(k) .LT. eos_smallt) then
              tempRow(k) = eos_smallt
              error(k)    = .1*eos_tol
           endif

        enddo

        ! Land here if too many iterations are needed

        print *, ' '
        print *, 'Newton-Raphson failed in Helmholtz Eos:'
        print *, '(p and rho as input)'
        print *, ' '
        print *, 'too many iterations'
        print *, ' '
        print *, ' k    = ', k,vecBegin,vecEnd
        print *, ' temp = ', tempRow(k)
        print *, ' dens = ', denRow(k)
        print *, ' pres = ', ptotRow(k)

        call Driver_abortFlash('[Eos] Error: too many Newton-Raphson iterations in eos_helmholtz')


        ! Land here if the Newton iteration converged

170     continue

     enddo

     ! Crank through the entire eos one last time

     call eos_helm(vecBegin,vecEnd)

     ! Fill the FLASH arrays with the results.  
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
     !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
     if (eos_forceConstantInput) then
        eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
     else
        eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
     end if


     !==============================================================================

     ! Unknown EOS mode selected

  else if (mode .NE. MODE_EOS_NOP) then
     if (eos_meshMe .EQ. MASTER_PE) print*, '[eos_helmholtz] Error: unknown input mode', mode
     call Driver_abortFlash('[Eos] Error: unknown input mode in subroutine eos_helmholtz')
  end if


  ! Get the optional values
  if(present(mask)) then
     ! Entropy derivatives
     if(mask(EOS_DST)) then
        dst = (EOS_DST-1)*vecLen
        eosData(dst+1:dst+vecLen) = dstRow(1:vecLen)
     end if
     if(mask(EOS_DSD)) then
        dsd = (EOS_DSD-1)*vecLen
        eosData(dsd+1:dsd+vecLen) = dsdRow(1:vecLen)
     end if
     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        eosData(dpt+1:dpt+vecLen) = dptRow(1:vecLen)
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        eosData(dpd+1:dpd+vecLen) = dpdRow(1:vecLen)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
        eosData(det+1:det+vecLen) = detRow(1:vecLen)
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+1:ded+vecLen) = dedRow(1:vecLen)
     end if
     if(mask(EOS_DEA))then 
        dea = (EOS_DEA-1)*vecLen
        eosData(dea+1:dea+vecLen) = deaRow(1:vecLen)
     end if
     if(mask(EOS_DEZ))then 
        dez = (EOS_DEZ-1)*vecLen
        eosData(dez+1:dez+vecLen) = dezRow(1:vecLen)
     end if
     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
        eosData(pel+1:pel+vecLen) = pelRow(1:vecLen)
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
        eosData(ne+1:ne+vecLen) = neRow(1:vecLen)
     end if
     if(mask(EOS_ETA))then 
        eta = (EOS_ETA-1)*vecLen
        eosData(eta+1:eta+vecLen) = etaRow(1:vecLen)
     end if

     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           eosData(c_v+1:c_v+vecLen) = cvRow(1:vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if

     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           eosData(c_p+1:c_p+vecLen) = cpRow(1:vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
     end if
  end if




  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

  return

end subroutine eos_helmholtz


