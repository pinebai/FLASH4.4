!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/PrimordialChemistry_init
!!
!! NAME
!!
!!  PrimordialChemistry_init
!!
!! SYNOPSIS
!!
!!  call PrimordialChemistry_init()
!!
!! 
!! DESCRIPTION
!!
!!  Initalizes various runtime paramters for PrimordialChemistry
!!
!!  ARGUMENTS
!!
!!  PARAMETERS
!!  
!!   usePrimordialChemistry -- Boolean, True. Turns on PrimordialChemistry module
!!   useShockBurn -- Boolean, FALSE. Prob don't need
!!   pchem_algebra -- Integer, 1,[1,2]. Controls choice of linear pchem_algebra package 
!!                       used. 1=Ma28 (don't use). 2=GIFT
!!   pchem_odeStepper -- Integer, 1, [1,2]. Controls time integration routines.
!!     1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!
!!  SEE ALSO
!!  pchem_initNetwork
!!***

subroutine PrimordialChemistry_init()
  
   use PrimordialChemistry_data
   use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use PhysicalConstants_interface, ONLY: PhysicalConstants_get

   implicit none

#include "constants.h"
#include "Flash.h"

   integer :: i
   
   print *, 'INITALIZING CHEMISTRY'

   call Driver_getMype(MESH_COMM,pchem_meshMe)

   call RuntimeParameters_get("usePrimordialChemistry", pchem_usePrimordialChemistry)
   call RuntimeParameters_get("pchem_odeStepper", pchem_odeStepper)
   call RuntimeParameters_get("pchem_algebra", pchem_algebra)
   call RuntimeParameters_get("pchem_doCool", pchem_doCool)
   call RuntimeParameters_get("pchem_mCool", pchem_mCool)
   call RuntimeParameters_get("pchem_ccCase", pchem_ccCase)
   call RuntimeParameters_get("pchem_rcCase", pchem_rcCase)   

   if((pchem_algebra .lt. 1) .or. (pchem_algebra .gt. 2)) then
     write(6,*)
     write(6,*) 'only pchem_algebra=1 = ma28'
     write(6,*) 'But do not use this one!'
     write(6,*) 'pchem_algebra=2=GIFT'
     write(6,*) 'pchem_algebra is: ', pchem_algebra
     call Driver_abortFlash('ERROR in PrimordialChemistry, wrong pchem_algebra')
   end if
   if ((pchem_odeStepper .lt. 1) .or. (pchem_odeStepper .gt. 2)) then
     write(6,*)
     write(6,*) 'only pchem_odeStepper=1 = bader-deuflhard'
     write(6,*) 'and  pchem_odeStepper=2 = rosenbrock integration are valid'
     write(6,*) 'and you have specified pchem_odeStepper=', pchem_odeStepper
     write(6,*) 'error in routine PrimordialChemistry'
     call Driver_abortFlash('ERROR in PrimordialChemistry, wrong integration type')
   end if
 
  ! call RuntimeParameters_get('useShockBurn', pchem_useShockBurn)
   call pchem_initNetwork()

   call RuntimeParameters_get('pchem_fracHydrogen', pchem_fracHydrogen)
   call RuntimeParameters_get('pchem_fracHelium', pchem_fracHelium)
   call RuntimeParameters_get('pchem_fracDeuterium',pchem_fracDeuterium)
   call RuntimeParameters_get('pchem_j21',pchem_j21)
   call RuntimeParameters_get('pchem_fshh2',pchem_fshh2)
   call RuntimeParameters_get('pchem_fshhd',pchem_fshhd)

   open(unit=2,file="./Total_metals_cooling.dat")
   open(unit=3,file="./Metal_free_cooling.dat")
   open(unit=4,file="./H_HE_COOLING.dat")

   do i=1,352
      read(2,*) metal_cooling(i)
      read(3,*) hhe_cooling(i)
      read(4,*) will_hhe_cooling(i)
   enddo

   close(2)
   close(3)
   close(4)

   call PhysicalConstants_get("proton mass", amu)
   call RuntimeParameters_get("pchem_tradmin", pchem_tradmin)
   call RuntimeParameters_get("pchem_tradmax", pchem_tradmax)
   call RuntimeParameters_get("pchem_dradmin", pchem_dradmin)
   call RuntimeParameters_get("pchem_dradmax", pchem_dradmax)
   call RuntimeParameters_get("pchem_massFracH", pchem_massFracH)
   call RuntimeParameters_get("pchem_noCool", pchem_noCool)


   return

end subroutine PrimordialChemistry_init
