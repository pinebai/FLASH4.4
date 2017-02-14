!!****if* source/Simulation/SimulationMain/unitTest/Opacity/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!! SYNOPSIS
!!
!!  Simulation_initSpecies ()
!!
!! DESCRIPTION
!!
!!  This routine initializes the species and species values needed for a 
!!  unit test of the Opacity unit. This routine is called from Multispecies_init.
!!  Just one species is defined here.
!!
!!***

subroutine Simulation_initSpecies ()

  use Multispecies_interface, ONLY : Multispecies_setProperty

  implicit none

  real, allocatable :: Zvalues   (:)
  real, allocatable :: Fractions (:)

# include "Flash.h"
# include "Multispecies.h"
!
!
!       TEST Nr. 1 -> only one species consisting of one elements
!
!    ...Define the single species as the Xe atomic element
!
!             Atomic elements  : Xe
!             Atomic numbers   : 54
!             Atomic weights   : will be determined by the opacity unit
!             Number fractions : 1/1
!
!
  allocate (Zvalues   (1:1))
  allocate (Fractions (1:1))

  Zvalues   (1) = 54.
  Fractions (1) = 1.

  call Multispecies_setProperty (SP1_SPEC, A,            131.29   ) ! atomic weight
  call Multispecies_setProperty (SP1_SPEC, Z,            54.      ) ! total # of protons
  call Multispecies_setProperty (SP1_SPEC, MS_NUMELEMS,  1        )
  call Multispecies_setProperty (SP1_SPEC, MS_ZELEMS,    Zvalues  )
  call Multispecies_setProperty (SP1_SPEC, MS_FRACTIONS, Fractions)
  call Multispecies_setProperty (SP1_SPEC, MS_OPLOWTEMP, 100.     )

  deallocate (Zvalues)
  deallocate (Fractions)
!
!
!       TEST Nr. 2 -> only one species consisting of many elements
!
!    ...Define the single species (potassium hexathiocyanoplatinate IV) with chemical formula
!
!                      K [Pt(CNS) ]
!                       2        6                        
!
!             Atomic elements  : C     N     S     K    Pt
!             Atomic numbers   : 6     7    16    19    78
!             Atomic weights   : will be determined by the opacity unit
!             Number fractions : 6/21  6/21  6/21  2/21  1/21
!
!
!  allocate (Zvalues   (1:5))
!  allocate (Fractions (1:5))
!
!  Zvalues (1) = 6.
!  Zvalues (2) = 7.
!  Zvalues (3) = 16.
!  Zvalues (4) = 19.
!  Zvalues (5) = 78.
!
!  Fractions (1) = 6./21.
!  Fractions (2) = 6./21.
!  Fractions (3) = 6./21.
!  Fractions (4) = 2./21.
!  Fractions (5) = 1./21.
!
!  call Multispecies_setProperty (SP1_SPEC, A,            621.7428 ) ! total molecular weight
!  call Multispecies_setProperty (SP1_SPEC, Z,            290.     ) ! total # of protons
!  call Multispecies_setProperty (SP1_SPEC, MS_NUMELEMS,  5        )
!  call Multispecies_setProperty (SP1_SPEC, MS_ZELEMS,    Zvalues  )
!  call Multispecies_setProperty (SP1_SPEC, MS_FRACTIONS, Fractions)
!  call Multispecies_setProperty (SP1_SPEC, MS_OPLOWTEMP, 100.     )
!
!  deallocate (Zvalues)
!  deallocate (Fractions)
!
!
!       TEST Nr. 3 -> many species consisting of many (overlapping) elements
!
!    ...Define the following (ficticious) species:
!
!                      C H N
!                       2 6 3                       
!
!             Atomic elements  : C     H     N
!             Atomic numbers   : 6     1     7
!             Atomic weights   : will be determined by the opacity unit
!             Number fractions : 2/11  6/11  3/11
!
!                     Fe S C
!                       2 2 3                       
!
!             Atomic elements  : Fe    S     C
!             Atomic numbers   : 26   16     6
!             Atomic weights   : will be determined by the opacity unit
!             Number fractions : 2/7   2/7   3/7
!
!                     Zr H  N O
!                       4 10 5 7                      
!
!             Atomic elements  : Zr    H     N     O
!             Atomic numbers   : 40    1     7     8
!             Atomic weights   : will be determined by the opacity unit
!             Number fractions : 4/26  10/26 5/26  7/26
!
!                     Fe Zr C O
!                       6  3 9 2                      
!
!             Atomic elements  : Fe    Zr    C     O
!             Atomic numbers   : 26    40    6     8
!             Atomic weights   : will be determined by the opacity unit
!             Number fractions : 6/20  3/20  9/20  2/20
!
!
!  allocate (Zvalues   (1:4))
!  allocate (Fractions (1:4))
!
!  Zvalues (1) = 6.
!  Zvalues (2) = 1.
!  Zvalues (3) = 7.
!
!  Fractions (1) = 2./11.
!  Fractions (2) = 6./11.
!  Fractions (3) = 3./11.
!
!  call Multispecies_setProperty (SP1_SPEC, A,            100.0000 ) ! dummy total molecular weight
!  call Multispecies_setProperty (SP1_SPEC, Z,            200.     ) ! dummy ionization #
!  call Multispecies_setProperty (SP1_SPEC, MS_NUMELEMS,  3        )
!  call Multispecies_setProperty (SP1_SPEC, MS_ZELEMS,    Zvalues  )
!  call Multispecies_setProperty (SP1_SPEC, MS_FRACTIONS, Fractions)
!  call Multispecies_setProperty (SP1_SPEC, MS_OPLOWTEMP, 100.     )
!
!  Zvalues (1) = 26.
!  Zvalues (2) = 16.
!  Zvalues (3) = 6.
!
!  Fractions (1) = 2./7.
!  Fractions (2) = 2./7.
!  Fractions (3) = 3./7.
!
!  call Multispecies_setProperty (SP2_SPEC, A,            100.0000 ) ! dummy total molecular weight
!  call Multispecies_setProperty (SP2_SPEC, Z,            200.     ) ! dummy ionization #
!  call Multispecies_setProperty (SP2_SPEC, MS_NUMELEMS,  3        )
!  call Multispecies_setProperty (SP2_SPEC, MS_ZELEMS,    Zvalues  )
!  call Multispecies_setProperty (SP2_SPEC, MS_FRACTIONS, Fractions)
!  call Multispecies_setProperty (SP2_SPEC, MS_OPLOWTEMP, 100.     )
!
!  Zvalues (1) = 40.
!  Zvalues (2) = 1.
!  Zvalues (3) = 7.
!  Zvalues (4) = 8.
!
!  Fractions (1) =  4./26.
!  Fractions (2) = 10./26.
!  Fractions (3) =  5./26.
!  Fractions (4) =  7./26.
!
!  call Multispecies_setProperty (SP3_SPEC, A,            100.0000 ) ! dummy total molecular weight
!  call Multispecies_setProperty (SP3_SPEC, Z,            200.     ) ! dummy ionization #
!  call Multispecies_setProperty (SP3_SPEC, MS_NUMELEMS,  4        )
!  call Multispecies_setProperty (SP3_SPEC, MS_ZELEMS,    Zvalues  )
!  call Multispecies_setProperty (SP3_SPEC, MS_FRACTIONS, Fractions)
!  call Multispecies_setProperty (SP3_SPEC, MS_OPLOWTEMP, 100.     )
!
!  Zvalues (1) = 26.
!  Zvalues (2) = 40.
!  Zvalues (3) = 6.
!  Zvalues (4) = 8.
!
!  Fractions (1) = 6./20.
!  Fractions (2) = 3./20.
!  Fractions (3) = 9./20.
!  Fractions (4) = 2./20.
!
!  call Multispecies_setProperty (SP4_SPEC, A,            100.0000 ) ! dummy total molecular weight
!  call Multispecies_setProperty (SP4_SPEC, Z,            200.     ) ! dummy ionization #
!  call Multispecies_setProperty (SP4_SPEC, MS_NUMELEMS,  4        )
!  call Multispecies_setProperty (SP4_SPEC, MS_ZELEMS,    Zvalues  )
!  call Multispecies_setProperty (SP4_SPEC, MS_FRACTIONS, Fractions)
!  call Multispecies_setProperty (SP4_SPEC, MS_OPLOWTEMP, 100.     )
!
!  deallocate (Zvalues)
!  deallocate (Fractions)
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_initSpecies
