!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picInit
!!
!! NAME
!!
!!  pt_picInit
!!
!! SYNOPSIS
!!
!!  call pt_picInit()
!!
!! DESCRIPTION
!!  
!!   Initialize data for the PIC implementation
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***


subroutine pt_picInit()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use pt_picData, ONLY : pt_picNsub, pt_picGam, &
       pt_picTe, pt_picResistivity, pt_picResistivityHyper, &
       pt_picPmass_1, pt_picPcharge_1, pt_picPdensity_1, &
       pt_picPtemp_1, pt_picPvelx_1, pt_picPvely_1, &
       pt_picPvelz_1, pt_picPpc_1, &
       pt_picPmass_2, pt_picPcharge_2, pt_picPdensity_2, &
       pt_picPtemp_2, pt_picPvelx_2, pt_picPvely_2, &
       pt_picPvelz_2, pt_picPpc_2, &
       pt_picRng_seed, pt_picPname_1, pt_picPname_2, pt_picCdensMin, &
       pt_picDomainBC, pt_picDomainBoundBox
  use Particles_data, ONLY : pt_meshMe,pt_meshNumProcs
  use pt_picTools, ONLY : pt_picRan_ini
  use Grid_interface, ONLY : Grid_getDomainBC, Grid_getDomainBoundBox

  implicit none
  
  call RuntimeParameters_get('pt_picPname_1', pt_picPname_1)
  call RuntimeParameters_get('pt_picPmass_1', pt_picPmass_1)
  call RuntimeParameters_get('pt_picPcharge_1', pt_picPcharge_1)
  call RuntimeParameters_get('pt_picPdensity_1', pt_picPdensity_1)
  call RuntimeParameters_get('pt_picPtemp_1', pt_picPtemp_1)
  call RuntimeParameters_get('pt_picPvelx_1', pt_picPvelx_1)
  call RuntimeParameters_get('pt_picPvely_1', pt_picPvely_1)
  call RuntimeParameters_get('pt_picPvelz_1', pt_picPvelz_1)
  call RuntimeParameters_get('pt_picPpc_1', pt_picPpc_1)



  call RuntimeParameters_get('pt_picPname_2', pt_picPname_2)
  call RuntimeParameters_get('pt_picPmass_2', pt_picPmass_2)
  call RuntimeParameters_get('pt_picPcharge_2', pt_picPcharge_2)
  call RuntimeParameters_get('pt_picPdensity_2', pt_picPdensity_2)
  call RuntimeParameters_get('pt_picPtemp_2', pt_picPtemp_2)
  call RuntimeParameters_get('pt_picPvelx_2', pt_picPvelx_2)
  call RuntimeParameters_get('pt_picPvely_2', pt_picPvely_2)
  call RuntimeParameters_get('pt_picPvelz_2', pt_picPvelz_2)
  call RuntimeParameters_get('pt_picPpc_2', pt_picPpc_2)


  call RuntimeParameters_get('pt_picTe', pt_picTe)
  call RuntimeParameters_get('pt_picResistivity', pt_picResistivity)
  call RuntimeParameters_get('pt_picResistivityHyper', pt_picResistivityHyper)
  call RuntimeParameters_get('pt_picGam', pt_picGam)

  call RuntimeParameters_get('pt_picCdensMin', pt_picCdensMin)
  call RuntimeParameters_get('pt_picNsub', pt_picNsub)

  call RuntimeParameters_get('pt_picRng_seed', pt_picRng_seed)

  call pt_picRan_ini(pt_picRng_seed, pt_meshMe, pt_meshNumProcs)

  call Grid_getDomainBC(pt_picDomainBC)
  call Grid_getDomainBoundBox(pt_picDomainBoundBox)

End subroutine pt_picInit

