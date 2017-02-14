!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_init
!!
!! NAME
!!
!!  IncompNS_init
!!
!!
!! SYNOPSIS
!!
!!  call IncompNS_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_init(restart)

  use IncompNS_data
  use ImBound_interface, ONLY : ImBound_setData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm, Driver_getNstep, Driver_abortFlash
  use Grid_interface, ONLY   : Grid_getDomainBoundBox, Grid_getDomainBC
  use ins_interface, only : ins_getBulkVelocity,    &
                            ins_statsInit

  implicit none

#include "constants.h"
#include "IncompNS.h"
  logical, intent(IN) :: restart

  call Driver_getMype(MESH_COMM, ins_meshMe)
  call Driver_getNumProcs(MESH_COMM, ins_meshNumProcs)
  call Driver_getComm(MESH_COMM, ins_meshComm)  

  call RuntimeParameters_get("useIncompNS", ins_useIncompNS)
  if (.NOT. ins_useIncompNS) RETURN

  call RuntimeParameters_get("ins_cflFlg", ins_cflflg)
  call RuntimeParameters_get("cfl", ins_cfl)
  call RuntimeParameters_get("ins_isgs",ins_isgs)
  call RuntimeParameters_get("ins_invRe",ins_invRe)
  call RuntimeParameters_get("ins_sigma",ins_sigma)
  call RuntimeParameters_get("ins_dtSpec",ins_dtspec)
  call RuntimeParameters_get("ins_intSchm",ins_intschm)
  ! set scheme type flag
  select case(ins_intSchm)
     case( AB2_SCHM )
        ins_intSchm_type = INS_INTSCHM_MULTISTEP
     case( AB2_SCHM_V )
        ins_intSchm_type = INS_INTSCHM_MULTISTEP
     case( RK3_SCHM )
        ins_intSchm_type = INS_INTSCHM_RK
     case default
        call Driver_abortFlash("ins_intSchm_type not known.")
  end select
  ins_vardt(:) = ins_dtSpec

  call Driver_getNstep(ins_nstep)
  ins_restart=restart

  call RuntimeParameters_get("ins_pressureCorrect",ins_prescorr)
  ins_prescoeff = 0.
  if (ins_prescorr) ins_prescoeff = 1.

  call RuntimeParameters_get("ins_velProlongMethod",ins_prol_method)

  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) 'ins_cfl   =',ins_cfl
     write(*,*) 'ins_isgs  =',ins_isgs
     write(*,*) 'ins_invRe =',ins_invRe
     write(*,*) 'ins_sigma =',ins_sigma
     write(*,*) 'ins_dtSpec=',ins_dtspec
     write(*,*) 'ins_intSchm=',ins_intschm
     write(*,*) 'ins_prescoeff=',ins_prescoeff
     write(*,*) 'ins_velProlongMethod=',ins_prol_method
  endif

  ! Read gravity acceleration components:
  call RuntimeParameters_get("ins_gravX",ins_gravX)
  call RuntimeParameters_get("ins_gravY",ins_gravY)
  call RuntimeParameters_get("ins_gravZ",ins_gravZ)
  
  ! Read pressure gradients if necessary, constant mass simulation data:
  call RuntimeParameters_get("ins_dpdx",ins_dpdx)
  call RuntimeParameters_get("ins_dpdy",ins_dpdy)
  call RuntimeParameters_get("ins_dpdz",ins_dpdz)

  call RuntimeParameters_get("ins_constantMass",ins_constmass)
  call RuntimeParameters_get("ins_WBREF",ins_WBREF)
  call RuntimeParameters_get("ins_areaSolids",ins_area_solids)

  ! Populate old bulk velocity:
  if (ins_constmass) call ins_getBulkVelocity(ins_WBold,KAXIS)

  if(ins_meshMe==MASTER_PE) then
     if ( ins_constmass .and. (ins_WBREF .eq. 0.)) print*, 'WARNING: Constant Mass selected, but WBREF = 0.0'
     print*,'Incmp_Navier_Stokes initialized'
  endif

  call Grid_getDomainBoundBox(ins_globalDomain)
  call Grid_getDomainBC      (ins_domainBC)

  ! Define value of kinematic viscosity for IB distributed forces:
  call ImBound_setData()

  ! Call initialization of ins statistics parameters
  call ins_statsInit()

end subroutine IncompNS_init

