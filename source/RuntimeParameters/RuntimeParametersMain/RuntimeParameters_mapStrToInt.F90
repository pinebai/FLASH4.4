!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_mapStrToInt
!!
!! NAME
!!
!!  RuntimeParameters_mapStrToInt
!!
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_mapStrToInt(character(in) :: inputString(:),
!!                                     integer(out)  :: constKey)
!!
!!
!! DESCRIPTION
!!
!!  Convert a string parameter into the corresponding integer constant.
!!  The strings are defined in Config files and  provided by the flash.par file.
!!  The integer constants are defined in the header file constants.h
!!
!!  This routine is often used when mapping boundary conditions or geometry
!!  type from a string given in the flash.par to a constant key which
!!  is used by the rest of the code.
!!
!! 
!! ARGUMENTS
!!   
!!  inputString - input character string 
!!  constKey -    output integer key corresponding to inputString
!!
!! EXAMPLE
!!
!!  !  Determine the geometry requested by the flash.par
!!  call RuntimeParameters_get("geometry",pt_str_geometry)
!!  call RuntimeParameters_mapStrToInt(pt_str_geometry, pt_geometry)
!!
!!  if (pt_geometry == CARTESIAN) then
!!     .... code for rectangular domain
!!  else
!!     .... code for non-rectangular
!!  endif
!!
!!***

subroutine RuntimeParameters_mapStrToInt (inputString, constKey)

implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

character(len=*), intent(in) :: inputString
integer, intent(inout) :: constKey

constKey = NONEXISTENT

select case (inputString)

case ("periodic", "PERIODIC")
#ifdef PERIODIC
constKey = PERIODIC
#endif

case ("reflect", "reflecting", "REFLECT", "REFLECTING")
#ifdef REFLECTING
constKey = REFLECTING
#endif

case ("axisymmetric", "axisymmetry", "AXISYMMETRIC", "AXISYMMETRY")
#ifdef AXISYMMETRIC
constKey = AXISYMMETRIC
#endif

case ("eqtsymmetric", "eqtsymmetry", "EQTSYMMETRIC", "EQTSYMMETRY")
#ifdef EQTSYMMETRIC
constKey = EQTSYMMETRIC
#endif

case ("OUTFLOW", "neumann", "zero-gradient", "outflow")
#ifdef OUTFLOW
constKey = OUTFLOW
#endif

case ("diode", "DIODE")
#ifdef DIODE
constKey = DIODE
#endif

case ("DIRICHLET", "Dirichlet", "dirichlet")
#ifdef DIRICHLET
constKey = DIRICHLET
#endif

  case("NEUMANN_INS","Neumann_ins","neumann_ins")
#ifdef NEUMANN_INS
     constKey = NEUMANN_INS
#endif

  case("OUTFLOW_INS","Outflow_ins","outflow_ins")
#ifdef OUTFLOW_INS
     constKey = OUTFLOW_INS
#endif

  case("NOSLIP_INS","Noslip_ins","noslip_ins")
#ifdef NOSLIP_INS
     constKey = NOSLIP_INS
#endif

  case("SLIP_INS","Slip_ins","slip_ins")
#ifdef SLIP_INS
     constKey = SLIP_INS
#endif

  case("INFLOW_INS","Inflow_ins","inflow_ins")
#ifdef INFLOW_INS
     constKey = INFLOW_INS
#endif

  case("MOVLID_INS","Movlid_ins","movlid_ins")
#ifdef MOVLID_INS
     constKey = MOVLID_INS
#endif

case ("OUTSTREAM", "Outstream", "outstream")
#ifdef OUTSTREAM
constKey = OUTSTREAM
#endif

case ("hydrostatic", "HYDROSTATIC")
#ifdef HYDROSTATIC
constKey = HYDROSTATIC
#endif

case ("hydrostatic+nvdiode")
#ifdef HYDROSTATIC_NVDIODE
constKey = HYDROSTATIC_NVDIODE
#endif

case ("hydrostatic+nvrefl")
#ifdef HYDROSTATIC_NVREFL
constKey = HYDROSTATIC_NVREFL
#endif

case ("hydrostatic+nvout")
#ifdef HYDROSTATIC_NVOUT
constKey = HYDROSTATIC_NVOUT
#endif

case ("hydrostatic+nvzero")
#ifdef HYDROSTATIC_NVOUT
constKey = HYDROSTATIC_NVZERO
#endif

case ("hydrostatic-f2", "hydrostatic-F2", "HYDROSTATIC-F2")
#ifdef HYDROSTATIC_F2
constKey = HYDROSTATIC_F2
#endif

case ("hydrostatic-f2+nvdiode", "hydrostatic-F2+nvdiode")
#ifdef HYDROSTATIC_F2_NVDIODE
constKey = HYDROSTATIC_F2_NVDIODE
#endif

case ("hydrostatic-f2+nvrefl", "hydrostatic-F2+nvrefl")
#ifdef HYDROSTATIC_F2_NVREFL
constKey = HYDROSTATIC_F2_NVREFL
#endif

case ("hydrostatic-f2+nvout", "hydrostatic-F2+nvout")
#ifdef HYDROSTATIC_F2_NVOUT
constKey = HYDROSTATIC_F2_NVOUT
#endif

case ("extrapolate")
constKey = GRIDBC_MG_EXTRAPOLATE

case ("extrapolate+nsc")
constKey = GRIDBC_EXTRAPOLATE_NSC

case ("user", "user-defined","USER","USER-DEFINED")
#ifdef USER_DEFINED
constKey = USER_DEFINED
#endif

case ("CARTESIAN", "cartesian", "Cartesian")
#ifdef CARTESIAN
constKey = CARTESIAN
#endif

case ("polar", "POLAR")
#ifdef POLAR
constKey = POLAR
#endif

case ("cylindrical", "CYLINDRICAL", "Cylindrical")
#ifdef CYLINDRICAL
constKey = CYLINDRICAL
#endif

case ("spherical", "SPHERICAL", "Spherical")
#ifdef SPHERICAL
constKey = SPHERICAL
#endif


case("dens_ie","DENS_IE")
#ifdef MODE_DENS_EI
constKey = MODE_DENS_EI
#endif

case("dens_pres","DENS_PRES")
#ifdef MODE_DENS_PRES
constKey = MODE_DENS_PRES
#endif

case("dens_temp","DENS_TEMP")
#ifdef MODE_DENS_TEMP
constKey = MODE_DENS_TEMP
#endif

#ifdef MODE_DENS_TEMP_COMP
case("dens_temp_comp","DENS_TEMP_COMP")
constKey = MODE_DENS_TEMP_COMP
#endif

#ifdef MODE_DENS_TEMP_ALL
case("dens_temp_all","DENS_TEMP_ALL")
constKey = MODE_DENS_TEMP_ALL
#endif

#ifdef MODE_DENS_TEMP_EQUI
case("dens_temp_equi","DENS_TEMP_EQUI")
constKey = MODE_DENS_TEMP_EQUI
#endif

#ifdef MODE_DENS_TEMP_GATHER
case("dens_temp_gather","DENS_TEMP_GATHER")
constKey = MODE_DENS_TEMP_GATHER
#endif

#ifdef MODE_DENS_EI_COMP
case("dens_ie_comp","DENS_IE_COMP","dens_ei_comp","DENS_EI_COMP")
constKey = MODE_DENS_EI_COMP
#endif

#ifdef MODE_DENS_EI_ALL
case("dens_ie_all","DENS_IE_ALL","dens_ei_all","DENS_EI_ALL")
constKey = MODE_DENS_EI_ALL
#endif

#ifdef MODE_DENS_EI_EQUI
case("dens_ie_equi","DENS_IE_EQUI","dens_ei_equi","DENS_EI_EQUI")
constKey = MODE_DENS_EI_EQUI
#endif

#ifdef MODE_DENS_EI_SCATTER
case("dens_ie_scatter","DENS_IE_SCATTER","dens_ei_scatter","DENS_EI_SCATTER")
constKey = MODE_DENS_EI_SCATTER
#endif

#ifdef MODE_DENS_EI_GATHER
case("dens_ie_gather","DENS_IE_GATHER","dens_ei_gather","DENS_EI_GATHER")
constKey = MODE_DENS_EI_GATHER
#endif

#ifdef MODE_DENS_EI_RECAL_GATHER
case("dens_ie_recal_gather","DENS_IE_RECAL_GATHER","dens_ei_recal_gather","DENS_EI_RECAL_GATHER")
constKey = MODE_DENS_EI_RECAL_GATHER
#endif

#ifdef MODE_DENS_EI_SELE_GATHER
case("dens_ie_sele_gather","DENS_IE_SELE_GATHER","dens_ei_sele_gather","DENS_EI_SELE_GATHER")
constKey = MODE_DENS_EI_SELE_GATHER
#endif

#ifdef MODE_DENS_EI_SHOCKSELE_GATHER
case("dens_ie_shocksele_gather","DENS_IE_SHOCKSELE_GATHER","dens_ei_shocksele_gather","DENS_EI_SHOCKSELE_GATHER")
constKey = MODE_DENS_EI_SHOCKSELE_GATHER
#endif

#ifdef MODE_DENS_EI_MAT_GATHER_PRADSCALE
case("dens_ie_mat_gather_pradscale","DENS_IE_MAT_GATHER_PRADSCALE",&
     "dens_ei_mat_gather_pradscale","DENS_EI_MAT_GATHER_PRADSCALE")
   constKey = MODE_DENS_EI_MAT_GATHER_PRADSCALE
#endif

#ifdef MODE_EOS_NOP
case("eos_nop","EOS_NOP")
constKey = MODE_EOS_NOP
#endif

case("rt","RT")
#ifdef MODE_RT
constKey = MODE_RT
#endif

case("rp","RP")
#ifdef MODE_RP
constKey = MODE_RP
#endif

case("re","RE")
#ifdef MODE_RE
constKey = MODE_RE
#endif

case ("MARSHAK")
#ifdef MARSHAK
     constKey = MARSHAK
#endif

case ("VACUUM", "vacuum")
#ifdef VACUUM
     constKey = VACUUM
#endif


  case ("fl_none")
#ifdef FL_NONE
     constKey = FL_NONE
#endif

  case ("fl_harmonic")
#ifdef FL_HARMONIC
     constKey = FL_HARMONIC
#endif

  case ("fl_minmax")
#ifdef FL_MINMAX
     constKey = FL_MINMAX
#endif

  case ("fl_larsen")
#ifdef FL_LARSEN
     constKey = FL_LARSEN
#endif

  case ("fl_levermorepomraning1981")
#ifdef FL_LEVPOM
     constKey = FL_LEVPOM
#endif

  case ("grbd_manual")
#ifdef GRBD_MANUAL
     constKey = GRBD_MANUAL
#endif

  case ("HYPRE_AMG", "hypre_amg")
#ifdef HYPRE_AMG
     constKey = HYPRE_AMG
#endif

  case ("HYPRE_ILU", "hypre_ilu")
#ifdef HYPRE_ILU
     constKey = HYPRE_ILU
#endif

  case ("HYPRE_PCG", "hypre_pcg")
#ifdef HYPRE_PCG
     constKey = HYPRE_PCG
#endif

  case ("HYPRE_BICGSTAB", "hypre_bicgstab")
#ifdef HYPRE_BICGSTAB
     constKey = HYPRE_BICGSTAB
#endif

  case ("HYPRE_GMRES", "hypre_gmres")
#ifdef HYPRE_GMRES
     constKey = HYPRE_GMRES
#endif

  case ("HYPRE_SPLIT", "hypre_split")
#ifdef HYPRE_SPLIT
     constKey = HYPRE_SPLIT
#endif

  case ("HYPRE_PARASAILS", "hypre_parasails")
#ifdef HYPRE_PARASAILS
     constKey = HYPRE_PARASAILS
#endif     

case ("HYPRE_HYBRID", "hypre_hybrid")
#ifdef HYPRE_SPLIT
     constKey = HYPRE_HYBRID
#endif   

  case ("HYPRE_NONE", "hypre_none")
#ifdef HYPRE_NONE
     constKey = HYPRE_NONE
#endif

  case ("EOS_GAM", "eos_gam")
#ifdef EOS_GAM
     constKey = EOS_GAM
#endif

  case ("EOS_TAB", "eos_tab")
#ifdef EOS_TAB
     constKey = EOS_TAB
#endif

  case ("IONMIX", "ionmix")
     constKey = 1

  case ("IONMIX4", "ionmix4")
     constKey = 4

  case ("IONMIX6", "ionmix6")
     constKey = 6

  case ("PROPACEOS", "propaceos")      !to accomodate Propaceos
     constKey = 5                      !to accomodate Propaceos

  case ("OPACPLOT", "opacplot")
     constKey = 7

  case ("OP_UNDEFINED", "op_undefined")
#ifdef OP_UNDEFINED
     constKey = OP_UNDEFINED
#endif

  case ("OP_CONSTANT", "op_constant")
#ifdef OP_CONSTANT
     constKey = OP_CONSTANT
#endif

  case ("OP_CONSTCM2G", "op_constcm2g")
#ifdef OP_CONSTCM2G
     constKey = OP_CONSTCM2G
#endif

  case ("OP_TABPA", "op_tabpa")
#ifdef OP_TABULAR_PA
     constKey = OP_TABULAR_PA
#endif

  case ("OP_TABPE", "op_tabpe")
#ifdef OP_TABULAR_PE
     constKey = OP_TABULAR_PE
#endif

  case ("OP_TABRO", "op_tabro")
#ifdef OP_TABULAR_RO
     constKey = OP_TABULAR_RO
#endif

  case ("SMOOTH_NONE", "smooth_none")
#ifdef SMOOTH_NONE
     constKey = SMOOTH_NONE
#endif

  case ("SMOOTH_3POINT", "smooth_3point")
#ifdef SMOOTH_3POINT
     constKey = SMOOTH_3POINT
#endif

  case ("SMOOTH_3CPOINT", "smooth_3cpoint")
#ifdef SMOOTH_3CPOINT
     constKey = SMOOTH_3CPOINT
#endif

  case ("SMOOTH_SOR", "smooth_sor")
#ifdef SMOOTH_SOR
     constKey = SMOOTH_SOR
#endif

  case ("SMOOTH_HARMONIC_SOR", "smooth_harmonic_sor")
#ifdef SMOOTH_HARMONIC_SOR
     constKey = SMOOTH_HARMONIC_SOR
#endif

  case DEFAULT
     constKey = NONEXISTENT
  end select
  return
end subroutine RuntimeParameters_mapStrToInt
