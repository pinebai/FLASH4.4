!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_computeDt
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans_computeDt(integer(IN) :: blockID,
!!                          integer(IN) :: blkLimits(2,MDIM),
!!                          integer(IN) :: blkLimitsGC(2,MDIM),
!!                     real(IN),pointer::  solnData(:,:,:,:),   
!!                     real(OUT)   :: dt_radtrans, 
!!                     real(OUT)   :: dt_minloc(5)) 
!!  DESCRIPTION 
!!    Compute radiative transfer time step
!!
!!  ARGUMENTS
!!    blockID       --  local block ID
!!    blkLimits     --  the indices for the interior endpoints of the block
!!    blkLimitsGC   --  the indices for endpoints including the guardcells
!!    solnData      --  the physical, solution data from grid
!!    dt_radtrans   --  variable to hold timestep constraint
!!    dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!
!!
!!  J.E.Morel and R.G. McClarren, Stability of Explicit Radiation-Material
!!  Coupling in Radiative Transfer Calculations, Journal of Quantitative Spec-
!!  troscopy and Radiative Transfer, submitted July 2010.
!!
!!
!!***
subroutine RadTrans_computeDt(blockID,  blkLimits,blkLimitsGC, &
     solnData, dt_radtrans, dt_minloc)

  use RadTrans_data, ONLY: rt_speedlt, rt_radconst, rt_meshMe, &
       rt_dtfactor, rt_meshCopyCount, rt_acrossme, rt_boltz

  use rt_data, ONLY: rt_mgdNumGroups, rt_useMGD, rt_mgdBounds, rt_computeDt

  use RadTrans_interface, ONLY: RadTrans_planckInt
  use Opacity_interface, ONLY: Opacity
  use Driver_interface, ONLY: Driver_abortFlash

  use Eos_interface, ONLY : Eos, Eos_getTempDataFromVec
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: blkLimits(2,MDIM)
  integer, intent(IN) :: blkLimitsGC(2,MDIM)
  real, pointer :: solnData(:,:,:,:) 
  real, intent(INOUT) :: dt_radtrans
  integer, intent(INOUT)  :: dt_minloc(5)

  real :: C, A, dt_temp, dt_loc, KB
  
  integer :: temploc(5)

  logical :: mask(EOS_VARS+1:EOS_NUM)
  real    :: eos_arr(EOS_NUM)
  integer :: mode, vecLen
  real    :: absorb_opac, emit_opac, trans_opac
  integer :: g, gloc, gvar, i, j, k,n
  real    :: massfrac(NSPECIES)
  real    :: cvele, cveleV

  real    :: xg, xgp1, pxg, pxgp1, tele, dBdT

  real    :: eps = 1.0e-8
  
  if (.not. rt_useMGD .or. .not. rt_computeDt ) return
  
  C  = rt_speedlt
  A  = rt_radconst
  KB = rt_boltz
  
  dt_temp       = HUGE(1.0)
  dt_loc        = 0.0
  temploc(:)    = 0

#ifndef EOS_CVELE
#define EOS_CVELE EOS_CV
#endif
  mask            = .FALSE.
  mask(EOS_CVELE) = .TRUE.
  mask(EOS_CV )   = .TRUE.
  mask(EOS_DET)   = .TRUE.
  
  mode = MODE_DENS_TEMP_ELE
                
  
  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     
     gvar = MGDR_NONREP_LOC2UNK(gloc)
     
     g = NONREP_LOC2GLOB(gloc, rt_acrossMe, rt_meshCopyCount)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              !! Compute opacity.
              call Opacity(solnData(:,i,j,k), g, absorb_opac, emit_opac, trans_opac)              

              ! Compute CVELE
              vecLen = 1
!!$              eos_arr(EOS_TEMP)    = solnData(TEMP_VAR,i,j,k)
!!$              eos_arr(EOS_TEMPELE) = solnData(TELE_VAR,i,j,k)             
!!$              eos_arr(EOS_TEMPION) = solnData(TION_VAR,i,j,k)
!!$              eos_arr(EOS_TEMPRAD) = solnData(TRAD_VAR,i,j,k)
              call Eos_getTempDataFromVec(solnData(:,i,j,k),eos_arr,mode)
              eos_arr(EOS_DENS)    = solnData(DENS_VAR,i,j,k)                          
              
              ! load the mass fractions
              do n = 1, NSPECIES
                 massfrac(n) = solnData(SPECIES_BEGIN-1+n,i,j,k)
              enddo             
              
              call Eos(mode,vecLen,eos_arr,massfrac,mask)              
              cvele = eos_arr(EOS_CVELE)
              if (cvele == 0.0) cvele = eos_arr(EOS_CV)
              cveleV = cvele * eos_arr(EOS_DENS)

              !! Compute numerical dB/dT
              
              tele  = solnData(TELE_VAR,i,j,k) 
              if (tele .NE. 0.0) then
                 xg    = rt_mgdBounds(g)   / (KB*tele)
                 xgp1  = rt_mgdBounds(g+1) / (KB*tele)
                 call RadTrans_planckInt(xg, pxg)
                 call RadTrans_planckInt(xgp1, pxgp1)
              
                 dBdT = A*tele**4 * 15/PI**4*(pxgp1 - pxg)              
              
                 tele  = solnData(TELE_VAR,i,j,k) + eps*solnData(TELE_VAR,i,j,k)
                 xg    = rt_mgdBounds(g)   / (KB*tele)
                 xgp1  = rt_mgdBounds(g+1) / (KB*tele)
                 call RadTrans_planckInt(xg, pxg)
                 call RadTrans_planckInt(xgp1, pxgp1)
              
                 dBdT =  ((A*tele**4 * 15/PI**4*(pxgp1 - pxg))-dBdT)/(eps*solnData(TELE_VAR,i,j,k))

              else
                 dBdT =  0.0
              end if
              
              !! a MGD variant
              dt_loc = abs((2.0*cveleV)/  & 
                   (absorb_opac*C*(dBdT - cveleV)))
                            
              !! dt for conditional stability
!!$              dt_loc = abs((2.0*cveleV)/  & 
!!$                   (absorb_opac*C*(4.0*A*solnData(TELE_VAR,i,j,k)**3-cveleV)))
              
!!$           dt for non-oscillatory solutions.
!!$           dt_loc = cveleV/(4.0*absorb_opac*A*C*solnData(TELE_VAR,i,j,k)**3)
              
              if (dt_loc < dt_temp) then                              
                 dt_temp = dt_loc
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = k
                 temploc(4) = blockID
                 temploc(5) = rt_meshMe           
              end if              
              
           end do
        end do
     end do
     
  end do ! gloc

  dt_temp = rt_dtFactor*dt_temp  

  if (dt_temp < dt_radtrans) then
     dt_radtrans = dt_temp
     dt_minloc = temploc
  endif
  
  if(dt_radtrans <= 0.0) call Driver_abortFlash("[ RadTrans]: computed dt is not positive! Aborting!")   


  
  return
  
end subroutine RadTrans_computeDt
