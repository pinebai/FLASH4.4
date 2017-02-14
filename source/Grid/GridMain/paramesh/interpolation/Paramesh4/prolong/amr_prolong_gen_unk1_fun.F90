#define DEBUG_DIMS 0
!*******************************************************************************
!
! Routine:      amr_prolong_gen_unk1_fun
!
! Description:
!
! This routine takes data from the array recv, originally extracted 
! from the solution array unk, and performs a prolongation operation 
! on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
! The data in recv is from a parent block and the
! result of the prolongation operation is written into unk1 for block
! isg, which is one of its children. The position of the child within the 
! parent block is specified by the ioff, joff and koff arguments.
!
! This particular prolongation is conservative, cell-averaged, triquadratic
! interpolation. It can be used for blocks with an even or odd number of grid
! cells.
!
!
! we want to find the value of each zone-average variable in the children
! of zone ip.
! This particular interpolation is monotonic and and the order can be chosen
! at run time. The default is parabolic.

!!REORDER(5): unk1
!!REORDER(4): recv

#include "Flash.h"

subroutine amr_prolong_gen_unk1_fun & 
     &       (recv, ia, ib, ja, jb, ka, kb, idest, & 
     &        ioff, joff, koff, mype, isg)
  
  
  use physicaldata, ONLY : unk1, gcell_on_cc
  use Grid_interface, ONLY : Grid_getCellCoords
  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_data,ONLY: gr_convertToConsvdForMeshCalls, gr_convertToConsvdInMeshInterp, &
       gr_vartypes, gr_dirGeom, gr_smallx, gr_intpol
  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, intent(IN), &
       dimension(NUNK_VARS,GRID_ILO_GC:GRID_IHI_GC,&
       GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC)::recv

  integer, parameter :: niver = NUNK_VARS
  integer,parameter :: lw = 4*(NXB+2*NGUARD+4) + K2D*(4*(NYB+2*NGUARD +4)+&
       5*(NXB+2*NGUARD+4)*NUNK_VARS) + K3D*(4*(NZB+2*NGUARD+4)+&
       5*(NXB+2*NGUARD+4)*(NYB+2*NGUARD+4)*NUNK_VARS)
  integer,parameter :: liw=  NXB+NYB+NZB+max(NDIM*2*NGUARD,12+2*NGUARD)
  
  integer,parameter :: ref_ratio = 2
  integer,parameter :: intpol_guard = 2
  integer,parameter :: twice_iguard = 2*intpol_guard

  integer,parameter :: lxiend = (NXB+2*NGUARD)/2 + twice_iguard
!   NGUARD  lxiend        lxiend(NXB=8) lxiend(NXB=16) lxoend lxoend(NXB=8) lxoend(NXB=16)
!   =========================================================================================
!     2     NXB/2 + 2 + 4    10             14         NXB+4     12           20
!     4     NXB/2 + 4 + 4    12             16         NXB+8     16           24
!     6     NXB/2 + 6 + 4    14             18         NXB+12    20           28
  integer,parameter :: linxu  = lxiend*niver
  integer,parameter :: lyiend = max(1,K2D*((NYB+2*NGUARD)/2 + twice_iguard))
  integer,parameter :: lziend = max(1,K3D*((NZB+2*NGUARD)/2 + twice_iguard))

  integer,parameter :: lxoend = NXB+2*NGUARD
  integer,parameter :: lonxu  = lxoend*niver
  integer,parameter :: lyoend = NYB+2*NGUARD
  integer,parameter :: lzoend = NZB+2*NGUARD

  integer,parameter :: linxu_s = lxiend*niver
  integer,parameter :: lonxu_s = lxoend*niver
  
  integer, intent(IN) :: idest, isg, mype, ia, ib, ja, jb, ka, kb,ioff,joff,koff

  integer :: i,j,k,ivar,last_var
  integer :: i1,j1,k1,offia,offib,offic,offoa,offob,offoc
  
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: xo, xi
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: yo, yi
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: zo, zi

  real, dimension(linxu,lyiend,lziend) :: pu
  real, dimension(lonxu,lyoend,lzoend) :: qu

  real, dimension(lw)  :: wmap     ! the two work arrays needed by the interpolation routine
  integer, dimension(liw) :: iwmap ! This used to erroneously say real, for a long time
  
  real          :: pdx, pdy, pdz, cdx, cdy, cdz ! half grid spacings for parent and child
  real          :: dens_inv
  ! variable to indicate whether mass fractions have been converted to conserved form
  ! (i.e., partial densities) when UMAPn is called - KW
  logical       :: conserved_var

  logical, dimension(NUNK_VARS) :: toUmapMask, toCallerMask
  integer, dimension(NUNK_VARS) :: varIndexPacked, varIndexUnpacked
  logical :: need_relatives
  integer :: numVarsPacked, ivarPacked, densIndexPacked, ui

  integer       :: imapcx, ierr  
  integer       :: xoend, yoend, zoend, xiend, yiend, ziend
  integer       :: ii, jj, kk, n
  integer       :: inxu, onxu
  integer, dimension(MDIM) :: geom
  logical       :: lumap(NUNK_VARS)
  integer       :: mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m
  integer       :: n_dens, nu_d, nud, iud, i_dens, nui, iu
  integer       :: ip_i_iiu, ip_ii_iiu, ip_jj_iiujj, ip_iiu_iiujj

  !===============================================================================

#include "umap.h"
  
  COMMON  /amrmapi_p/                   &
       n_dens, nud(mui_m), iud(0:mui_m,mui_m),     &
       ip_i_iiu     (mvxu_m),          &
       ip_ii_iiu    (mvxu_m),          &
       ip_jj_iiujj  (mvxmu_m),         &
       ip_iiu_iiujj (mvxmu_m)
  
  !=====================================================================
#ifdef DEBUG_GRID
  if ( max(NXB,NYB,NZB) + 2*NGUARD > mvx_m ) then
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: umap.h workspace too small; mvx_m < max(NXB,NYB,NZB) + 2*NGUARD")
  end if
  
  if ( mui_m < NUNK_VARS ) then
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: umap.h workspace too small; mui_m < FLASH_NUMBER_OF_VARIABLES ")
  end if
  
  if(NGUARD<twice_iguard) then
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: not enough guardcells to support asked for interpolation")
  end if
  
#endif

  conserved_var = (gr_convertToConsvdInMeshInterp .or. gr_convertToConsvdForMeshCalls)
  ! total densities are derived from partial densities
  
  ! identify mass-like species

  n_dens = 0
  geom=gr_dirGeom
#ifdef DENS_VAR
  if ( NSPECIES > 1 ) then

     ! in general we may carry many density-like groups of species

     n_dens = n_dens + 1
     i_dens = DENS_VAR
     iud(0,n_dens) = i_dens
     
     nu_d = 0
     
     ! this is the number of variables to be interpolated; this implementation
     ! does not make any distinction, all variables are interpolated (lumap=true)

     nui = NUNK_VARS
     lumap = .true.

     do iu = 1,nui
        !if ( lumap(iu) .and. i_uftype(iu).eq.i_dens ) then
        if ( iu >= SPECIES_BEGIN .and. iu <= SPECIES_END  ) then
           nu_d = nu_d + 1
           iud(nu_d,n_dens) = iu
        end if
     end do
     
     nud(n_dens) = nu_d
     
  end if
#endif
  
  !! initialize all values that may not be used in non-3D case
  pdy = 0.e0
  cdy = 0.e0
  pdz = 0.e0
  cdz = 0.e0
  zoend=1
  yoend=1
  ziend=1
  yiend=1
  offob = 0
  offoc = 0
  offib = 0
  offic = 0
  
!!$  offia  = ia/2 + NGUARD - twice_iguard
  offia  = (ia-1+NGUARD)/ref_ratio - intpol_guard
  offoa  = ia - 1
!!WAS:
!   NGUARD  offia         offia(ia=NGUARD+1) offia(ia=1) offoa offoa(ia=NGUARD+1) offoa(ia=1)
!   =========================================================================================
!     2     ia/2  + 2 - 4   -1                -2         ia-1      2               0
!     4     ia/2  + 4 - 4    2                 0         ia-1      4               0
!     6     ia/2  + 6 - 4    5                 2         ia-1      6               0
!!NOW:
!   NGUARD  offia         offia(ia=NGUARD+1) offia(ia=1) offoa offoa(ia=NGUARD+1) offoa(ia=1)
!   =========================================================================================
!     2     (ia+1)/2 - 2     0                -1         ia-1      2               0
!     3     (ia+2)/2 - 2     1                -1         ia-1      3               0
!     4     (ia+3)/2 - 2     2                 0         ia-1      4               0
!     5     (ia+4)/2 - 2     3                 0         ia-1      5               0
!     6     (ia+5)/3 - 2     4                 1         ia-1      6               0
  
  xoend = ib-ia+1
!!$  xiend = xoend/ref_ratio + twice_iguard
  xiend = (ib-ia+ref_ratio)/ref_ratio + min(1,mod(ib-ia,ref_ratio)*mod(ib+NGUARD,ref_ratio)) + twice_iguard
!!WAS:
!   NGUARD   xiend         xiend(len([ia..lb]=NXB)  xoend  xoend(ia/b=NGUARD+1..NGUARD+NXB) xoend(ia/b=1..NGUARD)
!   ================================================================================================================
!     2     xoend/2 + 4      NXB/2 + 4  (8 or 12)   ib-ia+1   NXB                           2
!     4     xoend/2 + 4      NXB/2 + 4  (8 or 12)   ib-ia+1   NXB                           4
!     6     xoend/2 + 4      NXB/2 + 4 ([8 or]12)   ib-ia+1   NXB                           6
!!NOW:...
  
  call Grid_getCellCoords(IAXIS,isg,CENTER,.true.,xi,GRID_IHI_GC)
  pdx = xi(2) - xi(1)
  cdx = 0.5e0*pdx

  do i = 1,xoend
     xo(i) = xi(offoa+i)
  end do
  
  xi(1) = xo(1) + cdx - twice_iguard*pdx 

  do i = 2,xiend
     xi(i) = xi(i-1)+2.e0*pdx
  end do
#if N_DIM > 1

!!$  offib = ja/2 + NGUARD - twice_iguard
  offib  = (ja-1+NGUARD)/ref_ratio - intpol_guard
  offob = ja - 1

  yoend = jb-ja+1
!!$  yiend = yoend/ref_ratio + twice_iguard
  yiend = (jb-ja+ref_ratio)/ref_ratio + min(1,mod(jb-ja,ref_ratio)*mod(jb+NGUARD,ref_ratio)) + twice_iguard

  call Grid_getCellCoords(JAXIS,isg,CENTER,.true.,yi,GRID_JHI_GC)
  pdy = yi(2) - yi(1)
  cdy = 0.5e0*pdy
  
  do i = 1,yoend
     yo(i) = yi(offob+i)
  end do

  yi(1) = yo(1) + cdy - twice_iguard*pdy
  
  do i = 2,yiend
     yi(i) = yi(i-1)+2.e0*pdy
  end do
#endif
#if N_DIM == 3

!!$  offic = ka/2 + NGUARD - twice_iguard
  offic  = (ka-1+NGUARD)/ref_ratio - intpol_guard
  offoc = ka - 1
  
  zoend = kb-ka+1
!!$  ziend = zoend/ref_ratio + twice_iguard
  ziend = (kb-ka+ref_ratio)/ref_ratio + min(1,mod(kb-ka,ref_ratio)*mod(kb+NGUARD,ref_ratio)) + twice_iguard

  call Grid_getCellCoords(KAXIS,isg,CENTER,.true.,zi,GRID_KHI_GC)
  pdz = zi(2) - zi(1)
  cdz = 0.5e0*pdz
  
  do i = 1,zoend
     zo(i) = zi(offoc+i)
  end do

  zi(1) = zo(1) + cdz - twice_iguard*pdz
  
  do i = 2,ziend
     zi(i) = zi(i-1)+2.e0*pdz
  end do
#endif


  ! mask for variables that need to be input to the UMAPn routine
  toUmapMask = gcell_on_cc
  !!! mask for variables for which UMAPn's output is needed
  !!! (subset of variables given by toUmapMask)
  !!fromUmapMask(ii) = gcell_on_cc(ii)
  ! mask for variables for which the user really wants interpolation output
  ! (subset of variables given by fromUmapMask)
  toCallerMask = gcell_on_cc

#ifdef DENS_VAR
  if(niver.gt.0) then
     do ivar = 1,niver
        if (conserved_var .and. toCallerMask(ivar) &
             .and. (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS)) then
           toUmapMask(DENS_VAR) = .TRUE.
!!           if (gr_convertToConsvdForMeshCalls) fromUmapMask(DENS_VAR) = .TRUE.
           if (gr_convertToConsvdForMeshCalls) toCallerMask(DENS_VAR) = .TRUE.
        end if
     end do
  end if
#endif

  
  do ii = 1,n_dens
     if (  nud(ii).gt.0 ) then
        i_dens = iud(0,ii)
        if (conserved_var) then
           if (toUmapMask(i_dens)) then
              nu_d = nud(ii)
              do jj = 1,nu_d
                 ivar = iud(jj,ii)
                 toUmapMask(ivar) = .TRUE.
              end do
           end if
        else
           nu_d = nud(ii)
           need_relatives = .FALSE.
           do jj = 1,nu_d
              ivar = iud(jj,ii)
              if (toUmapMask(ivar)) then
                 need_relatives = .TRUE.
                 exit
              end if
           end do
           if (need_relatives) then
              do jj = 1,nu_d
                 ivar = iud(jj,ii)
                 toUmapMask(ivar) = .TRUE.
              end do
           end if
        end if
     end if
  end do


  numVarsPacked = 0
  do ivar = 1,niver
     if (toUmapMask(ivar)) then
        numVarsPacked = numVarsPacked + 1
        varIndexPacked(ivar) = numVarsPacked
        varIndexUnpacked(numVarsPacked) = ivar
     else
        varIndexPacked(ivar) = 0
     end if
  end do

  do i_dens = 1,n_dens
     if (nud(i_dens).gt.0 ) then
        if (varIndexPacked(iud(0,i_dens)).LE.0) then
           nud(i_dens) = 0
        else
           do ui = 1,nud(i_dens)
              if (varIndexPacked(iud(ui,i_dens)).LE.0) then
                 nud(i_dens) = 0
                 exit
              end if
           end do
        end if
     end if
     if (  nud(i_dens).gt.0 ) then
        iud(0,i_dens) = varIndexPacked(iud(0,i_dens))
        do ui = 1,nud(i_dens)
           iud(ui,i_dens) = varIndexPacked(iud(ui,i_dens))
        end do
     end if
  end do

#ifdef DENS_VAR
  densIndexPacked = varIndexPacked(DENS_VAR)
#endif


  if (numVarsPacked == 0) then
     ! Nothing to do, RETURN prematurely!
     return
  end if


  inxu = numVarsPacked*xiend
  onxu = numVarsPacked*xoend  
#ifdef DEBUG_GRID

  ! optional check of array dimensioning

  if ( size(pu,1) < inxu ) then
     write(*,*) '[amr_prolong_gen_unk1_fun] ERROR: pu size too small along dim 1: ',size(pu,1),' < ',numVarsPacked*xiend
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: pu too small along dim 1")
  end if
  if ( size(pu,2) < yiend ) then
     write(*,*) '[amr_prolong_gen_unk1_fun] ERROR: pu size too small along dim 2: ',size(pu,2),' < ',yiend
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: pu too small along dim 2")
  end if
  if ( size(pu,3) < ziend ) then
     write(*,*) '[amr_prolong_gen_unk1_fun] ERROR: pu size too small along dim 3: ',size(pu,3),' < ',ziend
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: pu too small along dim 3")
  end if

  if ( size(qu,1) < onxu ) then
     write(*,*) '[amr_prolong_gen_unk1_fun] ERROR: qu size too small along dim 1: ',size(qu,1),' < ',numVarsPacked*xoend
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: qu too small along dim 1")
  end if
  if ( size(qu,2) < yoend ) then
     write(*,*) '[amr_prolong_gen_unk1_fun] ERROR: qu size too small along dim 2: ',size(qu,2),' < ',yoend
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: qu too small along dim 2")
  end if
  if ( size(qu,3) < zoend ) then
     write(*,*) '[amr_prolong_gen_unk1_fun] ERROR: qu size too small along dim 3: ',size(qu,3),' < ',zoend
     call Driver_abortFlash("[amr_prolong_gen_unk1_fun] ERROR: qu too small along dim 3")
  end if
#endif

  pu = 0.e0

  do k = 1,ziend
     k1 = (offic+koff)*K3D+k
     do j = 1,yiend
        j1 = (offib+joff)*K2D+j


        do ivarPacked = 1,numVarsPacked
           ivar = varIndexUnpacked(ivarPacked)
           if (toUmapMask(ivar)) then
              do i = 1,xiend
                 i1 = offia+ioff+i
#ifdef DENS_VAR
                 if (gr_convertToConsvdInMeshInterp .and. (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS)) then
                    if (recv(DENS_VAR, i1, j1, k1) .eq. 0.0) then
                       if (recv(ivar,i1,j1,k1) .ne. 0.0) then
                          ! This situation would probably lead to division by zero errors in the
                          ! unk1(ivar)/unk1(dens) operation when converting back from conserved form later,
                          ! if we did no check. Abort if recv(ivar)!=0 and recv(dens)==0, but let
                          ! recv(ivar)==recv(dens)==0 pass. - KW
99                        format ('[amr_prolong_gen_unk1_fun] PE=',I7,', ivar=',I3,', value=',1PG9.3)
                          print 99,mype,ivar,recv(ivar,i1,j1,k1)
                          print*,'Trying to convert non-zero mass-specific variable to per-volume form in monotonic&
                               & mesh interpolation, but dens is zero!'
                          call Driver_abortFlash &
                               ('Trying to convert non-zero mass-specific variable to per-volume form, but dens is zero!')
                       end if
                    end if
                    pu((ivarPacked-1)*xiend+i,j,k) = recv(ivar,i1,j1,k1) * recv(DENS_VAR, i1, j1, k1)
                 else
                    pu((ivarPacked-1)*xiend+i,j,k) = recv(ivar,i1,j1,k1)
                 end if
#else
                 pu((ivarPacked-1)*xiend+i,j,k) = recv(ivar,i1,j1,k1)
#endif
              end do
           end if
        end do

     end do
  end do

  imapcx = 1

  if(NDIM == 1) then
     call umap1 (xiend,xi,pdx,inxu,pu,&
          &      xoend,xo,cdx,onxu,qu,&
          &      numVarsPacked,gr_intpol,imapcx,&  
          &      geom(1),ref_ratio,&
          &      wmap,lw,iwmap,liw,ierr,conserved_var,gr_smallx)
  elseif(NDIM==2) then
     call umap2 (lxiend,lyiend,xiend,yiend,linxu,inxu,xi,pdx,yi,pdy,pu,&
          &      lxoend,lyoend,xoend,yoend,lonxu,onxu,xo,cdx,yo,cdy,qu,&
          &      numVarsPacked,gr_intpol,imapcx,&  
          &      geom(1),geom(2),ref_ratio,ref_ratio,&
          &      wmap,lw,iwmap,liw,ierr,conserved_var,gr_smallx)
  else
     call umap3 (lxiend,lyiend,lziend,xiend,yiend,ziend,&
          &      linxu,inxu,xi,pdx,yi,pdy,zi,pdz,pu,&
          &      lxoend,lyoend,lzoend,xoend,yoend,zoend,&
          &      lonxu,onxu,xo,cdx,yo,cdy,zo,cdz,qu,&
          &      numVarsPacked,gr_intpol,imapcx,&
          &      geom(1),geom(2),geom(3),ref_ratio,ref_ratio,ref_ratio,&
          &      wmap,lw,iwmap,liw,ierr,conserved_var,gr_smallx)
  end if

  ! Check whether the interpolation routines detected an error condition,
  ! and abort if that is the case. - KW
  if (ierr .ne. 0) then
     print*,'**ERROR** - ierr from UMAPn is ', ierr, ', isg=', isg
     call Driver_abortFlash('ERROR - ierr from UMAPn is not 0')
  end if

  do k = 1,zoend
     k1 = offoc*K3D+k
     do j = 1,yoend
        j1 = offob*K2D+j
        do ivarPacked = 1,numVarsPacked
           ivar = varIndexUnpacked(ivarPacked)
           if (toCallerMask(ivar)) then

              do i = 1,xoend
                 i1 = offoa+i
#ifdef DENS_VAR
                 if (gr_convertToConsvdInMeshInterp .and. (gr_vartypes(ivar) .eq. VARTYPE_PER_MASS)) then
                    if (qu((densIndexPacked-1)*xoend+i,j,k) == 0.0) then
                       ! If unk(dens)==0, assume that for all ivars of interest --- namely,
                       ! the PER_MASS type ones --- unk(ivar)==0 held before the forward
                       ! conversion to conserved form (and thus should be made true again
                       ! here after interpolation); otherwise, the program should have
                       ! aborted above in the forward conversion. - KW
                       unk1(ivar,i1,j1,k1,idest) = 0.0
                    else
                       dens_inv = 1. / qu((densIndexPacked-1)*xoend+i,j,k)
                       unk1(ivar,i1,j1,k1,idest) = dens_inv * qu((ivarPacked-1)*xoend+i,j,k)
                    end if
                 else
                    unk1(ivar,i1,j1,k1,idest) = qu((ivarPacked-1)*xoend+i,j,k)
                 end if
#else
                 unk1(ivar,i1,j1,k1,idest) = qu((ivarPacked-1)*xoend+i,j,k)
#endif
              end do
           end if
        end do
     end do
  end do
  
  return
end subroutine amr_prolong_gen_unk1_fun
