!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/multiTemp/hy_dbgWriteStates
!!
!! NAME
!!
!!  hy_dbgWriteStates
!!
!! SYNOPSIS
!!
!!  call hy_dbgWriteStates(character(len=*),intent(IN)  :: vname,
!!                         real,intent(IN),dimension(numCells)  :: var,
!!                         real,intent(IN),dimension(numCells)  :: varl,
!!                         real,intent(IN),dimension(numCells)  :: varr,
!!                         real,intent(IN),dimension(numCells)  :: varav,
!!                         real,intent(IN),dimension(numCells)  :: varflx,
!!                         real,intent(IN),dimension(numCells)  :: rho,
!!                         real,intent(IN),dimension(numCells)  :: rhol,
!!                         real,intent(IN),dimension(numCells)  :: rhor,
!!                         real,intent(IN),dimension(numCells)  :: rhoav,
!!                         real,intent(IN),dimension(numCells)  :: urell,
!!                         real,intent(IN)  :: dt,
!!                         integer,intent(IN)  :: numintcells,
!!                         integer,intent(IN)  :: numcells,
!!                         integer,intent(IN)  :: nguard,
!!                         real,intent(IN),dimension(numCells)  :: xcoords,
!!                         real,intent(IN),dimension(numCells)  :: xl,
!!                         real,intent(IN),dimension(numCells)  :: xr,
!!                         integer,intent(IN)  :: blockid)
!!
!! DESCRIPTION
!!
!! outputs states
!!
!! ARGUMENTS
!!
!!   vname : var name
!!
!!   var : variable
!!
!!   varl : left state
!!
!!   varr : right state
!!
!!   varav : avg. variable
!!
!!   varflx : variable flux 
!!
!!   rho : density
!!
!!   rhol : left state of density
!!
!!   rhor : right state of density
!!
!!   rhoav : avg. density
!!
!!   urell : 
!!
!!   dt : timestep
!!
!!   numintcells : number of internal cells 
!!
!!   numcells : number of cells
!!
!!   nguard : number of guard cells
!!
!!   xcoords : x coordinates
!!
!!   xl : left x
!!
!!   xr : right x
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!***

#include "Flash.h"

#ifdef SELE_MSCALAR
#define SELE_IND (SELE_MSCALAR-SPECIES_BEGIN+1)
#else
#define SELE_IND 0
#endif

subroutine hy_dbgWriteStates(vname,var,varl,varr, &
          varav, varflx,&
          rho,rhol,rhor, rhoav,&
          urell,dt,&
          numIntCells, numCells, nguard,&
          xCoords,xL,xR,blockID)
  use Driver_interface, ONLY: Driver_getNStep, Driver_getSimTime
  use Hydro_data,   ONLY:  hy_numXn, hy_smlrho, hy_smallx, hy_small, hy_gravl, &
                           hy_drho, hy_rho6, &
                           hy_du, hy_u6, hy_dut, hy_ut6, &
                           hy_dutt, hy_utt6, hy_dp, hy_dgame, hy_game6, &
                           hy_gravr, hy_p6, hy_dgamc, hy_gamc6, hy_dgrav, &
                           hy_grav6, hy_dxn, hy_xn6, &
                           hy_deint, hy_eint6, hy_detot,hy_etot6,hy_dpFrac,hy_pFrac6, &
                           hy_pwcubic, hy_pwl, hy_pwr, hy_dpw, hy_pw6l, hy_pw6r, &
                           hy_ppmModifystates, &
                           hy_useSteepening, hy_useCmaFlattening,&
                           hy_epsiln, hy_omg1, hy_omg2, hy_charLimiting, &
                           hy_dbgReconstConsvSele
  implicit none
#undef REAL_FORMAT
#define REAL_FORMAT "(ES20.13)"

  character(len=*),intent(IN) :: vname
  integer,intent(IN) :: numIntCells, numCells, nguard
  real,intent(IN),dimension(numCells) :: var,varl,varr
  real,intent(IN),dimension(numCells) :: varav,varflx
  real,intent(IN),dimension(numCells) :: rho,rhol,rhor,rhoav
  real,intent(IN),dimension(numCells) :: xCoords,xL,xR,urell
  real,intent(IN) :: dt
  integer,intent(IN) :: blockID

  real,dimension(numCells) :: dvar,var6
  character(len=80) :: ff1
  character(len=10) :: aindic
  integer,parameter :: nStart=5, nEnd=10
  integer :: nCurr, nstep, i, sGen
  real :: sTime
  integer,parameter :: nInner = 64, nInnerp = nInner+1
  integer,parameter :: fu=28
  real :: x,xm,xp, shiftL, shiftR, avBox, flxBox, flxRaw
  real :: dxInner, alph, fRec, fvarRec, frhoRec
  logical :: perVolVar

  if (blockID.GT.6) return
  call Driver_getNStep(nstep)
!!$  if (nstep < 16 .OR. nstep > 41) then
!!$     return
!!$  end if
  if (nstep < 0 .OR. nstep > 60) then
     return
  end if
  call Driver_getSimTime(sTime,sGen)
  select case (vname)
  case('rho')
     dvar(:) = hy_drho(1:numCells)
     var6(:) = hy_rho6(1:numCells)
  case('eint')
     dvar(:) = hy_deint(1:numCells,0)
     var6(:) = hy_eint6(1:numCells,0)
  case('eion')
     dvar(:) = hy_deint(1:numCells,1)
     var6(:) = hy_eint6(1:numCells,1)
  case('eele')
     dvar(:) = hy_deint(1:numCells,2)
     var6(:) = hy_eint6(1:numCells,2)
  case('etot')
     dvar(:) = hy_detot(1:numCells,0)
     var6(:) = hy_etot6(1:numCells,0)
  case('sele')
#if SELE_IND > 0
     dvar(:) = hy_dxn(1:numCells,SELE_IND)
     var6(:) = hy_xn6(1:numCells,SELE_IND)
#else
     dvar(:) = 0.0
     var6(:) = 0.0
#endif
  end select
  if (vname .EQ. 'rho') then
     perVolVar = .TRUE.
  else  if (hy_dbgReconstConsvSele .AND. vname .EQ. 'sele') then
     perVolVar = .TRUE.
  else  if (vname .EQ. 'sele') then
     perVolVar = .false.
  else
     perVolVar = .FALSE.
  end if

22 format('step',I4.4,'.',I1,'-',A,'.txt')
  write(ff1,22) nstep,sGen,trim(vname)
  open(fu,file=ff1,form='formatted',status='UNKNOWN')

  nCurr = nStart
  do while (nCurr .LE. nEnd)
!     print*,'nCurr,x:',nCurr,xCoords(nCurr)
     dxInner = (xR(nCurr)-xL(nCurr))/nInner
     if (urell(nCurr)>0.0) then
        aindic = '    < VVVV'
     else if (urell(nCurr)<0.0) then
        aindic = '    < ^^^^'
     else 
        aindic = '    < ----'
     end if
!     write(fu, 9998)'   ',xl(nCurr),varav(nCurr),aindic
     shiftL = min(urell(nCurr)*dt,0.0)
     shiftR = max(urell(nCurr+1)*dt,0.0)
     do i=0,nInner
        x = xL(nCurr) + i*dxInner
        alph = i*dxInner/(xR(nCurr)-xL(nCurr))
        fvarRec = varl(nCurr)+alph*(dvar(nCurr) + var6(nCurr)*(1-alph))
        frhoRec = rhol(nCurr)+alph*(hy_drho(nCurr) + hy_rho6(nCurr)*(1-alph))
        if (perVolVar) then
           fRec = fvarRec
        else
           fRec = fvarRec * frhoRec
        end if
        flxBox = 0.0
        if (x.LE.xL(nCurr)-shiftL) then
           avBox = varav(nCurr)
           if (.NOT.perVolVar) avBox = avBox * rhoav(nCurr)
           if (urell(nCurr) .NE. 0.0) flxBox = varflx(nCurr) / urell(nCurr)
        else if (x.GE.xR(nCurr)-shiftR) then
           avBox = varav(nCurr+1)
           if (.NOT.perVolVar) avBox = avBox * rhoav(nCurr+1)
           if (urell(nCurr+1) .NE. 0.0) flxBox = varflx(nCurr+1) / urell(nCurr+1)
        else
           avBox = 0.0
        end if
        flxRaw = 0.0
        if (x==xL(nCurr)) flxRaw = varflx(nCurr)
!        print 9999,'   ',x,fRec,varl(nCurr),var(nCurr),varr(nCurr)
        if (perVolVar) then
           write(fu, 9999)'   ',x,fRec,var(nCurr),avBox,flxBox, flxRaw
        else
           write(fu, 9999)'   ',x,fRec,var(nCurr)*rho(nCurr),avBox,flxBox, flxRaw
        end if
9999    format(a,F20.8,5(G20.12))
9998    format(a,F20.8,1(G20.12),' #',A)
     end do
     write(fu,*)
     nCurr = nCurr+1
  end do
  if (urell(nCurr)>0.0) then
     aindic = '    < VVV.'
  else if (urell(nCurr)<0.0) then
     aindic = '    < ^^^.'
  else 
     aindic = '    < ---.'
  end if
!  write(fu, 9998)'   ',xl(nCurr),varav(nCurr),aindic

  close(fu,status='KEEP')
end subroutine hy_dbgWriteStates
  
