#define CIP_
!!****if* source/hydro/explicit/split/ppm/hydro_1d
!!
!! NAME
!!
!!  ModuleHydro_1d
!!
!!
!! SYNOPSIS
!!  
!!  call hydro_1d(...)
!!
!!
!! DESCRIPTION
!!
!!  Compute the 1-d directionally split fluxes through the boundary
!!  of the computational zone using the PPM method.   
!!
!!***

module ModuleHydro_1d
  
  use runtime_parameters, ONLY: &
       get_parm_from_context, global_parm_context

!!Not sure what the access method is for these.. dBaseDeclarations
  use dBase, ONLY: & 
       ionmax,  bnd_reflect,                &
       nguard, nxb, nyb, nzb, k2d, k3d,     &
       ndim, numFluxes => nfluxes,          &
       i_strt => iLo_gc, i_end => iHi_gc,   &
       j_strt => jLo_gc, j_end => jHi_gc,   &
       k_strt => kLo_gc, k_end => kHi_gc,   &
       mfaces,                              &
       sweep_x,sweep_y,sweep_z,             &   
       dBaseGetPtrToXCoords,                &
       dBaseGetPtrToYCoords,                &
       dBaseGetPtrToZCoords,                &
       dBasePropertyReal,                   &
       dBaseNeighborBlock,                  &
       dBaseKeyNumber, dBasePropertyInteger

  use Gravity, ONLY: GravAccelOneRow

  use PPMModule, ONLY: PPMInit, &
                       intrfc,  &
                       states,  &
                       rieman                       

  use hydroModule, ONLY:   use_cma_advection
  
  implicit none

  private
  public :: hydro_1d

contains
  subroutine hydro_1d (j, k, xyzswp, block_no, q, qn, dt, dx, dtdx, &
       & u, ut, utt, rho, p, e, tmp, game, gamc, xn, grav, &
       & nzn, &
       & x, xl, xr, y, z, radial_coord, &
       & xbot, xtop, ybot, ytop, ylft, yrgt, zlft, zrgt, ugrid, &
       & utbt, uttp, utlt, utrt, &
       & igeom, lgrav, &
       & shock_multid, &
       & tempArea, tempGrav1d_o, tempGrav1d, tempDtDx, tempFict, &
       & temp_flx, temp_fly, temp_flz)


    implicit none

!--arguments-------------------------
    integer ::  j, k, xyzswp, block_no, q, qn, nzn, igeom,  &
                ntmpvar, iarea, igrav1d, idtdx

    logical :: lgrav

    real :: dt

    real, DIMENSION(i_strt:i_end, &
                    j_strt:j_end, &
                    k_strt:k_end  ) :: tempArea, tempGrav1d_o, tempGrav1d, &
                                       tempDtDx, tempFict

    real, DIMENSION(numFluxes,i_strt:i_end,j_strt:j_end,k_strt:k_end) :: &
         & temp_flx, temp_fly, temp_flz


    real, DIMENSION(q)    :: rho, u, ut, utt, p, e, ei, tmp
    real, DIMENSION(q,qn) :: xn, xnav, xnflx, xnl, xnr
    real, DIMENSION(q)    :: rhoav, uav, utav, uttav, pav, rhoflx, &
         &                   uflx, utflx, uttflx, eflx, eintflx, eint, &
         &                   rhol, rhor, ul, ur, utl, utr, uttl, uttr, &
         &                   pl, pr, vl, vr, gamcl, gamcr, c, ce, ugrid, &
         &                   urell, ugrdl, game, gamc, gameav, gamel, gamer, &
         &                   uttp, utbt, utlt, utrt, v, dvol, x, xr, xl, y, z, &
         &                   radial_coord, &
         &                   grav, ograv, hgrav, ngrav, fict

    real, DIMENSION(q) :: shock_multid

    real, DIMENSION(q) ::  scrch1
    real, DIMENSION(q) ::  avis
    real, DIMENSION(q) ::  xzn, yzn, zzn
    real, pointer, dimension(:,:,:) :: xCoord, yCoord, zCoord

!--locals-------------------------
    integer :: i, n, kk
    integer :: nzn4, nzn5, nzn8

    real, DIMENSION(q) :: xlold, xrold, dxold, dvold, &
         &                 alold, aold, &
         &                 area, areal, arear, dx, dtdx


    real ::  sum, suminv, dtfac, dg

    real :: xbot, xtop, ybot, ytop, zbot, ztop, xlft, xrgt, &
         &       ylft, yrgt, zlft, zrgt

    real, save :: smallu, small, smallp, smlrho, smallx, & 
         rieman_tol, omg1, omg2, epsiln, dp_sh, cvisc

    integer, save :: igodu, nriem, ishkbn, igeomx, igeomy, igeomz

    real :: dtold

    integer, DIMENSION(mfaces) :: neigh

    integer, save :: izn, myproc

    integer, save :: irhoflx, iuflx, ipflx, iutflx, iuttflx, &
         ieflx, ieintflx, inucflx_begin
    integer, save :: ixCoord, iyCoord, izCoord, &
                     lowerFace, upperFace

    real :: pres_jump
    
    logical, save :: hybrid_riemann

    real, save :: vgrid
    logical, save :: moving_grid

    logical, save :: firstcall = .true.
#ifdef CIP

    integer, save :: itrcr, inuc_begin
    real, save    :: pi
#endif

!BGL add arrays for vector math
    integer :: npts
    real :: recrho(q), recdx(q), tmpsqrt(q), recdvol(q)


    if (firstcall) then

       irhoflx  = dBaseKeyNumber('rhoflx')
       iuflx    = dBaseKeyNumber('uflx')
       ipflx    = dBaseKeyNumber('pflx')
       iutflx   = dBaseKeyNumber('utflx')
       iuttflx  = dBaseKeyNumber('uttflx')
       ieflx    = dBaseKeyNumber('eflx')
       ieintflx = dBaseKeyNumber('eintflx')

       if (irhoflx < 0) call abort_flash("[HYDRO_1D] ERROR: no such key rhoflx")
       if (iuflx < 0)   call abort_flash("[HYDRO_1D] ERROR: no such key uflx")
       if (ipflx < 0)   call abort_flash("[HYDRO_1D] ERROR: no such key pflx")
       if (iutflx < 0)  call abort_flash("[HYDRO_1D] ERROR: no such key utflx")
       if (iuttflx < 0) call abort_flash("[HYDRO_1D] ERROR: no such key uttflx")
       if (ieflx < 0)   call abort_flash("[HYDRO_1D] ERROR: no such key eflx")
       if (ieintflx < 0)call abort_flash("[HYDRO_1D] ERROR: no such key eintflx")

       inucflx_begin = dBaseKeyNumber('nucflx_begin')

       myproc = dBasePropertyInteger("MyProcessor")

       call get_parm_from_context("smallu", smallu)
       call get_parm_from_context("small", small)
       call get_parm_from_context("smlrho", smlrho)
       call get_parm_from_context("smallp", smallp)
       call get_parm_from_context("smallx", smallx)

       call get_parm_from_context("nriem", nriem)
       call get_parm_from_context("rieman_tol", rieman_tol)

       call get_parm_from_context("ishkbn", ishkbn)

       call get_parm_from_context("igodu", igodu)
       call get_parm_from_context("epsiln", epsiln)
       call get_parm_from_context("omg1", omg1)
       call get_parm_from_context("omg2", omg2)
       call get_parm_from_context("dp_sh", dp_sh)
       call get_parm_from_context("cvisc", cvisc)

       call get_parm_from_context("igeomx", igeomx)
       call get_parm_from_context("igeomy", igeomy)
       call get_parm_from_context("igeomz", igeomz)

       call get_parm_from_context("hybrid_riemann", hybrid_riemann)

       izn      = dBaseKeyNumber('zn')

       ixCoord   = dBaseKeyNumber("xCoord")
       iyCoord   = dBaseKeyNumber("yCoord")
       izCoord   = dBaseKeyNumber("zCoord")

       lowerFace = dBaseKeyNumber("lowerFace")
       upperFace = dBaseKeyNumber("upperFace")

       ograv(:) = 0.e0       ! initialize grav arrays to zero
       hgrav(:) = 0.e0       ! need to reconsider whether this
       ngrav(:) = 0.e0       ! is the best place to do this

       call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'vgrid', vgrid)

       if (vgrid /= 0.e0) then
          moving_grid = .true.
       else
          moving_grid = .false.
       endif

    endif

    dtold = dBasePropertyReal("OldTimeStep")

    nzn4 = nzn + 4; nzn5 = nzn + 5; nzn8 = nzn + 8

    xCoord  => dBaseGetPtrToXCoords()
    yCoord  => dBaseGetPtrToYCoords()
    zCoord  => dBaseGetPtrToZCoords()
!!    neigh(:) = dBaseNeighborBlockList(block_no)

    !---------------------------
    call geom (j, k, xyzswp, block_no, igeom, nzn, q, &
         &     areal, arear, area, dx, dvol, xl, xr)
    !---------------------------

    if (lgrav) then
       call GravAccelOneRow (j, k, xyzswp, block_no, 0, ograv, nzn8)
       call GravAccelOneRow (j, k, xyzswp, block_no, 1, grav, nzn8)

       dtfac = dt/dtold

       do i = 1,nzn8
          dg       = dtfac*(grav(i) - ograv(i))
          hgrav(i) = grav(i) + 0.5e0*dg
          ngrav(i) = grav(i) +       dg
       enddo

    else

       ograv(:) = 0.e0
       hgrav(:) = 0.e0
       ngrav(:) = 0.e0

    end if

    !---------------------------
    call force (j, k, block_no, q, nzn, igeom, x, u, ut, utt, fict)
    !---------------------------

    do i = 2, nzn8
       ugrdl(i) = 0.5e0 * (ugrid(i) + ugrid(i-1))
    end do
    ugrdl(1) = ugrdl(2)

    call vrec(recrho, rho, nzn8)
    do i = 1, nzn8
      tmpsqrt(i) = gamc(i) * p(i) * rho(i)
    end do
    call vsqrt(tmpsqrt, tmpsqrt, nzn8)

    do i = 1, nzn8
       v(i)  = recrho(i)
       c(i)  = tmpsqrt(i)
       ce(i) = c(i)*v(i)
    end do


    ! Initialization -- grid/geometry factors

    call vrec(recdx, dx, nzn8)
    do i = 1, nzn8
       dtdx(i) = dt * recdx(i)
    enddo

    ! Note - if igeom = 3, 4, or 5, the x coordinate will be a radial
    ! coordinate and the index j will always refer to the x direction.

    if (igeom >= 3) then
       do i = 1, nzn8
          dtdx(i) = dtdx(i) / xCoord(izn,j,block_no)
       enddo
    endif

!------------------------------------------------------------------------------

    if ( firstcall ) then

       call PPMInit(small, smallu, smallp, igodu, smlrho, smallx, &
                    dp_sh, nriem, rieman_tol)

       firstcall = .false.

    end if

!             Compute zone-edge, time-averaged fluxes using PPM.

! Obtain PPM interpolation coefficients.

    call intrfc(nzn,rho,u,ut,utt,p, &
         &      rhol,rhor,ul,ur, &
         &      utl,utr,uttl,uttr, &
         &      pl,pr,vl,vr,gamcl, &
         &      gamcr,game, &
         &      gamel,gamer,gamc,hgrav,xn, &
         &      xnl,xnr,v,dx, x, ishkbn, &
         &      epsiln, omg1, omg2)       

    ! Determine input states for the Riemann solver.

    call states (j,igeom,nzn,rho,u,rhol,rhor,ul,ur,utl,utr,uttl,uttr, &
         &       p, pl,pr,gamcl,gamcr,ugrid,ce,game,gamer,gamc,gamel,xnl,xnr, &
         &       dtdx, dt, x, xl, radial_coord, hgrav, fict)

    ! Solve Riemann problems at zone edges.

    call rieman(nzn,ei,rhoav,uav,utav,uttav,pav,urell,ugrdl,game, &
         &      gameav,xnav,x,xyzswp,block_no,myproc,j,k)

    ! no, I *really* meant reflecting

    select case (xyzswp)

    case (sweep_x)
       if (dBaseNeighborBlock(block_no,ixCoord,lowerFace) .eq. bnd_reflect) then
          uav  (nguard+1) = 0.e0
          urell(nguard+1) = 0.e0
       endif
       if (dBaseNeighborBlock(block_no,ixCoord,upperFace) .eq. bnd_reflect) then
          uav  (nguard+nxb+1) = 0.e0
          urell(nguard+nxb+1) = 0.e0
       endif

    case (sweep_y)
       if (dBaseNeighborBlock(block_no,iyCoord,lowerFace) .eq. bnd_reflect) then
          uav  (nguard*k2d+1) = 0.e0
          urell(nguard*k2d+1) = 0.e0
       endif
       if (dBaseNeighborBlock(block_no,iyCoord,upperFace) .eq. bnd_reflect) then
          uav  (nguard*k2d+nyb+1) = 0.e0
          urell(nguard*k2d+nyb+1) = 0.e0
       endif

    case (sweep_z)
       if (dBaseNeighborBlock(block_no,izCoord,lowerFace) .eq. bnd_reflect) then
          uav  (nguard*k3d+1) = 0.e0
          urell(nguard*k3d+1) = 0.e0
       endif
       if (dBaseNeighborBlock(block_no,izCoord,upperFace) .eq. bnd_reflect) then
          uav  (nguard*k3d+nzb+1) = 0.e0
          urell(nguard*k3d+nzb+1) = 0.e0
       endif
    end select

    ! -------------------------------------------------------------------
    ! Consistent Multi-fluid Advection (Plewa & Mueller 1999, CMA Eq. 13)

    ! mass fractions renormalization and optional limiting (not if CMA)

    if ( ionmax > 1  ) then

       scrch1(5:nzn5) = 0.e0

       do n = 1, ionmax

          ! renormalize and limit mass fractions: note that limiting introduces
          ! non-conservation of species

          if ( .not.use_cma_advection ) then
             do i = 5, nzn5
                xnav(i,n) = max(smallx, min(1.e0, xnav(i,n)))
             end do
          end if

          do i = 5, nzn5
             scrch1(i) = scrch1(i) + xnav(i,n)
          end do
       end do

       do i = 5, nzn5
          if ( scrch1(i) /= 0.e0 ) scrch1(i) = 1.e0 / scrch1(i)
       end do

       do n = 1, ionmax
          do i = 5, nzn5
             xnav(i,n) = xnav(i,n) * scrch1(i)
          end do
       end do

    end if

    !--------------------------------------------------------------------------
    ! Save old grid information and move the grid using the previously computed 
    ! grid velocities (for artificial dissipation purposes).

    do i = 1, nzn8
       xlold(i) = xl(i)
       xrold(i) = xr(i)
       dxold(i) = dx(i)
       dvold(i) = dvol(i)
       alold(i) = areal(i)
       aold(i)  = area(i)
    enddo

    if (moving_grid) then
       do i = 2, nzn8
          xl(i)   = xlold(i) + dt * ugrdl(i)
          xr(i-1) = xl(i)
       enddo
       
       do i = 2, nzn8
          x(i)  = 0.5e0 * (xr(i) + xl(i))    
          dx(i) = xr(i) - xl(i)
       enddo
    endif

    call geom (j, k, xyzswp, block_no, igeom, nzn, q, &
         &     areal, arear, area, dx, dvol, xl, xr)

    !------------------------------------------------------------------------------
    ! Compute unmodified fluxes for each of the conserved quantities.

    do i = 5, nzn5
       rhoflx(i) = rhoav (i) * urell(i)
       uflx(i)   = rhoflx(i) * uav  (i)
       utflx(i)  = rhoflx(i) * utav (i)
       uttflx(i) = rhoflx(i) * uttav(i)
       scrch1(i) = pav(i) / ( rhoav(i) * (gameav(i)-1.e0) )

       ! compute the internal energy flux
       eintflx(i) = rhoflx(i) * scrch1(i) + uav(i) * pav(i)

       ! add the kinetic energy 
       scrch1(i) = scrch1(i) &
            &        + 0.5e0 * (uav(i)**2 + utav(i)**2 + uttav(i)**2)

       ! compute the total energy flux
       eflx(i)   = rhoflx(i) * scrch1(i) + uav(i) * pav(i)

       ! initialize the artificial viscosity coefficient
       avis(i) = 0.e0
    enddo

    ! compute the internal energy over the whole structure
    call vrec(recrho, rho, nzn8)
    do i = 1,nzn8 
       eint(i) = e(i) - 0.5e0*(u(i)**2 + ut(i)**2 + utt(i)**2)
       eint(i) = max(eint(i),smallp*recrho(i))
    enddo

    !------------------------------------------------------------------------------
    ! Compute quantities needed for artificial viscosity.  Unfortunately in 
    ! multidimensions this requires some direction-dependent code.

    if (cvisc > 0.e0)  then 

       if (xyzswp == sweep_x)  then
          xzn(:) = x(:)
          yzn(:) = yCoord(izn, j, block_no)
          zzn(:) = zCoord(izn, k, block_no)
       endif
       if (xyzswp == sweep_y)  then
          xzn(:) = xCoord(izn, j, block_no)
          yzn(:) = x(:)
          zzn(:) = zCoord(izn, k, block_no)
       endif
       if (xyzswp == sweep_z)  then
          xzn(:) = xCoord(izn, j, block_no)
          yzn(:) = yCoord(izn, k, block_no)
          zzn(:) = x(:)
       endif

       call avisco( j, k , avis,                              &
            igeomx, igeomy, igeomz, xyzswp, 5, nzn5, ndim,    &
            xtop, xbot, ytop, ybot, ylft, yrgt, ztop, zbot,   &
            zlft, zrgt,                                       &
            x, xl, xzn, yzn, zzn,                             &
            u, uttp, utbt, utrt, utlt, cvisc  )
       
       ! there should be no flux through a reflecting boundary, so force the 
       ! artificial viscosity there to 0

       if (xyzswp == sweep_x) then
          if (dBaseNeighborBlock(block_no,ixCoord,lowerFace) == &
               bnd_reflect) avis(nguard+1) = 0.e0
          if (dBaseNeighborBlock(block_no,ixCoord,upperFace) == &
               bnd_reflect) avis(nguard+nxb+1) = 0.e0
       endif

       if (xyzswp == sweep_y) then
          if (dBaseNeighborBlock(block_no,iyCoord,lowerFace) == &
               bnd_reflect) avis(nguard*k2d+1) = 0.e0
          if (dBaseNeighborBlock(block_no,iyCoord,upperFace) == &
               bnd_reflect) avis(nguard*k2d+nyb+1) = 0.e0
       endif
       
       if (xyzswp == sweep_z) then
          if (dBaseNeighborBlock(block_no,izCoord,lowerFace) == &
               bnd_reflect) avis(nguard*k3d+1) = 0.e0
          if (dBaseNeighborBlock(block_no,izCoord,upperFace) == &
               bnd_reflect) avis(nguard*k3d+nzb+1) = 0.e0
       endif

       avis(5:nzn5) = avis(5:nzn5)*alold(5:nzn5)
    endif

    ! odd-even decoupling fix
    !
    ! now that we've called the accurate Riemann solver, loop over all the zones,
    ! and for any that is inside a shock, call the HLLE Riemann solver, and 
    ! replace the fluxes computed above with those.  This fixes the odd-even
    ! decoupling problem (see Quirk 1997)

    if (hybrid_riemann) then

       do i = 5, nzn5

          ! check for the presence of shocks by looking at the artificial viscosity
          ! and the pressure jump

          pres_jump = abs(p(i) - p(i-1))/min(p(i),p(i-1))

!          if (pres_jump >= dp_sh .AND. avis(i) /= 0.e0) then
          if (shock_multid(i) == 1 .AND. avis(i) /= 0.e0) then

             ! the interface between zones i and i-1 is a shock.  Use the HLLE solver
             ! and replace the fluxes computed above

             call riemann_hlle( &
                  rho(i-1), rho(i), &
                  u(i-1), u(i), &
                  ut(i-1), ut(i), &
                  utt(i-1), utt(i), &
                  p(i-1), p(i), &
                  e(i-1), e(i), &
                  eint(i-1), eint(i), &
                  game(i-1), game(i), &
                  gamc(i-1), gamc(i), &
                  hgrav(i-1), hgrav(i), &
                  dt, &
                  rhoflx(i), uflx(i), utflx(i), uttflx(i), eflx(i), eintflx(i))
             
          endif
          
       enddo

    endif

    ! Do the elemental abundance fluxes.

    do n = 1, qn
       do i = 5, nzn5
          xnflx(i,n) = xnav(i,n) * rhoflx(i)
       enddo
    enddo

    !---------------------------------------------------------------
    ! Apply the diffusive fluxes to the unmodified conserved fluxes.

    do i = 5, nzn5
       rhoflx(i) = rhoflx(i) * alold(i)
       uflx  (i) = uflx  (i) * alold(i)
       utflx (i) = utflx (i) * alold(i)
       uttflx(i) = uttflx(i) * alold(i)
       eflx  (i) = eflx  (i) * alold(i)

       eintflx(i) = eintflx(i)*alold(i)

       rhoflx(i)  = rhoflx(i)  + avis(i) * (rho(i-1)           - rho(i))  
       uflx  (i)  = uflx  (i)  + avis(i) * (rho(i-1)*u(i-1)    - rho(i)*u(i))
       utflx (i)  = utflx (i)  + avis(i) * (rho(i-1)*ut(i-1)   - rho(i)*ut(i))
       uttflx(i)  = uttflx(i)  + avis(i) * (rho(i-1)*utt(i-1)  - rho(i)*utt(i))
       eflx  (i)  = eflx  (i)  + avis(i) * (rho(i-1)*e(i-1)    - rho(i)*e(i))
       eintflx(i) = eintflx(i) + avis(i) * (rho(i-1)*eint(i-1) - rho(i)*eint(i))
    enddo

    do n = 1, qn
       do i = 5, nzn5
          xnflx(i,n) = xnflx(i,n) * alold(i)
          xnflx(i,n) = xnflx(i,n) + &
                       avis(i) * (rho(i-1)*xn(i-1,n) - rho(i)*xn(i,n))
       enddo
    enddo

    !------------------------------------------------------------------------------
    ! Store dt/dx, geometry factors, and the modified fluxes in 'global' arrays
    ! for use in updating the solution after all of the 1D sweeps in this direction
    ! are done. Note that dt/dx is not constant in non-Cartesian geometries, so it
    ! is saved for each zone in the tempDtDx() array.

    npts = nzn4 - 5 + 1
    call vrec(recdvol(5), dvol(5), npts)
    do i = 5, nzn4
       dtdx(i) = dt*recdvol(i)
    enddo

    select case (xyzswp)

    case (sweep_x)

       do i = 5, nzn4
          tempDtDx(i,j,k)     = dtdx(i)
          tempArea(i,j,k)     = area(i)
          tempGrav1d_o(i,j,k) = grav(i)
          tempGrav1d(i,j,k)   = ngrav(i)
          tempFict(i,j,k)     = fict(i)
       enddo

       do i = 5, nzn4 + 1
          temp_flx(irhoflx,i,j,k)  = rhoflx(i)
          temp_flx(iuflx,i,j,k)    = uflx(i)
          temp_flx(ipflx,i,j,k)    = pav(i)
          temp_flx(iutflx,i,j,k)   = utflx(i)
          temp_flx(iuttflx,i,j,k)  = uttflx(i)
          temp_flx(ieflx,i,j,k)    = eflx(i)
          temp_flx(ieintflx,i,j,k) = eintflx(i)
       enddo
#ifdef CIP

       ! recover proper flux for CIP

       if ( itrcr > 0 ) then

          do n = itrcr,itrcr
             kk = n - inuc_begin + 1
             do i = 5, nzn5, nzn5-5
                xnav(i,kk) = atan(xnav(i,kk))/(0.9999d0*pi) + 0.5e0
                xnflx(i,kk) = xnav(i,kk) * rhoflx(i)
             end do
          end do

       end if
#endif

       do kk = 1, qn
          do i = 5, nzn4 + 1
             temp_flx(inucflx_begin+kk-1,i,j,k) = xnflx(i,kk)
          end do
       enddo

    case (sweep_y)

       do i = 5, nzn4
          tempDtDx(j,i,k)     = dtdx(i)
          tempArea(j,i,k)     = area(i)
          tempGrav1d_o(j,i,k) = grav(i)
          tempGrav1d(j,i,k)   = ngrav(i)
          tempFict(j,i,k)     = fict(i)
       enddo

       do i = 5, nzn4 + 1
          temp_fly(irhoflx,j,i,k)  = rhoflx(i)
          temp_fly(iuflx,j,i,k)    = uflx(i)
          temp_fly(ipflx,j,i,k)    = pav(i)
          temp_fly(iutflx,j,i,k)   = utflx(i)
          temp_fly(iuttflx,j,i,k)  = uttflx(i)
          temp_fly(ieflx,j,i,k)    = eflx(i)
          temp_fly(ieintflx,j,i,k) = eintflx(i)
       enddo
#ifdef CIP

       ! recover proper flux for CIP

       if ( itrcr > 0 ) then

          do n = itrcr,itrcr
             kk = n - inuc_begin + 1
             do i = 5, nzn5, nzn5-5
                xnav(i,kk) = atan(xnav(i,kk))/(0.9999d0*pi) + 0.5e0
                xnflx(i,kk) = xnav(i,kk) * rhoflx(i)
             end do
          end do

       end if
#endif

       do kk = 1, qn
          do i = 5, nzn4 + 1
             temp_fly(inucflx_begin+kk-1,j,i,k) = xnflx(i,kk)
          end do
       enddo

    case (sweep_z)

       do i = 5, nzn4
          tempDtDx(j,k,i)     = dtdx(i)
          tempArea(j,k,i)     = area(i)
          tempGrav1d_o(j,k,i) = grav(i)
          tempGrav1d(j,k,i)   = ngrav(i)
          tempFict(j,k,i)     = fict(i)
       enddo

       do i = 5, nzn4 + 1
          temp_flz(irhoflx,j,k,i)  = rhoflx(i)
          temp_flz(iuflx,j,k,i)    = uflx(i)
          temp_flz(ipflx,j,k,i)    = pav(i)
          temp_flz(iutflx,j,k,i)   = utflx(i)
          temp_flz(iuttflx,j,k,i)  = uttflx(i)
          temp_flz(ieflx,j,k,i)    = eflx(i)
          temp_flz(ieintflx,j,k,i) = eintflx(i)
       enddo
#ifdef CIP

       ! recover proper flux for CIP

       if ( itrcr > 0 ) then

          do n = itrcr,itrcr
             kk = n - inuc_begin + 1
             do i = 5, nzn5, nzn5-5
                xnav(i,kk) = atan(xnav(i,kk))/(0.9999d0*pi) + 0.5e0
                xnflx(i,kk) = xnav(i,kk) * rhoflx(i)
             end do
          end do

       end if
#endif

       do kk = 1, qn
          do i = 5, nzn4 + 1
             temp_flz(inucflx_begin+kk-1,j,k,i) = xnflx(i,kk)
          end do
       enddo

    end select

!===================================================================
    
    return
  end subroutine hydro_1d
  
end module ModuleHydro_1d
