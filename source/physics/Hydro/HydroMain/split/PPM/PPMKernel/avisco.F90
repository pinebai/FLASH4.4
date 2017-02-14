!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/avisco
!!
!! NAME
!!
!!  avisco
!!
!!
!! SYNOPSIS
!!
!!  avisco(integer(IN) :: j, 
!!         integer(IN) :: k,
!!         real(OUT)   :: avis(*),
!!         integer(IN) :: dirGeom(MDIM),
!!         integer(IN) :: xyzswp,
!!         integer(IN) :: nzni,
!!         integer(IN) :: nznf,
!!         integer(IN) :: nsdim,
!!         real(IN)    :: xtop,
!!         real(IN)    :: xbot,
!!         real(IN)    :: ytop,
!!         real(IN)    :: ybot,
!!         real(IN)    :: ylft,
!!         real(IN)    :: yrgt, 
!!         real(IN)    :: zlft,
!!         real(IN)    :: zrgt,
!!         real(IN)    :: x(*),
!!         real(IN)    :: xl(*),
!!         real(IN)    :: xzn(*),
!!         real(IN)    :: yzn(*),
!!         real(IN)    :: zzn(*),
!!         real(IN)    :: u(*),
!!         real(IN)    :: uttp(*),
!!         real(IN)    :: utbt(*),
!!         real(IN)    :: utrt(*),
!!         real(IN)    :: utlt(*),
!!         real(IN)    :: cvisc)
!!
!!
!! DESCRIPTION
!!
!!  This routine computes the artifical viscosity that should be 
!!  applied inside shocks.  It is computed correctly for all the 
!!  geometries we support.  DEV: Verify this claim!
!!
!! ARGUMENTS
!!
!!  j -       index that indicates where we are in terms of the
!!            first transversal coordinate
!!  k -       index that indicates where we are in terms of the
!!            second transversal coordinate
!!  avis -    Output. The artifical viscosity, for use by hydro_1d.
!!
!!  xyzswp -  direction of current sweep, one of SWEEP_X, SWEEP_Y, or SWEEP_Z
!!  nzni -    initial index for the loop in sweep direction
!!  nznf -    final index   for the loop in sweep direction
!!  nsdim -   number of dimensions, should be NDIM
!!  x -   positions of cell centers of the 1d array (in sweep direction)
!!  xl -   positions of left interfaces of the 1d array (sweep direction)
!!  xzn,yzn,zzn -  IAXIS, JAXIS, and KAXIS coordinates of each cell in the 1d array
!!  dirGeom - geometry types of the coordinate directions - probably as
!!            set in Hydro_init.
!!  xtop -
!!  xbot - 
!!  ytop - 
!!  ybot - 
!!  ylft - 
!!  yrgt - 
!!  zlft - 
!!  zrgt - 
!!
!!  u - 
!!  uttp -
!!  utbt -
!!  utlt - 
!!  utrt -
!!  cvisc -   Artificial viscosity constant. Probably as set
!!            in Hydro_init from a runtime parameter.
!!
!!***

subroutine avisco(j, k, avis, dirGeom, xyzswp,                 &
                  nzni, nznf, nsdim,                                          &
                  xtop, xbot, ytop, ybot, ylft, yrgt, zlft, zrgt, &
                  x, xl, xzn, yzn, zzn,                                       &
                  u, uttp, utbt, utrt, utlt, cvisc  )
  use Driver_interface, ONLY : Driver_abortFlash

  implicit NONE
#include "constants.h"

  integer,intent(IN)  :: j, k 
  integer,intent(IN) :: dirGeom(MDIM)
  integer,intent(IN)  :: xyzswp
  real,intent(OUT)     :: avis(*)
  integer,intent(IN)  :: nzni, nznf, nsdim
  real,intent(IN)     :: xtop, xbot, ytop, ybot
!! The following also used to be passed in, but were actually unused.
!! I therefore removed them from the argument list. KW 
!  real     ::                         ztop, zbot

  real,intent(IN)     ::             ylft, yrgt, zlft, zrgt
  real,intent(IN)     :: x(*), xl(*), xzn(*), yzn(*), zzn(*)
  real,intent(IN)     :: u(*), uttp(*), utbt(*), utrt(*), utlt(*)
  real,intent(IN)     :: cvisc
      
  real     :: dxtb, dytb, dztb
  real     :: sintop, sinbot, sinth, dyrl, dzrl
  real     :: dx3, xlst, sthdth
  integer  :: i
!
!------------------------------
  if (nsdim == 1)   then
!------------------------------
!

     if (dirGeom(IAXIS) == XYZ) then

        do i = nzni, nznf
           avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))
           avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
        enddo

     end if

     if (dirGeom(IAXIS) == RAD_CYL) then

        do i = nzni, nznf
           dxtb = x(i)**2 - x(i-1)**2
           if ( dxtb /= 0.e0 ) then
              avis(i) = (x(i) * u(i) - x(i-1) * u(i-1)) * 2.e0 / dxtb
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           end if
        enddo

     end if

     if (dirGeom(IAXIS) == RAD_SPH) then

        do i = nzni, nznf
           avis(i) = (x(i)**2 * u(i) - x(i-1)**2 * u(i-1)) * 3.0 /         &
                     (x(i)**3 - x(i-1)**3) 
           avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
        enddo

     end if
!
!------------------------------
  end if
!------------------------------
!
!                                    
!------------------------------
  if (nsdim == 2)   then
!------------------------------
!
!
!-------  (x, y)
!

     if (dirGeom(IAXIS) == XYZ  .and.  dirGeom(JAXIS) == XYZ) then    

        if (xyzswp == SWEEP_X) then
           dytb = 0.5 / (ytop - ybot)
           do i = nzni, nznf
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +               &
                   (uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1)) *           &
                   dytb
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        else
           dxtb = 0.5 / (xtop - xbot)
           do i = nzni, nznf
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +               &
                   (uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1)) *           &
                   dxtb
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        end if
!
!
!-------  (r_cyl, z)
!
     else if (dirGeom(IAXIS) == RAD_CYL  .and.  dirGeom(JAXIS) == XYZ) then 
!
        if (xyzswp == SWEEP_X) then 

           dytb = 0.5e0 / (ytop - ybot)

           do i = nzni, nznf
              dxtb = x(i)**2 - x(i-1)**2
              if ( dxtb /= 0.e0 ) then
                 avis(i) = (x(i) * u(i) - x(i-1) * u(i-1)) * 2.e0 / dxtb   &
                          +(uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1))     &
                          *dytb
                 avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
              else
                 avis(i) = 0.0     ! Should not happen?
              end if
           enddo

        else

           do i = nzni, nznf
              dxtb = xtop**2 - xbot**2
              if ( dxtb /= 0.e0 ) then
                 avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))               &
                          +( xtop * (uttp(i) + uttp(i-1))                  &
                            -xbot * (utbt(i) + utbt(i-1))) / dxtb
                 avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
              else
                 avis(i) = 0.0     ! Should not happen?
              end if
           enddo

        end if
!
!
!-------  (r_cyl, phi)
!
     else if (dirGeom(IAXIS) == RAD_CYL  .and.  dirGeom(JAXIS) == PHI_CYL) then 

        if (xyzswp == SWEEP_X) then

           dytb = 0.5e0 / (ytop - ybot)

           do i = nzni, nznf
              dxtb = x(i)**2 - x(i-1)**2
              avis(i) = 0.0
              if ( dxtb /= 0.e0 )                                          &
              avis(i) = (x(i) * u(i) - x(i-1) * u(i-1)) * 2.e0 / dxtb
              if ( xl(i) /= 0.e0 )                                         &
              avis(i) = avis(i)                                            &
                       +(uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1))        &
                       *dytb/xl(i)
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo

        else

           dxtb = xtop**2 - xbot**2

           if ( dxtb /= 0.e0 ) then

              dxtb = 1.e0 / dxtb

              do i = nzni, nznf
                 avis(i) = (u(i) - u(i-1)) / (xzn(j)*(x(i) - x(i-1)))      &
                          +( xtop * (uttp(i) + uttp(i-1))                  &
                            -xbot * (utbt(i) + utbt(i-1)) ) * dxtb
                 avis(i) = - cvisc * avis(i) * xzn(j) * (x(i) - x(i-1))
              enddo

           else
              avis(nzni:nznf) = 0.0     ! Should not happen?
           end if

        end if
!
!
!-------  (r_sph, theta)
!
     else if (dirGeom(IAXIS) == RAD_SPH  .and.  dirGeom(JAXIS) == THETA) then 
!
        if (xyzswp == SWEEP_X) then

           sintop = sin(ytop)
           sinbot = sin(ybot)
           sinth  = sin(yzn(j))
           dytb   = 2.e0 * sinth * (ytop - ybot)

           do i = nzni,nznf
              dx3 = x(i)**3 - x(i-1)**3

              avis(i) = 0.0

              if ( dx3 /= 0.e0 ) &
                 avis(i) = (x(i)**2 * u(i) - x(i-1)**2 * u(i-1)) * 3.e0 / dx3

              xlst = xl(i) * dytb

              if ( xlst /= 0.e0 )                                          &
                 avis(i) = avis(i) +                                       &
                   (sintop * (uttp(i) + uttp(i-1)) -                       &
                    sinbot * (utbt(i) + utbt(i-1))  )  / xlst

              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           end do
        else

           dxtb = 3.e0 / ( 2.e0 * (xtop**3 - xbot**3) )

           do i = nzni,nznf
              sthdth = sin(xl(i)) * (x(i) - x(i-1))

              avis(i) = 0.0

              if ( sthdth /= 0.e0 ) &
                 avis(i) = (sin(x(i)) * u(i) - sin(x(i-1)) * u(i-1))       &
                          /(xzn(j) * sthdth)

              avis(i) = avis(i) +                                          &
                   (xtop**2 * (uttp(i) + uttp(i-1)) -                      &
                    xbot**2 * (utbt(i) + utbt(i-1))  ) * dxtb

              avis(i) = - cvisc * avis(i) * xzn(j) * (x(i) - x(i-1))
           end do
        end if
!
!
!-------  (r_sph, phi) ;  theta = pi/2
!
     else if (dirGeom(IAXIS) == RAD_SPH  .and.  dirGeom(JAXIS) == PHI_SPH) then  
!
        if (xyzswp == SWEEP_X) then
           dytb = 0.5 / (ytop - ybot)
           do i = nzni, nznf
              avis(i) = (x(i)**2 * u(i) - x(i-1)**2 * u(i-1)) * 3.0 /      &
                   (x(i)**3 - x(i-1)**3)  +                                &
                   (uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1)) /           &
                   xl(i) * dytb
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        else
           do i = nzni, nznf
              dxtb = 3.0 / ( 2.0 * (xtop**3 - xbot**3) )
              avis(i) = (u(i) - u(i-1)) / (xzn(j)*(x(i) - x(i-1))) +       &
                   (xtop**2 * (uttp(i) + uttp(i-1)) -                      &
                   xbot**2 * (utbt(i) + utbt(i-1))   ) * dxtb
              avis(i) = - cvisc * avis(i) * xzn(j) * (x(i) - x(i-1))
           enddo
        end if
!
!
     else 
 
        print *, ' Error: geometry not implemented in avisco()'
        call Driver_abortFlash("Error: geometry not implemented in avisco()")
!
     end if
!
!------------------------------
  end if
!------------------------------
!
!                                    
!------------------------------
  if (nsdim == 3)   then
!------------------------------
!
!
!-------  (x, y, z)
!
     if (dirGeom(IAXIS)==XYZ .and. dirGeom(JAXIS)==XYZ .and. dirGeom(KAXIS)==XYZ) then 
!
        if (xyzswp == SWEEP_X) then
           dytb = 0.5 / (ytop - ybot)
           dzrl = 0.5 / (zrgt - zlft)
           do i = nzni, nznf
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +               &
                   (uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1)) *           &
                   dytb +                                                  &
                   (utrt(i) + utrt(i-1) - utlt(i) - utlt(i-1)) *           &
                   dzrl
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        else if (xyzswp == SWEEP_Y) then
           dxtb = 0.5 / (xtop - xbot)
           dzrl = 0.5 / (zrgt - zlft)
           do i = nzni, nznf
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +               &
                   (uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1)) *           &
                   dxtb +                                                  &
                   (utrt(i) + utrt(i-1) - utlt(i) - utlt(i-1)) *           &
                   dzrl
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        else
           dxtb = 0.5 / (xtop - xbot)
           dyrl = 0.5 / (yrgt - ylft)
           do i = nzni, nznf
              avis(i) = (u(i) - u(i-1)) / (x(i) - x(i-1))  +               &
                   (uttp(i) + uttp(i-1) - utbt(i) - utbt(i-1)) *           &
                   dxtb +                                                  &
                   (utrt(i) + utrt(i-1) - utlt(i) - utlt(i-1)) *           &
                   dyrl
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        end if
!
!
!-------  (r_sph, theta, phi)
!
     else if (dirGeom(IAXIS)==RAD_SPH .and. &
          dirGeom(JAXIS)==THETA .and. dirGeom(KAXIS)==PHI_SPH) then 
!
        if (xyzswp == SWEEP_X) then
           sintop = sin(ytop)
           sinbot = sin(ybot)
           sinth  = sin(yzn(j))
           dytb   = 0.5 / ( sinth * (ytop - ybot) )
           dzrl   = 0.5 / ( sinth * (zrgt - zlft) )
           do i = nzni, nznf
              avis(i) = (x(i)**2 * u(i) - x(i-1)**2 * u(i-1)) * 3.0 /      &
                   (x(i)**3 - x(i-1)**3)  +                                &
                   (sintop * (uttp(i) + uttp(i-1)) -                       &
                   sinbot * (utbt(i) + utbt(i-1))  )  /                    &
                   xl(i) * dytb  +                                         &
                   (utrt(i) + utrt(i-1) - utlt(i) - utlt(i-1)) /           &
                   xl(i) * dzrl                                  
              avis(i) = - cvisc * avis(i) * (x(i) - x(i-1))
           enddo
        else if (xyzswp == SWEEP_Y)  then
           dxtb  = 3.0 / ( 2.0 * (xtop**3 - xbot**3) )
           dzrl  = 0.5 / ( xzn(j) * (zrgt - zlft) ) 
           do i = nzni, nznf
              sinth   = max(sin( xl(i) ),1e-6)
              avis(i) = (sin(x(i)) * u(i) - sin(x(i-1)) * u(i-1)) /        &
                   (xzn(j) * sinth * (x(i) - x(i-1)) ) +                   &
                   (xtop**2 * (uttp(i) + uttp(i-1)) -                      &
                   xbot**2 * (utbt(i) + utbt(i-1))  ) * dxtb +             &
                   (utrt(i) + utrt(i-1) - utlt(i) - utlt(i-1)) /           &
                   sinth * dzrl
              avis(i) = - cvisc * avis(i) * xzn(j) * (x(i) - x(i-1))
           enddo
        else
           sintop = sin(yrgt)
           sinbot = sin(ylft)
           sinth  = sin(yzn(k))
           dxtb   = 3.0 / ( 2.0 * (xtop**3 - xbot**3) )
           dyrl   = 0.5 / ( xzn(j) * sinth * (yrgt - ylft) )
           do i = nzni, nznf
              avis(i) = (u(i) - u(i-1)) /                                  &
                   (xzn(j) * sinth * (x(i) - x(i-1)) ) +                   &
                   (xtop**2 * (uttp(i) + uttp(i-1)) -                      &
                   xbot**2 * (utbt(i) + utbt(i-1))  ) * dxtb +             &
                   (sintop  * (utrt(i) + utrt(i-1)) -                      &
                   sinbot  * (utlt(i) + utlt(i-1))  ) * dyrl  
              avis(i) = - cvisc * avis(i) * xzn(j) * sinth *               &
                   (x(i) - x(i-1))
           enddo
        end if
!
!
     else
!
        print *, '[AVISCO] ERROR: geometry not implemented.'
        call Driver_abortFlash("[AVISCO] ERROR: geometry not implemented.")
!
     end if
!
!------------------------------
  end if
!------------------------------
!
!
!------- No Viscosity if zone has positive divergence  or
!                     if zone is outside shock.
!        Reflecting boundaries are handled by the caller.
!

  do i = nzni, nznf
     avis(i) = max (avis(i),  0.e0)
  enddo

  return
end subroutine avisco
