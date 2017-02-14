 


        subroutine mapcx
     &               ( nx,x,hdx,mx,xx,hdxx,
     &                 ix,dxm,dx,dxp,dxpp,
     &                 intp,
     &                 igx,igrx,ier )
*
*
*       SUBROUTINE MAPCX
*       -----------------
*
*       set x grid pointers for xx grid and interpolation scale terms
*       for linear or quadratic interpolation
*
*       x  - input grid
*       xx - output grid
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        REAL
     &  a2_b2, a, b, a3_b3, xr, xhl, xhr, xl, dvz, dvl, dvr,  
     &  cxl, cxr, cxhl, cxhr, dvrli

        INTEGER
     &  mx,nx,i,ii,isrt,igx,igrx,intp,ier

        INTEGER
     &  ix(mx)

        REAL
     &  x(nx),hdx,xx(mx),hdxx,dxm(mx),dx(mx),dxp(mx),dxpp(mx)

        a2_b2 ( a,b ) = ( a - b )*( a + b )
        a3_b3 ( a,b ) = ( a - b )*( a*a + a*b + b*b )
*
*
*
*       ----------------------------------------------------------
*       find right interface xr=0.5(x(i)+x(i+1)) s.t. xr.ge.xx(ii)
*       i.e. point xx(ii) belongs to zone "i"
        if ( max(nx,mx).eq.1 ) then
           ix(1) = 1
        else
           isrt = 1
           do ii = 1,mx
              ix(ii) = 0
              do i = isrt,nx-1
                 xr = x(i) + hdx
                 if ( xr.ge.xx(ii) ) then
                    isrt = i
                    ix(ii) = i
                    go to 1
                 end if
              end do
1             continue
              if ( ix(ii).eq.0 .and. isrt.eq.nx-1 ) then
                 if ( x(nx)+hdx.ge.xx(ii) ) ix(ii) = nx
              end if
           end do
        end if
*       ---------------------
*       set geometrical terms
        if ( intp.gt.0 ) then
*          -------------------------
*          for averaging over volume
           do ii = 1,mx
              i = ix(ii)
              xhl = x(i)   - hdx
              xhr = x(i)   + hdx
              xl  = xx(ii) - hdxx
              xr  = xx(ii) + hdxx
              xl = max( xl,xhl )
              xr = min( xr,xhr )
              dxm (ii) = xhl
              dx  (ii) = xl
              dxp (ii) = xr
              dxpp(ii) = xhr
           end do
           if ( igx.eq.0 ) then
*             -------------
*             cartesian (x)
              do ii = 1,mx
                 dvz = dxpp(ii) - dxm(ii)
                 dvl = dx(ii)   - dxm(ii)
                 dvr = dxp(ii)  - dxm(ii)
                 dxm(ii) = dvl
                 dxp(ii) = dvr
                 dx(ii)  = dvz
              end do
           else if ( igx.eq.1 ) then
*             ----------------------
*             cylindrical radial (r)
              do ii = 1,mx
                 dvz = a2_b2( dxpp(ii) , dxm(ii) )
                 dvl = a2_b2( dx(ii)   , dxm(ii) )
                 dvr = a2_b2( dxp(ii)  , dxm(ii) )
                 dxm(ii) = dvl
                 dxp(ii) = dvr
                 dx(ii)  = dvz
              end do
           else if ( igx.eq.2 ) then
*             --------------------
*             spherical radial (r)
              do ii = 1,mx
                 dvz = a3_b3( dxpp(ii) , dxm(ii) )
                 dvl = a3_b3( dx(ii)   , dxm(ii) )
                 dvr = a3_b3( dxp(ii)  , dxm(ii) )
                 dxm(ii) = dvl
                 dxp(ii) = dvr
                 dx(ii)  = dvz
              end do
           else if ( igx.eq.3 ) then
*             -----------------------
*             cylindrical angle (phi)
              do ii = 1,mx
                 dvz = dxpp(ii) - dxm(ii)
                 dvl = dx(ii)   - dxm(ii)
                 dvr = dxp(ii)  - dxm(ii)
                 dxm(ii) = dvl
                 dxp(ii) = dvr
                 dx(ii)  = dvz
              end do
           else if ( igx.eq.4 ) then
*             -----------------------------
*             spherical polar angle (theta)
              do ii = 1,mx
                 cxl  = cos(dx(ii))
                 cxr  = cos(dxp(ii))
                 cxhl = cos(dxm(ii))
                 cxhr = cos(dxpp(ii))
                 dvz = cxhl - cxhr
                 dvl = cxhl - cxl
                 dvr = cxhl - cxr
                 dxm(ii) = dvl
                 dxp(ii) = dvr
                 dx(ii)  = dvz
              end do
           else if ( igx.eq.5 ) then
*             -------------------------------
*             spherical azimuthal angle (phi)
              do ii = 1,mx
                 dvz = dxpp(ii) - dxm(ii)
                 dvl = dx(ii)   - dxm(ii)
                 dvr = dxp(ii)  - dxm(ii)
                 dxm(ii) = dvl
                 dxp(ii) = dvr
                 dx(ii)  = dvz
              end do
           end if
              
           do ii = 1,mx
              dvl = dxm(ii)
              dvr = dxp(ii)
              dvz = dx(ii)
              dvrli = one/(dvr - dvl)
              dxm (ii) = dvl / dvz ! dxl
              dx  (ii) = dvr / dvz ! dxr
              dxp (ii) = dvl * dvrli ! dvfl
              dxpp(ii) = dvr * dvrli ! dvfr
           end do
        end if
        return
        end
        
*
* --------------- UMAP1.F ---------------
*
 
        subroutine umap1
     &(
     &    nx,x ,hdx ,nxu,p1u,
     &    mx,xx,hdxx,mxu,q1u,
     &    nui,intp,imapcx,
     &    igx,igrx,
     &    w,lw,iw,liw,ier,conserved_var,smallx
     &)
*
*
*       SUBROUTINE UMAP1
*       ----------------
*
*       interpolates
*       p(i), x(i), i=1,...,nx
*       onto
*       q(ii), xx(ii), ii=1,...,mx
*
*       see UMAP3 for more information
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

        INTEGER
     &  i1, i2, i3, i4, ier, igrx, igx,
     &  ii, ii1, ii2, ii3, ii4, iiu, imapcx,
     &  intp, inu, iui, liw, lw, mx, mxu, nui, nx, nxu, iw

        REAL
     &  dxl, dxr, hdx, hdxx, smallx

        REAL
     &          p1u,q1u,
     &          x,xx,
     &          w

        DIMENSION
     &          p1u(nxu),q1u(mxu),
     &          x(nx),xx(mx),
     &          w(lw),iw(liw)
*
*
*
        ier = 0
*       -------------
*       working space
        i1 = 1
        i2 = i1 + mx
        i3 = i2 + mx
        i4 = i3 + mx
*       -------------------
*       geometrical factors
        if ( imapcx.eq.1 )
     &        call mapcx
     &       ( nx,x,hdx,mx,xx,hdxx,
     &         iw,w(i1),w(i2),w(i3),w(i4),
     &         intp,
     &         igx,igrx,ier )

*       ----------------------
*       interpolation pointers
        inu = 0
        iiu = 0
        do iui = 1,nui
           do ii = 1,mx
              iiu = (iui-1)*mx + ii
              ip_i_iiu (iiu) = (iui-1)*nx + iw(ii)
              ip_ii_iiu(iiu) = ii
           end do
        end do
*       ---------------
*       partial volumes
        if ( intp.eq.1 ) then
           do ii = 1,mx
              ii1 = i1 + ii - 1
              ii2 = i2 + ii - 1
              ii3 = i3 + ii - 1
              ii4 = i4 + ii - 1
              dxl         = w(ii1)
              dxr         = w(ii2)
              dvfl_ii(ii) = w(ii3)
              dvfr_ii(ii) = w(ii4)
              hdl_ii(ii)  = half*dxl
              hdr_ii(ii)  = half*dxr
           end do
        else if ( intp.eq.2 ) then
           do ii = 1,mx
              ii1 = i1 + ii - 1
              ii2 = i2 + ii - 1
              ii3 = i3 + ii - 1
              ii4 = i4 + ii - 1
              dxl         = w(ii1)
              dxr         = w(ii2)
              dvfl_ii(ii) = w(ii3)
              dvfr_ii(ii) = w(ii4)
              hdl_ii(ii)  = half*dxl
              hdr_ii(ii)  = half*dxr
              hd6l_ii(ii) = one - c2by3*dxl
              hd6r_ii(ii) = one - c2by3*dxr
           end do
        end if
        call map1
     &            ( nxu,p1u,mxu,q1u,
     &              mx,intp,igx,igrx,1,ier,
     &              conserved_var,smallx )
*
* return
*
        return
*
* end of UMAP1
*
        end
 


        subroutine map1
     &                 ( nxu,p1u,mxu,q1u,
     &                   mx,intp,igx,igrx,isw,ier,
     &                   conserved_var,smallx )
*
*
*       SUBROUTINE MAP1
*       ---------------
*
*       interpolate p1u(nx,x) onto q1u(mx,xx)
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

	INTEGER
     &  ier, isw, igrx, igx, intp, mx, nxu, mxu, iiu, i, ii

	REAL
     &  dupl, dupr, dubs, dup, pl, pr, sigi p1u3,
     &  dpmapc_i, p6mapc_i, dp_p6, dp_dp, hd_dv_l, hd_dv_r, sigi, p1u3,
     &  smallx

        REAL
     &          p1u, q1u
        DIMENSION
     &          p1u(nxu), q1u(mxu)
*
*
*
        if ( intp.eq.0 .or. igrx.eq.1 ) then
*          -----------------------
*          no interpolation needed
*          -----------------------
           do iiu = 1,mxu
              i   = ip_i_iiu(iiu)
              q1u(iiu) = p1u(i)
           end do
        else if ( intp.eq.1 ) then
*          --------------------
*          linear interpolation
*          --------------------
           do i = 2,nxu-1
              dupl = p1u(i  )-p1u(i-1)
              dupr = p1u(i+1)-p1u(i  )
              dubs = 1.e-3*(abs(dupl)+abs(dupr)+1.e-2)**2
              dup  = half*( dupl+dupr )
     &              *( (two*dupl*dupr + dubs)
     &                /(dupl*dupl + dupr*dupr + dubs) )
              dup = sign( min(abs(dup),
     &                    min(abs(dupl),abs(dupr))),
     &                    dup )
              if ( dupl*dupr.le.zero ) dup = zero
              pl = p1u(i) - half*dup
              pr = p1u(i) + half*dup
              pl = max( min(p1u(i-1),p1u(i)),
     &             min( max(p1u(i-1),p1u(i)), pl ))
              pr = max( min(p1u(i+1),p1u(i)),
     &             min( max(p1u(i+1),p1u(i)), pr ))
              plmapc(i) = pl
              dpmapc(i) = pr - pl
           end do
*          --------------
*          volume average
           do iiu = 1,mxu
              i   = ip_i_iiu(iiu)
              ii  = ip_ii_iiu(iiu)
              q1u(iiu) =
     &                   plmapc(i)
     &                  *(dvfr_ii(ii) - dvfl_ii(ii))
     &                  +dpmapc(i)
     &                  *(hdr_ii(ii)*dvfr_ii(ii)-hdl_ii(ii)*dvfl_ii(ii))
           end do
           if ( n_dens.gt.0 )
     &     call map1_f2d ( q1u,mxu,mx,conserved_var,smallx )
        else if ( intp.eq.2 ) then
*          -----------------------
*          quadratic interpolation
*          -----------------------
*          zone to zone differences
           do i = 2,nxu
              wvmap1(i) = p1u(i) - p1u(i-1)
              wvmap2(i) = abs(wvmap1(i) + wvmap1(i))
              wvmap3(i) = sign(one , wvmap1(i))
           end do
!          ! fixed to start at 2 - KW
           do i = 2,nxu-1
*             parabola mean slope
              wvmap4(i) = half*(wvmap1(i+1) + wvmap1(i))
*             slope limiting : CW84 (1.8)
              if ( wvmap4(i).lt.zero ) then
                 sigi = -one
              else
                 sigi = +one
              end if
              wvmap4(i) = min(abs(wvmap4(i)),wvmap2(i),wvmap2(i+1))
     &                   *sigi
              if ( -wvmap3(i)*wvmap3(i+1).ge.zero ) wvmap4(i) = zero
           end do
*          just do it
           wvmap4(nxu) = zero
*          interface values : CW84 (1.6)
!          ! fixed loop range - KW
           do i = 2,nxu-2
!             ! prmapc(i) = a_{i+1/2}
              !    = a_i + 0.5 (a_{i+1}-a_i) -1/6 delta_m a_{i+1} + 1/6 delta_m a_i
              prmapc(i)   =   p1u(i)
     &                      + half  * wvmap1(i+1)
     &                      - c1by6 * wvmap4(i+1)
     &                      + c1by6 * wvmap4(i)
              plmapc(i+1) = prmapc(i)
           end do
           do i = 3,nxu-2       ! fixed loop range - KW
              
*             monotonization : CW84 (1.10)
              if ((prmapc(i)-p1u(i))*(p1u(i)-plmapc(i)).le.zero) then
                 plmapc(i) = p1u(i)
                 prmapc(i) = p1u(i)
              end if
           end do
           do i = 3,nxu-2       ! fixed loop range - KW
              p1u3      = three*p1u(i)
              wvmap1(i) = p1u3 - two*prmapc(i)
              wvmap2(i) = p1u3 - two*plmapc(i)
           end do
           do i = 3,nxu-2       ! fixed loop range - KW
*             mean slope and quadratic coefficient
              dpmapc_i = prmapc(i) - plmapc(i)
              p6mapc_i = six*( p1u(i) - half*(plmapc(i)+prmapc(i)) )
*             monotonize quadratic coefficient
              dp_p6 = dpmapc_i*p6mapc_i
              dp_dp = dpmapc_i*dpmapc_i
              if ( +dp_p6 - dp_dp.gt.zero ) then
                 plmapc(i) = wvmap1(i)
              else
     &              if ( -dp_dp - dp_p6.gt.zero ) then
                 prmapc(i) = wvmap2(i)
              end if
           end do
*          final profile
           do i = 3,nxu-2       ! fixed loop range - KW
!             ! dpmapc(i) = a_{R,i}-a_{L_i}
              dpmapc(i) = prmapc(i) - plmapc(i)
              p6mapc(i) = six*( p1u(i) - half*(plmapc(i)+prmapc(i)) )
           end do
*          --------------
*          volume average
           do iiu = 1,mxu
              i   = ip_i_iiu(iiu)
              ii  = ip_ii_iiu(iiu)
              hd_dv_l = hdl_ii(ii)*dvfl_ii(ii)
              hd_dv_r = hdr_ii(ii)*dvfr_ii(ii)
              q1u(iiu) =
     &                   plmapc(i)
     &                  *(dvfr_ii(ii) - dvfl_ii(ii))
     &                  +dpmapc(i)
     &                  *(hd_dv_r - hd_dv_l)
     &                  +p6mapc(i)
     &                  *(hd6r_ii(ii)*hd_dv_r - hd6l_ii(ii)*hd_dv_l)
           end do
           if ( n_dens.gt.0 )
     &     call map1_f2d ( q1u,mxu,mx,conserved_var,smallx )
        else
*          ------------
*          unknown intp
*          ------------
           ier = 11
        end if
*
* return
*
        return
*
* end of MAP1
*
        end



        subroutine map1_f2d ( q1u,mxu,mx,conserved_var,smallx )
*
*
*       SUBROUTINE MAP1_F2D
*       -------------------
*
*       calculate total density from partial densities or renormalize
*       mass fractions (if .not.conserved_var)
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

	INTEGER
     &  mx, mxu, i_dens, ixdi, ixdf, ixd, ju, iu, ixui, ixuf, ixu

        REAL
     &  smallx, smallxdp

        REAL
     &          q1u
        DIMENSION
     &          q1u(mxu)
*
*
*
        smallxdp = smallx

        do i_dens = 1,n_dens
           if ( nud(i_dens).gt.0 ) then
              if ( conserved_var ) then

                 ixdi = (iud(0,i_dens)-1)*mx + 1
                 ixdf = ixdi + mx - 1
                 do ixd = ixdi,ixdf
                    wvmap1(ixd) = zero
                 end do
                 do ju = 1,nud(i_dens)
                    iu = iud(ju,i_dens)
                    ixui = (iu-1)*mx + 1
                    ixuf = ixui + mx - 1
                    do ixu = ixui,ixuf
                       ixd = ixu - ixui + ixdi
                       wvmap1(ixd) = wvmap1(ixd) + q1u(ixu)
                    end do
                 end do
                 do ixd = ixdi,ixdf
                    q1u(ixd) = wvmap1(ixd)
                 end do

              else

                 ixdi = (iud(0,i_dens)-1)*mx + 1
                 ixdf = ixdi + mx - 1
                 do ixd = ixdi,ixdf
                    wvmap1(ixd) = zero
                 end do
                 do ju = 1,nud(i_dens)
                    iu = iud(ju,i_dens)
                    ixui = (iu-1)*mx + 1
                    ixuf = ixui + mx - 1
                    do ixu = ixui,ixuf
                       ixd = ixu - ixui + ixdi
                       wvmap1(ixd) = wvmap1(ixd)
     &                              +max(smallxdp,min(one, q1u(ixu) ))
                    end do
                 end do
                 do ixd = ixdi,ixdf
                    wvmap1(ixd) = 1.e0 / wvmap1(ixd)
                 end do
                 do ju = 1,nud(i_dens)
                    iu = iud(ju,i_dens)
                    ixui = (iu-1)*mx + 1
                    ixuf = ixui + mx - 1
                    do ixu = ixui,ixuf
                       ixd = ixu - ixui + ixdi
                       q1u(ixu) = q1u(ixu) * wvmap1(ixd)
                    end do
                 end do

              end if
           end if
        end do
*
* return
*
        return
*
* end of MAP1_F2D
*
        end
        
*
* --------------- UMAP2.F ---------------
*
 
        subroutine umap2
     &(
     &    mnx,mny,nx,ny,mnxu,nxu,x ,hdx ,y ,hdy ,p2u,
     &    mmx,mmy,mx,my,mmxu,mxu,xx,hdxx,yy,hdyy,q2u,
     &    nui,intp,imapcx,
     &    igx,igy,igrx,igry,
     &    w,lw,iw,liw,ier,conserved_var,smallx
     &)
*
*
*       SUBROUTINE UMAP2
*       ----------------
*
*       interpolates
*       p(i,j), (x(i),y(j)), i=1,...,nx, j=1,...,ny
*       onto
*       q(ii,jj), (xx(ii),yy(jj)), ii=1,...,mx and jj=1,...,my
*
*       see UMAP3 for more information
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

	INTEGER
     &  ier, igrx, igry, igx, igy, imapcx, intp, nui,
     &  mxu, mmy, mmx, nxu, mny, mnx, liw, lw,
     &  nx, ny, mx, my, mnxu, mmxu, iui, iiu,
     &  jy, j2, j3, j4, j5, j6, j7, j8, j9, j10,
     &  i2, i3, i4, i5,
     &  ii, ii2, ii3, ii4, ii5,
     &  jj, jj2, jj3, jj4, jj5

        REAL
     &  hdxx, hdyy, hdy, hdx, dxl, dxr, smallx

        INTEGER
     &          iw
        REAL
     &          w,
     &          x,y,
     &          xx,yy,
     &          p2u,
     &          q2u
        DIMENSION
     &          iw(liw),w(lw),
     &          x(nx),y(ny),
     &          xx(mx),yy(my),
     &          p2u(mnxu,ny),
     &          q2u(mmxu,my)
*
*
*
        ier = 0
*       ---------------------------------
*       set pointer in integer work space
        jy = mx + 1
*       ------------------
*       set index pointers
        j2  = 1
        j3  = j2 + my
        j4  = j3 + my
        j5  = j4 + my
        j6  = j5 + my
        j7  = j6 + mxu
        j8  = j7 + mxu
        j9  = j8 + mxu
        j10 = j9 + mxu
        if ( imapcx.eq.1 )
     &        call mapcx
     &  ( ny,y,hdy,my,yy,hdyy,
     &    iw(jy),w(j2),w(j3),w(j4),w(j5),
     &    intp,
     &    igy,igry,ier )
        i2 =  j10 + mxu
*       ----------------------------------
*       set work space portion and indices
*       which depend on x interpolation
        i3 = i2 + mx
        i4 = i3 + mx
        i5 = i4 + mx
        if ( imapcx.eq.1 )
     &        call mapcx
     &  ( nx,x,hdx,mx,xx,hdxx,
     &    iw    ,w(i2),w(i3),w(i4),w(i5),
     &    intp,
     &    igx,igrx,ier )
*       ---------------------
*       interpolation pointers
        do iui = 1,nui
           do ii = 1,mx
              iiu = (iui-1)*mx + ii
              ip_i_iiu (iiu) = (iui-1)*nx + iw(ii)
              ip_ii_iiu(iiu) = ii
           end do
        end do
*       ---------------
*       partial volumes
        if ( intp.eq.1 ) then
           do ii = 1,mx
              ii2 = i2 + ii - 1
              ii3 = i3 + ii - 1
              ii4 = i4 + ii - 1
              ii5 = i5 + ii - 1
              dxl         = w(ii2)
              dxr         = w(ii3)
              dvfl_ii(ii) = w(ii4)
              dvfr_ii(ii) = w(ii5)
              hdl_ii(ii)  = half*dxl
              hdr_ii(ii)  = half*dxr
           end do
           do jj = 1,my
              jj2 = j2 + jj - 1
              jj3 = j3 + jj - 1
              jj4 = j4 + jj - 1
              jj5 = j5 + jj - 1
              dxl         = w(jj2)
              dxr         = w(jj3)
              dvfl_jj(jj) = w(jj4)
              dvfr_jj(jj) = w(jj5)
              hdl_jj(jj)  = half*dxl
              hdr_jj(jj)  = half*dxr
           end do
        else if ( intp.eq.2  ) then
           do ii = 1,mx
              ii2 = i2 + ii - 1
              ii3 = i3 + ii - 1
              ii4 = i4 + ii - 1
              ii5 = i5 + ii - 1
              dxl         = w(ii2)
              dxr         = w(ii3)
              dvfl_ii(ii) = w(ii4)
              dvfr_ii(ii) = w(ii5)
              hdl_ii(ii)  = half*dxl
              hdr_ii(ii)  = half*dxr
              hd6l_ii(ii) = one - c2by3*dxl
              hd6r_ii(ii) = one - c2by3*dxr
           end do
           do jj = 1,my
              jj2 = j2 + jj - 1
              jj3 = j3 + jj - 1
              jj4 = j4 + jj - 1
              jj5 = j5 + jj - 1
              dxl         = w(jj2)
              dxr         = w(jj3)
              dvfl_jj(jj) = w(jj4)
              dvfr_jj(jj) = w(jj5)
              hdl_jj(jj)  = half*dxl
              hdr_jj(jj)  = half*dxr
              hd6l_jj(jj) = one - c2by3*dxl
              hd6r_jj(jj) = one - c2by3*dxr
           end do
        end if
        call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p2u,
     &                    mmx,mmy,mx,my,mmxu,mxu,q2u,
     &                    iw(jy),w(j6),w(j7),w(j8),w(j9),w(j10),
     &                          iw(1),
     &                    nui,intp,
     &                    igx,igy,igrx,igry,1,ier,
     &                    conserved_var,smallx )
*
* return
*
        return
*
* end of UMAP2
*
        end
 


        subroutine map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p2u,
     &                    mmx,mmy,mx,my,mmxu,mxu,q2u,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,isw,ier,
     &                    conserved_var,smallx )
*
*
*       SUBROUTINE MAP2
*       ---------------
*
*       interpolate in Y direction
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var, l_new_j

        INTEGER
     &  jy, ix, iw, liw, lw, nx, ny, mx, my, nxu, mny, mnx,
     &  ier, isw, igrx, igry, igx, igy, intp,
     &  nui, mxu, mnxu, mmxu, jj, j, iiu, jsave,
     &  mmy, mmx

        REAL
     &  hdx, hdy, dupl, dupr, dubs, dup, pl, pr,
     &  dv_lr_jj, hv_lr_jj, c1_jm1, c2_jm1, c4_jm1,
     &  c1_j, c2_j, c4_j, c1_jp1, c2_jp1, c4_jp1,
     &  c1_jp2, c2_jp2, c4_jp2, c5_jm1, c5_j, c5_jp1,
     &  c3_jm1, c3_j, c3_jp1,
     &  plmap_j, prmap_j, dpmap_j, p6map_j, dp_p6, dp_dp,
     &  hd_dv_l, hd_dv_r, v6_lr_jj, smallx

        REAL
     &          pjmm,pjm,pj,pjp,pjpp,
     &          p2u,
     &          q2u
        DIMENSION
     &          jy(my),ix(mx),
     &          pjmm(mxu),pjm(mxu),pj(mxu),pjp(mxu),pjpp(mxu),
     &          p2u(mnxu,ny),
     &          q2u(mmxu,my)
*
*
*
        if ( intp.eq.0 .or. igry.eq.1 ) then
*          -----------------------
*          no interpolation needed
*          -----------------------
           do jj = 1,my
              j = jy(jj)
              call map1 ( nxu,p2u(1,j),mxu,pj    ,
     &                    mx,intp,igx,igrx,0,ier ,
     &                    conserved_var,smallx )
              do iiu = 1,mxu
                 q2u(iiu,jj) = pj(iiu)
              end do
           end do
        else if ( intp.eq.1 ) then
*          --------------------
*          linear interpolation
*          --------------------
           jsave = -1
           do jj = 1,my
*             -------------------------------------------------------
*             load closest three j lines containing interpolate on xx
*             mesh for j-1,j,j+1 in pjm,pj,pjp
              j = jy(jj)
              if ( j.eq.jsave ) then
*                ---------------------------------------
*                j pointer has not moved since last pass
*                (no updates or interpolation)
                 l_new_j = .false.
              else if ( j.eq.jsave+1 ) then
*                --------------------------------
*                update j-1,j and interpolate j+1
                 l_new_j = .true.
                 do iiu = 1,mxu
                    pjm(iiu) = pj  (iiu)
                    pj (iiu) = pjp (iiu)
                 end do
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              else if ( j.eq.jsave+2 ) then
*                --------------------------------
*                update j-1 and interpolate j,j+1
                 l_new_j = .true.
                 do iiu = 1,mxu
                    pjm(iiu) = pjp (iiu)
                 end do
                 call map1 ( nxu,p2u(1,j  ),mxu,pj  ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              else
*                -------------------------------
*                interpolate all three j-1,j,j+1
                 l_new_j = .true.
                 call map1 ( nxu,p2u(1,j-1),mxu,pjm ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j  ),mxu,pj  ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              end if
*             ----------------------------
*             save j pointer for next pass
              jsave = j
*             ---------------------------------------------
*             linearly interpolate q(ii,jj) from pjm,pj,pjp
*             in y direction
              if ( l_new_j ) then
                 do iiu = 1,mxu
                    dupl = pj (iiu)-pjm(iiu)
                    dupr = pjp(iiu)-pj (iiu)
                    dubs = 1.e-3*(abs(dupl)+abs(dupr)+1.e-2)**2
                    dup  = half*( dupl+dupr )
     &                    *( (two*dupl*dupr + dubs)
     &                      /(dupl*dupl + dupr*dupr + dubs) )
                    dup  = sign( min(abs(dup),min(abs(dupl),abs(dupr))),
     &                     dup )
                    if ( dupl*dupr.le.zero ) dup = zero
                    pl = pj(iiu) - half*dup
                    pr = pj(iiu) + half*dup
                    pl = max( min(pjm(iiu),pj(iiu)),
     &                   min( max(pjm(iiu),pj(iiu)), pl) )
                    pr = max( min(pjp(iiu),pj(iiu)),
     &                   min( max(pjp(iiu),pj(iiu)), pr) )
                    plmap_iiu(iiu) = pl
                    dpmap_iiu(iiu) = pr - pl
                 end do
              end if
*             --------------
*             volume average
              dv_lr_jj = dvfr_jj(jj) - dvfl_jj(jj)
              hv_lr_jj = hdr_jj(jj)*dvfr_jj(jj) - hdl_jj(jj)*dvfl_jj(jj)
              do iiu = 1,mxu
                 q2u(iiu,jj) =
     &                         plmap_iiu(iiu)*dv_lr_jj
     &                        +dpmap_iiu(iiu)*hv_lr_jj
              end do                ! ii (iiu)
              if ( n_dens.gt.0 )
     &        call map2_f2d ( q2u,mmxu,mx,my,conserved_var,smallx, jj )
           end do                ! jj
        else if ( intp.eq.2 ) then
*          -----------------------
*          quadratic interpolation
*          -----------------------
           jsave = -4
           do jj = 1,my

*             ------------------------------------------------------
*             load closest five j lines containing interpolate on xx
*             mesh for j-2,j-1,j,j+1,j+2 in pjmm,pjm,pj,pjp,pjpp
              j = jy(jj)
              if ( j.eq.jsave ) then
*                ---------------------------------------
*                j pointer has not moved since last pass
*                (no updates or interpolation)
                 l_new_j = .false.
              else if ( j.eq.jsave+1 ) then
*                ----------------------------------------
*                update j-2,j-1,j,j-1 and interpolate j+2
                 l_new_j = .true.
                 do iiu = 1,mxu
                    pjmm (iiu) = pjm  (iiu)
                    pjm  (iiu) = pj   (iiu)
                    pj   (iiu) = pjp  (iiu)
                    pjp  (iiu) = pjpp (iiu)
                 end do
                 call map1 ( nxu,p2u(1,j+2),mxu,pjpp,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              else if ( j.eq.jsave+2 ) then
*                ----------------------------------------
*                update j-2,j-1,j and interpolate j+1,j+2
                 l_new_j = .true.
                 do iiu = 1,mxu
                    pjmm (iiu) = pj   (iiu)
                    pjm  (iiu) = pjp  (iiu)
                    pj   (iiu) = pjpp (iiu)
                 end do
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+2),mxu,pjpp,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              else if ( j.eq.jsave+3 ) then
*                ----------------------------------------
*                update j-2,j-1 and interpolate j,j+1,j+2
                 l_new_j = .true.
                 do iiu = 1,mxu
                    pjmm (iiu) = pjp  (iiu)
                    pjm  (iiu) = pjpp (iiu)
                 end do
                 call map1 ( nxu,p2u(1,j  ),mxu,pj  ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+2),mxu,pjpp,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              else if ( j.eq.jsave+4 ) then
*                ----------------------------------------
*                update j-2 and interpolate j-1,j,j+1,j+2
                 l_new_j = .true.
                 do iiu = 1,mxu
                    pjmm (iiu) = pjpp (iiu)
                 end do
                 call map1 ( nxu,p2u(1,j-1),mxu,pjm ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j  ),mxu,pj  ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+2),mxu,pjpp,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              else
*                --------------------------------------
*                interpolate all five j-2,j-1,j,j+1,j+2
                 l_new_j = .true.
                 call map1 ( nxu,p2u(1,j-2),mxu,pjmm,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j-1),mxu,pjm ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j  ),mxu,pj  ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+1),mxu,pjp ,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
                 call map1 ( nxu,p2u(1,j+2),mxu,pjpp,
     &                       mx,intp,igx,igrx,0,ier,
     &                       conserved_var,smallx )
              end if
*             ----------------------------
*             save j pointer for next pass
              jsave = j
*             ------------------------------------------
*             quadratically interpolate q2u(iiu,jj) from
*             pjmm,pjm,pj,pjp,pjpp in y direction
              if ( l_new_j ) then
                 do iiu = 1,mxu
*                   ------------------------
*                   zone to zone differences
                    c1_jm1 = pjm(iiu) - pjmm(iiu)
                    c2_jm1 = abs(c1_jm1 + c1_jm1)
                    c4_jm1 = sign(one , c1_jm1)
                    c1_j   = pj(iiu) - pjm(iiu)
                    c2_j   = abs(c1_j   + c1_j  )
                    c4_j   = sign(one , c1_j)
                    c1_jp1 = pjp(iiu) - pj(iiu)
                    c2_jp1 = abs(c1_jp1 + c1_jp1)
                    c4_jp1 = sign(one , c1_jp1)
                    c1_jp2 = pjpp(iiu) - pjp(iiu)
                    c2_jp2 = abs(c1_jp2 + c1_jp2)
                    c4_jp2 = sign(one , c1_jp2)
*                   --------------------
*                   parabola mean slopes
                    c5_jm1 = half*(c1_jm1 + c1_j  )
                    c5_j   = half*(c1_j   + c1_jp1)
                    c5_jp1 = half*(c1_jp1 + c1_jp2)
*                   ---------------------------
*                   slope limiting : CW84 (1.8)
*                   left zone
                    if ( c5_jm1.lt.zero ) then
                       c3_jm1 = -one
                    else
                       c3_jm1 = +one
                    end if
                    c5_jm1 = min( abs(c5_jm1),c2_jm1,c2_j ) * c3_jm1
                    if ( -c4_jm1*c4_j.ge.zero ) c5_jm1 = zero
*                   central zone
                    if ( c5_j.lt.zero ) then
                       c3_j = -one
                    else
                       c3_j = +one
                    end if
                    c5_j = min( abs(c5_j),c2_j,c2_jp1 ) * c3_j
                    if ( -c4_j*c4_jp1.ge.zero ) c5_j = zero
*                   right zone
                    if ( c5_jp1.lt.zero ) then
                       c3_jp1 = -one
                    else
                       c3_jp1 = +one
                    end if
                    c5_jp1 = min( abs(c5_jp1),c2_jp1,c2_jp2 ) * c3_jp1
                    if ( -c4_jp1*c4_jp2.ge.zero ) c5_jp1 = zero
*                   -----------------------------
*                   interface values : CW84 (1.6)
                    plmap_j =   pjm(iiu)
     &                        + half  * c1_j
     &                        - c1by6 * c5_j
     &                        + c1by6 * c5_jm1
                    prmap_j =   pj(iiu)
     &                        + half  * c1_jp1
     &                        - c1by6 * c5_jp1
     &                        + c1by6 * c5_j
*                   ----------------------------
*                   monotonization : CW84 (1.10)
                    if ( (prmap_j-pj(iiu))
     &                  *(pj(iiu)-plmap_j).le.zero ) then
                       prmap_j = pj(iiu)
                       plmap_j = pj(iiu)
                    end if
*                   ------------------------------------
*                   mean slope and quadratic coefficient
                    dpmap_j = prmap_j - plmap_j
                    p6map_j = six*( pj(iiu) - half*(plmap_j+prmap_j) )
*                   --------------------------------
*                   monotonize quadratic coefficient
                    dp_p6 =  dpmap_j*p6map_j
                    dp_dp =  dpmap_j*dpmap_j
                    if ( +dp_p6 - dp_dp.gt.zero ) then
                       plmap_j = three*pj(iiu) - two*prmap_j
                    else
     &                          if ( -dp_dp - dp_p6.gt.zero ) then
                       prmap_j = three*pj(iiu) - two*plmap_j
                    end if
*                   -------------
*                   final profile
                    plmap_iiu(iiu) = plmap_j
                    dpmap_iiu(iiu) = prmap_j - plmap_j
                    p6map_iiu(iiu) = six*( pj(iiu)
     &                                    -half*(prmap_j+plmap_j) )
                 end do                ! ii (iiu)
              end if                ! l_new_j
*             --------------
*             volume average
              hd_dv_l  = hdl_jj(jj)*dvfl_jj(jj)
              hd_dv_r  = hdr_jj(jj)*dvfr_jj(jj)
              dv_lr_jj = dvfr_jj(jj) - dvfl_jj(jj)
              hv_lr_jj = hd_dv_r - hd_dv_l
              v6_lr_jj = hd6r_jj(jj)*hd_dv_r - hd6l_jj(jj)*hd_dv_l
              do iiu = 1,mxu
                 q2u(iiu,jj) =
     &                         plmap_iiu(iiu)*dv_lr_jj
     &                        +dpmap_iiu(iiu)*hv_lr_jj
     &                        +p6map_iiu(iiu)*v6_lr_jj
              end do                ! ii (iiu)
              if ( n_dens.gt.0 )
     &        call map2_f2d ( q2u,mmxu,mx,my,conserved_var,smallx, jj )
           end do                ! jj
        else
*          ------------
*          unknown intp
*          ------------
           ier = 21
        end if
*
* return
*
        return
*
* end of MAP2
*
        end
 


        subroutine map2_f2d ( q2u,mmxu,mx,my,conserved_var,smallx, jj )
*
*
*       SUBROUTINE MAP2_F2D
*       -------------------
*
*       calculate total density from partial densities or renormalize
*       mass fractions (if .not.conserved_var)
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

        INTEGER
     &  jj, mx, mmxu, my, i_dens, ixdi, ixdf, ixd, ju, iu,
     &  ixui, ixuf, ixu

        REAL
     &  smallx, smallxdp

        REAL
     &          q2u
        DIMENSION
     &          q2u(mmxu,my)
*
*
*
        smallxdp = smallx

        do i_dens = 1,n_dens
           if ( nud(i_dens).gt.0 ) then
              if ( conserved_var ) then

                 ixdi = (iud(0,i_dens)-1)*mx + 1
                 ixdf = ixdi + mx - 1
                 do ixd = ixdi,ixdf
                    wvmap1(ixd) = zero
                 end do
                 do ju = 1,nud(i_dens)
                    iu = iud(ju,i_dens)
                    ixui = (iu-1)*mx + 1
                    ixuf = ixui + mx - 1
                    do ixu = ixui,ixuf
                       ixd = ixu - ixui + ixdi
                       wvmap1(ixd) = wvmap1(ixd) + q2u(ixu,jj)
                    end do
                 end do
                 do ixd = ixdi,ixdf
                    q2u(ixd,jj) = wvmap1(ixd)
                 end do

              else

                 ixdi = (iud(0,i_dens)-1)*mx + 1
                 ixdf = ixdi + mx - 1
                 do ixd = ixdi,ixdf
                    wvmap1(ixd) = zero
                 end do
                 do ju = 1,nud(i_dens)
                    iu = iud(ju,i_dens)
                    ixui = (iu-1)*mx + 1
                    ixuf = ixui + mx - 1
                    do ixu = ixui,ixuf
                       ixd = ixu - ixui + ixdi
                       wvmap1(ixd) = wvmap1(ixd)
     &                    +max(smallxdp,min(one, q2u(ixu,jj) ))
                    end do
                 end do
                 do ixd = ixdi,ixdf
                    wvmap1(ixd) = 1.e0 / wvmap1(ixd)
                 end do
                 do ju = 1,nud(i_dens)
                    iu = iud(ju,i_dens)
                    ixui = (iu-1)*mx + 1
                    ixuf = ixui + mx - 1
                    do ixu = ixui,ixuf
                       ixd = ixu - ixui + ixdi
                       q2u(ixu,jj) = q2u(ixu,jj) * wvmap1(ixd)
                    end do
                 end do

              end if
           end if
        end do
*
* return
*
        return
*
* end of MAP2_F2D
*
        end
        
*
* --------------- UMAP3.F ---------------
*
 
        subroutine umap3
     &(
     &    mnx,mny,mnz,nx,ny,nz,mnxu,nxu,x ,hdx ,y ,hdy ,z ,hdz ,p3u,
     &    mmx,mmy,mmz,mx,my,mz,mmxu,mxu,xx,hdxx,yy,hdyy,zz,hdzz,q3u,
     &    nui,intp,imapcx,
     &    igx,igy,igz,igrx,igry,igrz,
     &    w,lw,iw,liw,ier,conserved_var,smallx
     &)
*
*
*       SUBROUTINE UMAP3
*       ----------------
*
*       interpolates
*       p(i,j,k), (x(i),y(j),z(k))
*       i=1,...,nx; j=1,...,ny; k=1,...,nz
*
*       onto
*       q(ii,jj,kk), (xx(ii),yy(jj),zz(kk))
*       ii=1,...,mx; jj=1,...,my; kk=1,...,mz.
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPMACH.INC
*
*
*       machine-dependent floating-point constants
*
        REAL
     &
     &  bmach,emach,hemach,smach,hsmach,
     &        beexp,seexp,bexp,sexp,
     &  log_10_i
        COMMON /fpmach/
     &
     &  bmach,emach,hemach,smach,hsmach,
     &        beexp,seexp,bexp,sexp,
     &  log_10_i
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

	INTEGER
     &  liw, lw, nx, ny, nz, mx, my, mz, mnxu, mny, mmxu, mmy,
     &  k3, k4, k5, k6, k7, k8, k9, k10,
     &  j2, j3, j4, j5, j6, j7, j8, j9, j10,
     &  i2, i3, i4, i5,
     &  iui, iiu, iiujj, ii2, ii3, ii4, ii5,
     &  jj2, jj3, jj4, jj5,
     &  kk2, kk3, kk4, kk5,
     &  ier, igrz, igry, igrx, igz, igy, igx, imapcx, intp, nui,
     &  mxu, mmz, mmx, nxu, mnz, mnx, mxumy, nli,
     &  nlw1d, nlw2d, nlw3d, nlw, nliw,
     &  i, j, k, ii, jj, kk, jy, kz, k2

        REAL
     &  dxl, dxr, hdzz, hdyy, hdxx, hdx, hdy, hdz, edx, edy, edz,
     &  smallx

        INTEGER
     &          iw
        REAL
     &          w,
     &          x,y,z,
     &          xx,yy,zz,
     &          p3u,
     &          q3u
        DIMENSION
     &          iw(liw),w(lw),
     &          x(nx),y(ny),z(nz),
     &          xx(mx),yy(my),zz(mz),
     &          p3u(mnxu,mny,nz),
     &          q3u(mmxu,mmy,mz)
*
*       check input arguments
*
        ier = 1
        if ( min(mx,my,mz).lt.1 ) return
        ier = 2
        if ( min(nx,ny,nz).lt.1 ) return
        ier = 6
        if ( intp.lt.0 .or. intp.gt.2 ) return
        ier = 5
        mxumy = mxu*my
        nliw  = mx + my + mz
        nlw1d = 4*mx
        nlw2d = 4*my + 5*mx*nui
        nlw3d = 4*mz + 5*mx*my*nui
        nlw   = nlw1d + nlw2d + nlw3d
        if ( liw .lt. nliw ) return
        if ( lw  .lt. nlw  ) return
        ier = 3
        edx = hemach*(x(2)-x(1))
        edy = hemach*(y(2)-y(1))
        edz = hemach*(z(2)-z(1))
        if ( xx(1).lt.x(1)-edx .or. xx(mx).gt.x(nx)+edx ) return
        if ( yy(1).lt.y(1)-edy .or. yy(my).gt.y(ny)+edy ) return
        if ( zz(1).lt.z(1)-edz .or. zz(mz).gt.z(nz)+edz ) return
        ier = 411
        do i = 2,nx
           if ( x(i-1).ge.x(i) ) return
        end do
        ier = 412
        do j = 2,ny
           if ( y(j-1).ge.y(j) ) return
        end do
        ier = 413
        do k = 2,nz
           if ( z(k-1).ge.z(k) ) return
        end do
        ier = 421
        do ii = 2,mx
           if ( xx(ii-1).gt.xx(ii) ) return
        end do
        ier = 422
        do jj = 2,my
           if ( yy(jj-1).gt.yy(jj) ) return
        end do
        ier = 423
        do kk = 2,mz
           if ( zz(kk-1).gt.zz(kk) ) return
        end do
*       --------------
*       arguments o.k.
        ier = 0
        jy = mmx + 1
        kz = mmx + mmy + 1
*       -------------------------------------------
*       set indices which depend on z interpolation
        k2  = 1
        k3  = k2 + mmz
        k4  = k3 + mmz
        k5  = k4 + mmz
        k6  = k5 + mmz
        k7  = k6 + mxumy
        k8  = k7 + mxumy
        k9  = k8 + mxumy
        k10 = k9 + mxumy
        if ( imapcx.eq.1 )
     &        call mapcx
     &  ( nz,z,hdz,mz,zz,hdzz,
     &    iw(kz),w(k2),w(k3),w(k4),w(k5),
     &    intp,
     &    igz,igrz,ier )
        j2  = k10 + mxumy
*       -------------------------------------------
*       set indices which depend on y interpolation
        j3  = j2 + mmy
        j4  = j3 + mmy
        j5  = j4 + mmy
        j6  = j5 + mmy
        j7  = j6 + mxu
        j8  = j7 + mxu
        j9  = j8 + mxu
        j10 = j9 + mxu
        if ( imapcx.eq.1 )
     &        call mapcx
     &  ( ny,y,hdy,my,yy,hdyy,
     &    iw(jy),w(j2),w(j3),w(j4),w(j5),
     &    intp,
     &    igy,igry,ier )
        i2  = j10 + mxu
*       ----------------------------------
*       set work space portion and indices
*       which depend on x interpolation
        i3 = i2 + mmx
        i4 = i3 + mmx
        i5 = i4 + mmx
        if ( imapcx.eq.1 )
     &        call mapcx
     &  ( nx,x,hdx,mx,xx,hdxx,
     &    iw    ,w(i2),w(i3),w(i4),w(i5),
     &    intp,
     &    igx,igrx,ier )
*       ----------------------
*       interpolation pointers
        do iui = 1,nui
           do ii = 1,mx
              iiu = (iui-1)*mx + ii
              ip_i_iiu (iiu) = (iui-1)*nx + iw(ii)
              ip_ii_iiu(iiu) = ii
           end do
        end do
        do iiujj = 1,mxumy
           jj = int( ( iiujj+mxu  -1 ) / mxu   )
           ip_jj_iiujj(iiujj)  = jj
           ip_iiu_iiujj(iiujj) = iiujj - (jj-1) * mxu
        end do
*       ---------------
*       partial volumes
        if ( intp.eq.1 ) then
           do ii = 1,mx
              ii2 = i2 + ii - 1
              ii3 = i3 + ii - 1
              ii4 = i4 + ii - 1
              ii5 = i5 + ii - 1
              dxl         = w(ii2)
              dxr         = w(ii3)
              dvfl_ii(ii) = w(ii4)
              dvfr_ii(ii) = w(ii5)
              hdl_ii(ii)  = half*dxl
              hdr_ii(ii)  = half*dxr
           end do
           do jj = 1,my
              jj2 = j2 + jj - 1
              jj3 = j3 + jj - 1
              jj4 = j4 + jj - 1
              jj5 = j5 + jj - 1
              dxl         = w(jj2)
              dxr         = w(jj3)
              dvfl_jj(jj) = w(jj4)
              dvfr_jj(jj) = w(jj5)
              hdl_jj(jj)  = half*dxl
              hdr_jj(jj)  = half*dxr
           end do
           do kk = 1,mz
              kk2 = k2 + kk - 1
              kk3 = k3 + kk - 1
              kk4 = k4 + kk - 1
              kk5 = k5 + kk - 1
              dxl         = w(kk2)
              dxr         = w(kk3)
              dvfl_kk(kk) = w(kk4)
              dvfr_kk(kk) = w(kk5)
              
              hdl_kk(kk)  = half*dxl
              hdr_kk(kk)  = half*dxr
           end do
        else if ( intp.eq.2  ) then
           do ii = 1,mx
              ii2 = i2 + ii - 1
              ii3 = i3 + ii - 1
              ii4 = i4 + ii - 1
              ii5 = i5 + ii - 1
              dxl         = w(ii2)
              dxr         = w(ii3)
              dvfl_ii(ii) = w(ii4)
              dvfr_ii(ii) = w(ii5)
              hdl_ii(ii)  = half*dxl
              hdr_ii(ii)  = half*dxr
              hd6l_ii(ii) = one - c2by3*dxl
              hd6r_ii(ii) = one - c2by3*dxr
           end do
           do jj = 1,my
              jj2 = j2 + jj - 1
              jj3 = j3 + jj - 1
              jj4 = j4 + jj - 1
              jj5 = j5 + jj - 1
              dxl         = w(jj2)
              dxr         = w(jj3)
              dvfl_jj(jj) = w(jj4)
              dvfr_jj(jj) = w(jj5)
              hdl_jj(jj)  = half*dxl
              hdr_jj(jj)  = half*dxr
              hd6l_jj(jj) = one - c2by3*dxl
              hd6r_jj(jj) = one - c2by3*dxr
           end do
           do kk = 1,mz
              kk2 = k2 + kk - 1
              kk3 = k3 + kk - 1
              kk4 = k4 + kk - 1
              kk5 = k5 + kk - 1
              dxl         = w(kk2)
              dxr         = w(kk3)
              dvfl_kk(kk) = w(kk4)
              dvfr_kk(kk) = w(kk5)
              hdl_kk(kk)  = half*dxl
              hdr_kk(kk)  = half*dxr
              hd6l_kk(kk) = one - c2by3*dxl
              hd6r_kk(kk) = one - c2by3*dxr
           end do
        end if
        call map3
     &         ( mnx,mny,mnz,nx,ny,      nz,mnxu,nxu,p3u,
     &           mmx,mmy,mmz,mx,my,mxumy,mz,mmxu,mxu,q3u,
     &           iw(kz),w(k6),w(k7),w(k8),w(k9),w(k10),
     &           iw(jy),w(j6),w(j7),w(j8),w(j9),w(j10),
     &           iw    ,
     &           nui,intp,
     &           igx,igy,igz,igrx,igry,igrz,1,ier,
     &           conserved_var,smallx )
*
* return
*
        return
*
* end of UMAP3
*
        end
 

        subroutine map3
     &(
     &           mnx,mny,mnz,nx,ny,      nz,mnxu,nxu,p3u,
     &           mmx,mmy,mmz,mx,my,mxumy,mz,mmxu,mxu,q3u,
     &           kz,pkmm,pkm,pk,pkp,pkpp,
     &           jy,pjmm,pjm,pj,pjp,pjpp,
     &           ix,
     &           nui,intp,
     &           igx,igy,igz,igrx,igry,igrz,isw,ier,conserved_var,smallx
     &)
*
*
*       SUBROUTINE MAP3
*       ---------------
*
*       interpolate in Z direction
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var, l_new_k

	INTEGER
     &  ier, isw, igrz, igry, igrx, igz, igx, intp,
     &  kk, k, iiujj, jj, iiu, ksave,
     &  igy, nui, mmz, mmx, nxu, ny, nx, mnz, mnx

        REAL
     &  dupl, dupr, dubs, dup, pl, pr, dv_lr_kk, hv_lr_kk,
     &  c1_km1, c2_km1, c4_km1, c1_k, c2_k, c4_k,
     &  c1_kp1, c2_kp1, c4_kp1, c1_kp2, c2_kp2, c4_kp2,
     &  c5_km1, c5_k, c5_kp1, c3_km1, c3_k, c3_kp1,
     &  plmap_k, prmap_k, dpmap_k, p6map_k,
     &  dp_p6, dp_dp, hd_dv_l, hd_dv_r, v6_lr_kk, smallx

	INTEGER
     &          mxumy, mxu, mz, my, mx, kz, jy, ix,
     &          mnxu, mny, nz, mmxu, mmy

        REAL
     &          pkmm, pkm, pk, pkp, pkpp,
     &          pjmm, pjm, pj, pjp, pjpp,
     &          p3u, q3u

        DIMENSION
     &          pkmm(mxumy),pkm(mxumy),pk(mxumy),pkp(mxumy),pkpp(mxumy),
     &          pjmm(mxu),pjm(mxu),pj(mxu),pjp(mxu),pjpp(mxu),
     &          kz(mz),jy(my),ix(mx),
     &          p3u(mnxu,mny,nz),
     &          q3u(mmxu,mmy,mz)
*
*
*
        if ( intp.eq.0 .or. igrz.eq.1 ) then
           do kk = 1,mz
              k = kz(kk)
              call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k),
     &                    mmx,my,mx,my,mxu,mxu,pk ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              do iiujj = 1,mxumy
                 jj  = ip_jj_iiujj(iiujj)
                 iiu = ip_iiu_iiujj(iiujj)
                 q3u(iiu,jj,kk) = pk(iiujj)
              end do
           end do
        else if ( intp.eq.1 ) then
           ksave = -2
           do kk = 1,mz
              k = kz(kk)
              if ( k.eq.ksave ) then
*                ---------------------------------------
*                k pointer has not moved since last pass
*                (no updates or interpolation)
                 l_new_k = .false.
              else if ( k.eq.ksave+1 ) then
*                --------------------------------
*                update k-1,k and interpolate k+1
                 l_new_k = .true.
                 do iiujj = 1,mxumy
                    pkm(iiujj) = pk (iiujj)
                    pk (iiujj) = pkp(iiujj)
                 end do
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              else if ( k.eq.ksave+2 ) then
*                --------------------------------
*                update k-1 and interpolate k,k+1
                 l_new_k = .true.
                 do iiujj = 1,mxumy
                    pkm(iiujj) = pkp(iiujj)
                 end do
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k  ),
     &                    mmx,my,mx,my,mxu,mxu,pk  ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              else
*                ---------------------------------------------
*                interpolate all three k-1,k,k+1 in pkm,pk,pkp
                 l_new_k = .true.
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k-1),
     &                    mmx,my,mx,my,mxu,mxu,pkm ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k  ),
     &                    mmx,my,mx,my,mxu,mxu,pk  ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              end if
*             ----------------------------
*             save k pointer for next pass
              ksave = k
*             ----------------------------------------------------
*             linearly interpolate q(iiujj,k) from pkm,pk,pkp in z
              if ( l_new_k ) then
                 do iiujj = 1,mxumy
                    dupl = pk (iiujj)-pkm(iiujj)
                    dupr = pkp(iiujj)-pk (iiujj)
                    dubs = 1.e-3*(abs(dupl)+abs(dupr)+1.e-2)**2
                    dup  = half*( dupl+dupr )
     &                    *( (two*dupl*dupr + dubs)
     &                      /(dupl*dupl + dupr*dupr + dubs) )
                    dup  = sign( min(abs(dup),min(abs(dupl),abs(dupr))),
     &                     dup )
                    if ( dupl*dupr.le.zero ) dup = zero
                    pl = pk(iiujj) - half*dup
                    pr = pk(iiujj) + half*dup
                    pl = max( min(pkm(iiujj),pk(iiujj)),
     &                   min( max(pkm(iiujj),pk(iiujj)), pl) )
                    pr = max( min(pkp(iiujj),pk(iiujj)),
     &                   min( max(pkp(iiujj),pk(iiujj)), pr) )
                    plmap_iuj(iiujj) = pl
                    dpmap_iuj(iiujj) = pr - pl
                 end do
              end if
*             --------------
*             volume average
              dv_lr_kk = dvfr_kk(kk) - dvfl_kk(kk)
              hv_lr_kk = hdr_kk(kk)*dvfr_kk(kk) - hdl_kk(kk)*dvfl_kk(kk)
              do iiujj = 1,mxumy
                 jj  = ip_jj_iiujj(iiujj)
                 iiu = ip_iiu_iiujj(iiujj)
                 q3u(iiu,jj,kk) =
     &                            plmap_iuj(iiujj)*dv_lr_kk
     &                           +dpmap_iuj(iiujj)*hv_lr_kk
              end do                ! iiujj
              if ( n_dens.gt.0 )
     &              call map3_f2d ( q3u,mmxu,mmy,mz,mxumy,mxu,mx,my,
     &                        conserved_var,smallx, kk )
           end do                ! kk
        else if ( intp.eq.2 ) then
           ksave = -4
           do kk = 1,mz
*             ------------------------------------------------------
*             load closest five k lines containing interpolate on xy
*             mesh for k-2,k-1,k,k+1,k+2 in pkmm,pkm,pk,pkp,pkpp
              k = kz(kk)
              if ( k.eq.ksave ) then
*                ---------------------------------------
*                k pointer has not moved since last pass
*                (no updates or interpolation)
                 l_new_k = .false.
              else if ( k.eq.ksave+1 ) then
*                ----------------------------------------
*                update k-2,k-1,k,k+1 and interpolate k+2
                 l_new_k = .true.
                 do iiujj = 1,mxumy
                    pkmm (iiujj) = pkm  (iiujj)
                    pkm  (iiujj) = pk   (iiujj)
                    pk   (iiujj) = pkp  (iiujj)
                    pkp  (iiujj) = pkpp (iiujj)
                 end do
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+2),
     &                    mmx,my,mx,my,mxu,mxu,pkpp,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              else if ( k.eq.ksave+2 ) then
*                ----------------------------------------
*                update k-2,k-1,k and interpolate k+1,k+2
                 l_new_k = .true.
                 do iiujj = 1,mxumy
                    pkmm (iiujj) = pk   (iiujj)
                    pkm  (iiujj) = pkp  (iiujj)
                    pk   (iiujj) = pkpp (iiujj)
                 end do
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+2),
     &                    mmx,my,mx,my,mxu,mxu,pkpp,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              else if ( k.eq.ksave+3 ) then
*                ----------------------------------------
*                update k-2,k-1 and interpolate k,k+1,k+2
                 l_new_k = .true.
                 do iiujj = 1,mxumy
                    pkmm (iiujj) = pkp  (iiujj)
                    pkm  (iiujj) = pkpp (iiujj)
                 end do
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k  ),
     &                    mmx,my,mx,my,mxu,mxu,pk  ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+2),
     &                    mmx,my,mx,my,mxu,mxu,pkpp,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              else if ( k.eq.ksave+4 ) then
*                ----------------------------------------
*                update k-2 and interpolate k-1,k,k+1,k+2
                 l_new_k = .true.
                 do iiujj = 1,mxumy
                    pkmm (iiujj) = pkpp (iiujj)
                 end do
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k-1),
     &                    mmx,my,mx,my,mxu,mxu,pkm ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k  ),
     &                    mmx,my,mx,my,mxu,mxu,pk  ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+2),
     &                    mmx,my,mx,my,mxu,mxu,pkpp,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              else
*                --------------------------------------
*                interpolate all five k-2,k-1,k,k+1,k+2
                 l_new_k = .true.
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k-2),
     &                    mmx,my,mx,my,mxu,mxu,pkmm,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k-1),
     &                    mmx,my,mx,my,mxu,mxu,pkm ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k  ),
     &                    mmx,my,mx,my,mxu,mxu,pk  ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+1),
     &                    mmx,my,mx,my,mxu,mxu,pkp ,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
                 call map2
     &                  ( mnx,mny,nx,ny,mnxu,nxu,p3u(1,1,k+2),
     &                    mmx,my,mx,my,mxu,mxu,pkpp,
     &                    jy,pjmm,pjm,pj,pjp,pjpp,
     &                          ix,
     &                    nui,intp,
     &                    igx,igy,igrx,igry,0,ier,
     &                    conserved_var,smallx )
              end if
*             ----------------------------
*             save k pointer for next pass
              ksave = k
*             ------------------------------------------
*             quadratically interpolate q(iiujj,kk) from
*             pkmm,pkm,pk,pkp,pkpp in z direction
              if ( l_new_k ) then
                 do iiujj = 1,mxumy
*                   ------------------------
*                   zone to zone differences
                    c1_km1 = pkm(iiujj) - pkmm(iiujj)
                    c2_km1 = abs(c1_km1 + c1_km1)
                    c4_km1 = sign(one , c1_km1)
                    c1_k   = pk(iiujj) - pkm(iiujj)
                    c2_k   = abs(c1_k   + c1_k  )
                    c4_k   = sign(one , c1_k)
                    c1_kp1 = pkp(iiujj) - pk(iiujj)
                    c2_kp1 = abs(c1_kp1 + c1_kp1)
                    c4_kp1 = sign(one , c1_kp1)
                    c1_kp2 = pkpp(iiujj) - pkp(iiujj)
                    c2_kp2 = abs(c1_kp2 + c1_kp2)
                    c4_kp2 = sign(one , c1_kp2)
*                   --------------------
*                   parabola mean slopes
                    c5_km1 = half*(c1_km1 + c1_k  )
                    c5_k   = half*(c1_k   + c1_kp1)
                    c5_kp1 = half*(c1_kp1 + c1_kp2)
*                   ---------------------------
*                   slope limiting : CW84 (1.8)
*                   left zone
                    if ( c5_km1.lt.zero ) then
                       c3_km1 = -one
                    else
                       c3_km1 = +one
                    end if
                    c5_km1 = min( abs(c5_km1),c2_km1,c2_k ) * c3_km1
                    if ( -c4_km1*c4_k.ge.zero ) c5_km1 = zero
*                   central zone
                    if ( c5_k.lt.zero ) then
                       c3_k = -one
                    else
                       c3_k = +one
                    end if
                    c5_k = min( abs(c5_k),c2_k,c2_kp1 ) * c3_k
                    if ( -c4_k*c4_kp1.ge.zero ) c5_k = zero
*                   right zone
                    if ( c5_kp1.lt.zero ) then
                       c3_kp1 = -one
                    else
                       c3_kp1 = +one
                    end if
                    c5_kp1 = min( abs(c5_kp1),c2_kp1,c2_kp2 ) * c3_kp1
                    if ( -c4_kp1*c4_kp2.ge.zero ) c5_kp1 = zero
*                   -----------------------------
*                   interface values : CW84 (1.6)
                    plmap_k = pkm(iiujj)
     &                       + half  * c1_k
     &                       - c1by6 * c5_k
     &                       + c1by6 * c5_km1
                    prmap_k = pk(iiujj)
     &                       + half  * c1_kp1
     &                       - c1by6 * c5_kp1
     &                       + c1by6 * c5_k
*                   ----------------------------
*                   monotonization : CW84 (1.10)
                    if ( (prmap_k-pk(iiujj))
     &               *(pk(iiujj)-plmap_k).le.zero) then
                       prmap_k = pk(iiujj)
                       plmap_k = pk(iiujj)
                    end if
*                   ------------------------------------
*                   mean slope and quadratic coefficient
                    dpmap_k = prmap_k - plmap_k
                    p6map_k = six*( pk(iiujj) - half*(plmap_k+prmap_k) )
*                   --------------------------------
*                   monotonize quadratic coefficient
                    dp_p6 =  dpmap_k*p6map_k
                    dp_dp =  dpmap_k*dpmap_k
                    if ( +dp_p6 - dp_dp.gt.zero ) then
                       plmap_k = three*pk(iiujj) - two*prmap_k
                    else
     &                    if ( -dp_dp - dp_p6.gt.zero ) then
                       prmap_k = three*pk(iiujj) - two*plmap_k
                    end if
*                   -------------
*                   final profile
                    plmap_iuj(iiujj) = plmap_k
                    dpmap_iuj(iiujj) = prmap_k - plmap_k
                    p6map_iuj(iiujj) = six*( pk(iiujj)
     &                                      -half*(prmap_k+plmap_k) )
                 end do
              end if                ! l_new_k
*             --------------
*             volume average
              hd_dv_l  = hdl_kk(kk)*dvfl_kk(kk)
              hd_dv_r  = hdr_kk(kk)*dvfr_kk(kk)
              dv_lr_kk = dvfr_kk(kk) - dvfl_kk(kk)
              hv_lr_kk = hd_dv_r - hd_dv_l
              v6_lr_kk = hd6r_kk(kk)*hd_dv_r - hd6l_kk(kk)*hd_dv_l
              do iiujj = 1,mxumy
                 jj  = ip_jj_iiujj(iiujj)
                 iiu = ip_iiu_iiujj(iiujj)
                 q3u(iiu,jj,kk) =
     &                            plmap_iuj(iiujj)*dv_lr_kk
     &                           +dpmap_iuj(iiujj)*hv_lr_kk
     &                           +p6map_iuj(iiujj)*v6_lr_kk
              end do                ! iiujj
              if ( n_dens.gt.0 )
     &              call map3_f2d ( q3u,mmxu,mmy,mz,mxumy,mxu,mx,my,
     &                        conserved_var,smallx, kk )
           end do                ! kk
        else
*          ------------
*          unknown intp
*          ------------
           ier = 31
        end if
*
* return
*
        return
*
* end of MAP3
*
        end
 


        subroutine map3_f2d ( q3u,mmxu,mmy,mz,mxumy,mxu,mx,my,
     &                        conserved_var,smallx, kk )
*
*
*       SUBROUTINE MAP3_F2D
*       -------------------
*
*       calculate total density from partial densities or renormalize
*       mass fractions (if .not.conserved_var)
*
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*        PRECIS.INC
*
        IMPLICIT NONE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       FPCONST.INC
*
*
*       replacements for floating point literal constants
*
        REAL
     &
     &  zero, half, one, two, three, four, six, ten,
     &  c1by3, c1by6, c2by3, c4by3, c5by3, c7by12, cm1by12,
     &  qrtr, c1by5, c1by7, c3by4, c1by12

        PARAMETER (  zero    =  0.e0       )
        PARAMETER (  qrtr    =  0.25d0     )
        PARAMETER (  half    =  0.5d0      )
        PARAMETER (  one     =  1.e0       )
        PARAMETER (  two     =  2.e0       )
        PARAMETER (  three   =  3.e0       )
        PARAMETER (  four    =  4.e0       )
        PARAMETER (  six     =  6.e0       )
        PARAMETER (  ten     =  10.e0      )
        PARAMETER (  c1by3   =  1.e0/3.e0  )
        PARAMETER (  c1by5   =  1.e0/5.e0  )
        PARAMETER (  c1by6   =  1.e0/6.e0  )
        PARAMETER (  c1by7   =  1.e0/7.e0  )
        PARAMETER (  c2by3   =  2.e0/3.e0  )
        PARAMETER (  c4by3   =  4.e0/3.e0  )
        PARAMETER (  c5by3   =  5.e0/3.e0  )
        PARAMETER (  c3by4   =  3.e0/4.e0  )
        PARAMETER (  c7by12  =  7.e0/12.e0 )
        PARAMETER (  c1by12  =  1.e0/12.e0 )
        PARAMETER (  cm1by12 = -1.e0/12.e0 )
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*
*
*       AMRMAPS.INC
*
*
*       AMR working storage for MAPPING module (umap subroutines)
*
*
        INTEGER
     &
     &  mvx_m, mvm_m, mvxu_m, mvxmu_m, mui_m

#include "umap.h"

        INTEGER
     &
     &  n_dens, nud, iud,
     &  ip_i_iiu,
     &  ip_ii_iiu,
     &  ip_jj_iiujj,
     &  ip_iiu_iiujj

        COMMON  /amrmapi_p/
     &
     &  n_dens, nud(mui_m), iud(0:mui_m,mui_m),
     &  ip_i_iiu     (mvxu_m),
     &  ip_ii_iiu    (mvxu_m),
     &  ip_jj_iiujj  (mvxmu_m),
     &  ip_iiu_iiujj (mvxmu_m)

        REAL
     &
     &  hdl_ii         , hdr_ii         ,
     &  hd6l_ii        , hd6r_ii        ,
     &  dvfl_ii        , dvfr_ii        ,
     &  plmapc         , prmapc         , 
     &  dpmapc         , p6mapc         ,
     &  wvmap1         , wvmap2         ,
     &  wvmap3         , wvmap4

        COMMON  /amrmapr1d_p/
     &
     &  hdl_ii  (mvx_m), hdr_ii  (mvx_m),
     &  hd6l_ii (mvx_m), hd6r_ii (mvx_m),
     &  dvfl_ii (mvx_m), dvfr_ii (mvx_m),
     &  plmapc (mvxu_m), prmapc (mvxu_m), 
     &  dpmapc (mvxu_m), p6mapc (mvxu_m),
     &  wvmap1 (mvxu_m), wvmap2 (mvxu_m),
     &  wvmap3 (mvxu_m), wvmap4 (mvxu_m)

        REAL
     &
     &  hdl_jj            ,  hdr_jj         ,
     &  hd6l_jj           ,  hd6r_jj        ,
     &  dvfl_jj           ,  dvfr_jj        ,
     &  plmap_iiu         ,
     &  dpmap_iiu         ,
     &  p6map_iiu

        COMMON  /amrmapr2d_p/
     &
     &  hdl_jj     (mvx_m),  hdr_jj  (mvx_m),
     &  hd6l_jj    (mvx_m),  hd6r_jj (mvx_m),
     &  dvfl_jj    (mvx_m),  dvfr_jj (mvx_m),
     &  plmap_iiu (mvxu_m),
     &  dpmap_iiu (mvxu_m),
     &  p6map_iiu (mvxu_m)

        REAL
     &
     &  hdl_kk             , hdr_kk         ,
     &  hd6l_kk            , hd6r_kk        ,
     &  dvfl_kk            , dvfr_kk        ,
     &  plmap_iuj          ,
     &  dpmap_iuj          ,
     &  p6map_iuj

        COMMON  /amrmapr3d_p/
     &
     &  hdl_kk  (mvx_m)    , hdr_kk  (mvx_m),
     &  hd6l_kk (mvx_m)    , hd6r_kk (mvx_m),
     &  dvfl_kk (mvx_m)    , dvfr_kk (mvx_m),
     &  plmap_iuj (mvxmu_m),
     &  dpmap_iuj (mvxmu_m),
     &  p6map_iuj (mvxmu_m)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
        LOGICAL
     &  conserved_var

	INTEGER
     &  kk, my, mx, mxy, mxumy, mxu,
     &  i_dens, jj, ixdi, ixdf, ixd, ju, iu, ixui, ixuf, ixu

        REAL
     &  smallx, smallxdp

	INTEGER
     &  mmxu, mmy, mz

        REAL
     &  q3u

        DIMENSION
     &          q3u(mmxu,mmy,mz)
*
*
*
        smallxdp = smallx

        do i_dens = 1,n_dens
           if (  nud(i_dens).gt.0 ) then
              if ( conserved_var ) then

                 do jj = 1,my
                    ixdi = (iud(0,i_dens)-1)*mx + 1
                    ixdf = ixdi + mx - 1
                    do ixd = ixdi,ixdf
                       wvmap1(ixd) = zero
                    end do
                    do ju = 1,nud(i_dens)
                       iu = iud(ju,i_dens)
                       ixui = (iu-1)*mx + 1
                       ixuf = ixui + mx - 1
                       do ixu = ixui,ixuf
                          ixd = ixu - ixui + ixdi
                          wvmap1(ixd) = wvmap1(ixd) + q3u(ixu,jj,kk)
                       end do
                    end do
                    do ixd = ixdi,ixdf
                       q3u(ixd,jj,kk) = wvmap1(ixd)
                    end do
                 end do

              else

                 do jj = 1,my
                    ixdi = (iud(0,i_dens)-1)*mx + 1
                    ixdf = ixdi + mx - 1
                    do ixd = ixdi,ixdf
                       wvmap1(ixd) = zero
                    end do
                    do ju = 1,nud(i_dens)
                       iu = iud(ju,i_dens)
                       ixui = (iu-1)*mx + 1
                       ixuf = ixui + mx - 1
                       do ixu = ixui,ixuf
                          ixd = ixu - ixui + ixdi
                          wvmap1(ixd) = wvmap1(ixd)
     &                                 +max(smallxdp,min(one,
     &                                                 q3u(ixu,jj,kk) ))
                       end do
                    end do
                    do ixd = ixdi,ixdf
                       wvmap1(ixd) = 1.e0 / wvmap1(ixd)
                    end do
                    do ju = 1,nud(i_dens)
                       iu = iud(ju,i_dens)
                       ixui = (iu-1)*mx + 1
                       ixuf = ixui + mx - 1
                       do ixu = ixui,ixuf
                          ixd = ixu - ixui + ixdi
                          q3u(ixu,jj,kk) = q3u(ixu,jj,kk) * wvmap1(ixd)
                       end do
                    end do
                 end do

              end if
           end if
        end do
*
* return
*
        return
*
* end of MAP3_F2D
*
        end
