#include "Flash.h"
      subroutine calc_shielding(yn, dl, temp, abh2, abco, npixels,
     $              cell_cdto_map, cell_cdH2_map, cell_cdCO_map,
     $                          fshield_H2, fshield_CO, AV_mean, 
     $                          chi_mean)

      implicit NONE
      INTEGER npixels
      REAL yn, dl, abh2, abco, temp
      REAL cell_cdto_map(npixels)
      REAL cell_cdH2_map(npixels), cell_cdCO_map(npixels)
      REAL fshield_H2, fshield_CO, AV_mean, chi_mean
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5

#include "param.h"
#include "cool.h"
#include "constants.h"
#include "cons.h"
      REAL AV, chi
      REAL sigma_dust
      REAL local_NH, AV_local, local_NH2, local_NCO, f_dust
      REAL max_cdto_map
      external chi
      integer I, nray, nray_start, nray_end
c
c Total column density of HI nuclei within characteristic scale-length dl
c


      local_NH = dl * yn
      AV_local = AV_conversion_factor*dust_to_gas_ratio*local_NH

      sigma_dust=2e-21*dust_to_gas_ratio
c      NH_ext= AV_ext/AV_conversion_factor
c
c H2 self-shielding
c
#if CHEMISTRYNETWORK == 4 
      local_NH2 = local_NH * abh2

      call calc_H2_shielding_factor(temp, local_NH2, fshield_H2,
     $     npixels, cell_cdH2_map)

#endif
c
c CO self-shielding, shielding of CO by H2
c
#if CHEMISTRYNETWORK == 5
      local_NH2 = local_NH * abh2

!      print*, "local NH2=", local_NH2
      call calc_H2_shielding_factor(temp, local_NH2, fshield_H2,
     $     npixels, cell_cdH2_map)


      local_NCO = local_NH * abco
      call calc_CO_shielding_factor(local_NCO, local_NH2, fshield_CO,
     $     npixels, cell_cdH2_map, cell_cdCO_map)
#endif
c
c Dust shielding -- we compute a mean AV, assuming a weighting factor exp(-2.5AV); this
c is appropriate for photoelectric heating and CO photodissociation, shouldn't be too
c much in error for other things
c
      nray     = npixels
      f_dust   = 0e0
      chi_mean = 0e0
      AV       = 0e0
      max_cdto_map = 0e0

c      AV=AV_local + AV_ext
c      f_dust=exp(-2.5e0*AV)

      do I = 1, nray
c        AV = AV_conversion_factor *cell_cdto_map(I)!!+ AV_local + AV_ext
c        f_dust=f_dust+(1e0/real(nray))*exp(-sigma_dust* 
c     $        (cell_cdto_map(I)+NH_ext)) 

        AV = AV_local+ AV_ext + AV_conversion_factor
     $        *dust_to_gas_ratio*cell_cdto_map(I)
        f_dust   = f_dust + (1e0 / real(nray)) * exp(-2.5e0 * AV)
        chi_mean = chi_mean + (1e0 / real(nray)) * chi(AV)


        max_cdto_map = max(max_cdto_map,cell_cdto_map(I))

      enddo
      AV_mean  = -log(f_dust) / 2.5e0

!      print*,'[calc_shielding:] AV_mean,chi_mean=',AV_mean, chi_mean
!      print*,'AV_local, AV_ext, local_NH, dl, yn='
!      print*, AV_local,AV_ext,local_NH,dl,yn
!      print*, "abh2, abco =", abh2, abco
!      print*,'NH2, NCO=',local_NH2, local_NCO
!      print*,'max(cell_cdto_map)=',max_cdto_map
!      print*,'fshield_H2, fshield_CO=',fshield_H2, fshield_CO

#endif
      return
      end subroutine calc_shielding
c=======================================================================
c
c    \\\\\\\\\\          E N D     F U N C T I O N        //////////
c    //////////             C A L C _ P H O T O           \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////                C H I                      \\\\\\\\\\
c
c=======================================================================
c
      REAL function chi(AV)
      implicit NONE
      REAL AV
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
#include "param.h"
#include "cool.h"
#include "isrf.h"
      REAL dchi, dAV
      integer index

      if (AV .ne. AV) then
        print*, "Error: AV is NaN"
        stop
      endif

      if (AV .eq. 0) then
        chi = 1e0
        return
      elseif (AV .gt. 2e2) then
c Assume chi remains constant for AV > 200
        if (isrf_option .eq. 1) then
          chi = chi_tab_mathis(200)
        else
          chi = chi_tab_black(200)
        endif
        return
      else
        index = int(AV)
        if ( index.eq. 200 ) then
           print *, "Index = 200... Whoops."
           stop
        end if 
        if (isrf_option .eq. 1) then
          if (index .eq. 0) then
            dchi = chi_tab_mathis(1) - 1e0
            dAV  = AV - index
            chi  = 1e0 + dAV * dchi
          else
            dchi = chi_tab_mathis(index+1) - chi_tab_mathis(index)
            dAV  = AV - index
            chi  = chi_tab_mathis(index) + dAV * dchi
          endif
        else
          if (index .eq. 0) then
            dchi = chi_tab_black(1) - 1e0
            dAV  = AV - index
            chi  = 1e0 + dAV * dchi
          else
            dchi = chi_tab_black(index+1) - chi_tab_black(index)
            dAV  = AV - index
            chi  = chi_tab_black(index) + dAV * dchi
          endif
        endif
        return
      endif
#endif
      return
      end function chi
c
c=======================================================================
c
c
c    \\\\\\\\\\          E N D   F U N C T I O N          //////////
c    //////////                  C H I                    \\\\\\\\\\
c
c=======================================================================
c
      subroutine calc_H2_shielding_factor(temp, local_NH2, fshield_H2,
     $     npixels, cell_cdH2_map)


      use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
      implicit none
      INTEGER npixels
      REAL fshield_H2, temp, local_NH2
      REAL cell_cdH2_map(npixels)
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
#include "param.h"
#include "cool.h"
!      REAL column_density_projection_H2(npixels)
      REAL b, b5, h2dis0, h2dis1
      REAL kboltzmann, mh
!      REAL fshield_cell

      integer I

      call PhysicalConstants_get("Boltzmann", kboltzmann)
      call PhysicalConstants_get("proton mass", mh)

c      print*,'kb, mh=',kboltzmann,mh
c
c To compute H2 self-shielding, we assume that the local b param. is 
c representative of the whole flow. This isn't true in general, but 
c is still a minor approximation compared to the others we are making
c
      b  = sqrt(kboltzmann * temp / mh)
      b5 = b / 1e5
      fshield_H2 = 0e0
!      fshield_cell = 0.965/(1e0+local_NH2/5e14/b5)**2.e0
!     $     +3.5e-2*exp(-8.5e-4*sqrt(1.e0+local_NH2/5e14))/
!     $     sqrt(1e0+local_NH2/5e14)


      do I = 1, npixels
!!!        cell_cdH2_map(I)=0.e0
        h2dis0  = (local_NH2 + cell_cdH2_map(I)) / 5e14
        h2dis1  = sqrt(1e0 + h2dis0)
        fshield_H2 = fshield_H2 + (1e0 / real(npixels)) * (9.65e-1 /
     $             ((1e0 + h2dis0 / b5) * (1e0 + h2dis0 / b5))
     $             + 3.5e-2 * exp(-8.5e-4 * h2dis1) / h2dis1)


!        print*,I," fshield_H2=",fshield_H2,b5,local_NH2,cell_cdh2_map(I)
      enddo

!      if (fshield_cell .gt. fshield_H2) then
!         print*,"map contribution",fshield_cell, fshield_H2
!         print*, cell_cdH2_map
!      endif
!      if (fshield_cell .lt. fshield_H2) then 
!         print*,"Wrong shielding?",fshield_cell, fshield_H2
!         print*, cell_cdH2_map
!      endif

      return
#endif
      end subroutine calc_H2_shielding_factor
c
c
c
c =============================================================== 
      subroutine calc_CO_shielding_factor(local_NCO,local_NH2,
     $     fshield_CO, npixels, cell_cdH2_map, cell_cdCO_map)

      implicit none
      INTEGER npixels
      REAL cell_cdH2_map(npixels)
      REAL cell_cdCO_map(npixels)
      REAL fshield_CO, local_NCO, local_NH2
#if CHEMISTRYNETWORK == 5
#include "param.h"
#include "cool.h"
#include "shield_data.h"
!      REAL column_density_projection_h2(npixels)
!      REAL column_density_projection_co(npixels)
      REAL theta_CO, theta_H2, NCO, NH2
      REAL dth, dNact
      integer nray, I, ipos
c
c CO self-shielding, shielding by H2.
c
c Our treatment here is based on Lee et al (1996)
c
!      print *, "cs: ", local_NCO,local_NH2,
!     $     fshield_CO, npixels, cell_cdH2_map, cell_cdCO_map
      fshield_CO = 1e0

      nray = npixels
      theta_CO   = 0e0
      theta_H2   = 0e0
      do I = 1, nray
        NCO = cell_cdCO_map(I) + local_NCO
        if (NCO .le. NCO_shield_min) then
          theta_CO = theta_CO + 1e0 / real(nray)
        elseif (NCO .ge. NCO_shield_max) then
          theta_CO = theta_CO + (1e0 / real(nray))
     $             * CO_self_shielding(imax_COss)
        else
          ipos  = int(log10(NCO) / dNshield) - ioff_co
          if (ipos .eq. 0) then
            theta_CO = theta_CO + 1e0 / real(nray)
          else
            dth = CO_self_shielding(ipos+1) - 
     $            CO_self_shielding(ipos)
            dNact = log10(NCO) - (ipos + ioff_co) * dNshield
            theta_CO = theta_CO + (1e0/real(nray))
     $               * (CO_self_shielding(ipos) 
     $               + dth * (dNact / dNshield))
          endif
        endif
c
        NH2 = cell_cdH2_map(I) + local_NH2
        if (NH2 .le. NH2_shield_min) then
          theta_H2 = theta_H2 + 1e0 / real(nray)
        elseif (NH2 .ge. NH2_shield_max) then
          theta_H2 = theta_H2 + (1e0 / real(nray) )
     $             * CO_H2_shielding(imax_COH2)
        else
          ipos  = int(log10(NH2) / dNshield) - ioff_h2
          if (ipos .eq. 0) then
            theta_H2 = theta_H2 + 1e0 /real(nray)
          else
            dth   = CO_H2_shielding(ipos+1) - 
     $              CO_H2_shielding(ipos)
            dNact = log10(NH2) - (ipos + ioff_h2) * dNshield
            theta_H2 = theta_H2 + (1e0 / real(nray))
     $               * (CO_H2_shielding(ipos) 
     $               + dth * (dNact / dNshield))
          endif
        endif
       enddo

       fshield_CO = theta_H2 * theta_CO
#endif
       return
       end subroutine calc_CO_shielding_factor
