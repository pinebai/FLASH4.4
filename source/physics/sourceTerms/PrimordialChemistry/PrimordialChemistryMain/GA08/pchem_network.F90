!!***if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_network
!!
!! NAME
!!
!! pchem_network
!!
!! SYNOPSIS
!!
!! pchem_network(real, intent(IN)               :: tt,
!!              real, intent(IN), dimension(:)  :: y,
!!              real, intent(OUT), dimentsion(:) :: dydt
!!
!! DESCRIPTION
!!
!! this routine sets up the system of ode's for the globchem chemistry 
!! reaction, this is all the reactions given in Glover and Abel
!!
!! Ions: H,H+,H-,D,D+,D-,He++,He+,He,H2+,H2,HD+,HD,ELEC
!!
!! ARGUMENTS
!!
!! tt -- Not used
!! y  -- Molar mass fractions
!! dydt -- ODES
!!
!!***

subroutine pchem_network(tt,ys,dydt)

   use PrimordialChemistry_data
   use pchem_data

   implicit none

#include "constants.h"
#include "Flash.h"

   real, intent(IN)   :: tt,ys(*)
   real, intent(OUT)  :: dydt(*)
   integer i, ie

!! setup the system of odes:


do i=1,nisotp
   dydt(i) = 0.0e0
enddo

ie = iELEC

dydt(iHP) = -ratraw(iR003)*ys(iH)*ys(iHP)+ratraw(iR004)*ys(iH)*ys(iH2P) &
          & -ratraw(iR005)*ys(iHM)*ys(iHP)-ratraw(iR007)*ys(iH2)*ys(iHP) &
          & +ratraw(iR012)*ys(iH)*ys(ie)-ratraw(iR013)*ys(iHP)*ys(ie) &
          & -ratraw(iR016)*ys(iHP)*ys(iHM)+ratraw(iR024)*ys(iH2)*ys(iHEP) &
          & +ratraw(iR026)*ys(iHEP)*ys(iH)-ratraw(iR027)*ys(iHE)*ys(iHP) &
          & -ratraw(iR034)*ys(iD)*ys(iHP)+ratraw(iR035)*ys(iH)*ys(iDP) &
          & +ratraw(iR038)*ys(iHDP)*ys(iH)+ratraw(iR039)*ys(iH2)*ys(iDP) &
          & -ratraw(iR041)*ys(iHD)*ys(iHP)-ratraw(iR042)*ys(iD)*ys(iHP) &
          & -ratraw(iR060)*ys(iHP)*ys(iDM)-ratraw(iR067)*ys(iHP)*ys(iDM) &
          & +ratraw(iR082)*ys(iH2P)*ys(iD)+ratraw(iR085)*ys(iHDP)*ys(iD) &
          & +ratraw(iR087)*ys(iH)*ys(iD2P)-ratraw(iR092)*ys(iHD)*ys(iHP) &
          & -ratraw(iR093)*ys(iHD)*ys(iHP)+ratraw(iR095)*ys(iHD)*ys(iDP) &
          & -ratraw(iR097)*ys(iD2)*ys(iHP)-ratraw(iR098)*ys(iD2)*ys(iHP) &
          & -ratraw(iR099)*ys(iD2)*ys(iHP)+ratraw(iR102)*ys(iHD)*ys(iHEP) &
          & +ratraw(iR118)*ys(iH2P)+ratraw(iR120)*ys(iHDP)
 
dydt(iH) =  -ratraw(iR001)*ys(iH)*ys(ie)-ratraw(iR002)*ys(iHM)*ys(iH) &
          & -ratraw(iR003)*ys(iH)*ys(iHP)-ratraw(iR004)*ys(iH)*ys(iH2P) &
          & +ratraw(iR005)*ys(iHM)*ys(iHP)+ratraw(iR005)*ys(iHM)*ys(iHP) &
          & +ratraw(iR006)*ys(iH2P)*ys(ie)+ratraw(iR006)*ys(iH2P)*ys(ie) &
          & +ratraw(iR007)*ys(iH2)*ys(iHP)+ratraw(iR008)*ys(iH2)*ys(ie) &
          & +ratraw(iR008)*ys(iH2)*ys(ie)-ratraw(iR009)*ys(iH2)*ys(iH) &
          & +ratraw(iR009)*ys(iH2)*ys(iH)+ratraw(iR009)*ys(iH2)*ys(iH) &
          & +ratraw(iR009)*ys(iH2)*ys(iH)+ratraw(iR010)*ys(iH2)*ys(iH2) &
          & +ratraw(iR010)*ys(iH2)*ys(iH2)+ratraw(iR011)*ys(iH2)*ys(iHE) &
          & +ratraw(iR011)*ys(iH2)*ys(iHE)-ratraw(iR012)*ys(iH)*ys(ie) &
          & +ratraw(iR013)*ys(iHP)*ys(ie)+ratraw(iR014)*ys(iHM)*ys(ie) &
          & -ratraw(iR015)*ys(iHM)*ys(iH)+ratraw(iR015)*ys(iHM)*ys(iH) &
          & +ratraw(iR015)*ys(iHM)*ys(iH)+ratraw(iR021)*ys(iHM)*ys(iH2P) &
          & +ratraw(iR022)*ys(iHM)*ys(iH2P)+ratraw(iR022)*ys(iHM)*ys(iH2P) &
          & +ratraw(iR022)*ys(iHM)*ys(iH2P)+ratraw(iR023)*ys(iH2)*ys(ie) &
          & +ratraw(iR024)*ys(iH2)*ys(iHEP)-ratraw(iR026)*ys(iHEP)*ys(iH) &
          & +ratraw(iR027)*ys(iHE)*ys(iHP)+ratraw(iR028)*ys(iHEP)*ys(iHM) &
          & +ratraw(iR029)*ys(iHE)*ys(iHM)-ratraw(iR030)*ys(iH)*ys(iH)*ys(iH) &
          & -ratraw(iR030)*ys(iH)*ys(iH)*ys(iH)-ratraw(iR030)*ys(iH)*ys(iH)*ys(iH) &
          & +ratraw(iR030)*ys(iH)*ys(iH)*ys(iH)-ratraw(iR031)*ys(iH)*ys(iH)*ys(iH2) &
          & -ratraw(iR031)*ys(iH)*ys(iH)*ys(iH2)-ratraw(iR032)*ys(iH)*ys(iH)*ys(iHE) &
          & -ratraw(iR032)*ys(iH)*ys(iH)*ys(iHE)+ratraw(iR034)*ys(iD)*ys(iHP) &
          & -ratraw(iR035)*ys(iH)*ys(iDP)-ratraw(iR036)*ys(iH)*ys(iD) &
          & +ratraw(iR037)*ys(iH2)*ys(iD)-ratraw(iR038)*ys(iHDP)*ys(iH) &
          & -ratraw(iR040)*ys(iHD)*ys(iH)-ratraw(iR043)*ys(iH)*ys(iDP) &
          & +ratraw(iR044)*ys(iHDP)*ys(ie)+ratraw(iR048)*ys(iH2P)*ys(iD) &
          & -ratraw(iR050)*ys(iHDP)*ys(iH)-ratraw(iR052)*ys(iH)*ys(iDM) &
          & +ratraw(iR053)*ys(iD)*ys(iHM)-ratraw(iR055)*ys(iH)*ys(iDM) &
          & +ratraw(iR057)*ys(iHD)*ys(ie)-ratraw(iR064)*ys(iDM)*ys(iH) &
          & +ratraw(iR064)*ys(iDM)*ys(iH)+ratraw(iR066)*ys(iDP)*ys(iHM) &
          & +ratraw(iR067)*ys(iHP)*ys(iDM)+ratraw(iR070)*ys(iH2P)*ys(iDM) &
          & +ratraw(iR070)*ys(iH2P)*ys(iDM)+ratraw(iR071)*ys(iHDP)*ys(iHM) &
          & +ratraw(iR072)*ys(iHDP)*ys(iHM)+ratraw(iR072)*ys(iHDP)*ys(iHM) &
          & +ratraw(iR074)*ys(iHDP)*ys(iDM)+ratraw(iR075)*ys(iD2P)*ys(iHM) &
          & +ratraw(iR076)*ys(iD2P)*ys(iHM)-ratraw(iR083)*ys(iHDP)*ys(iH) &
          & +ratraw(iR084)*ys(iHDP)*ys(iD)-ratraw(iR087)*ys(iH)*ys(iD2P) &
          & -ratraw(iR088)*ys(iD2P)*ys(iH)-ratraw(iR089)*ys(iD2P)*ys(iH) &
          & +ratraw(iR091)*ys(iH2)*ys(iDP)+ratraw(iR092)*ys(iHD)*ys(iHP) &
          & +ratraw(iR096)*ys(iHD)*ys(iDP)+ratraw(iR099)*ys(iD2)*ys(iHP) &
          & +ratraw(iR103)*ys(iHD)*ys(iHEP)+ratraw(iR106)*ys(iHD)*ys(iD) &
          & -ratraw(iR107)*ys(iD2)*ys(iH)-ratraw(iR108)*ys(iHD)*ys(iH) &
          & +ratraw(iR108)*ys(iHD)*ys(iH)+ratraw(iR108)*ys(iHD)*ys(iH) &
          & +ratraw(iR109)*ys(iHD)*ys(iH2)+ratraw(iR110)*ys(iHD)*ys(iHE) &
          & +ratraw(iR111)*ys(iHD)*ys(ie)-ratraw(iR112)*ys(iD2)*ys(iH) &
          & +ratraw(iR112)*ys(iD2)*ys(iH)+ratraw(iR116)*ys(iHM) &
          & +ratraw(iR118)*ys(iH2P)+ratraw(iR119)*ys(iHDP) &
          & +ratraw(iR122)*ys(iH2)+ratraw(iR122)*ys(iH2)+ratraw(iR123)*ys(iHD)

dydt(iHM) =  ratraw(iR001)*ys(iH)*ys(ie)-ratraw(iR002)*ys(iHM)*ys(iH) &
          & -ratraw(iR005)*ys(iHM)*ys(iHP)-ratraw(iR014)*ys(iHM)*ys(ie) &
          & -ratraw(iR015)*ys(iHM)*ys(iH)-ratraw(iR016)*ys(iHP)*ys(iHM) &
          & -ratraw(iR021)*ys(iHM)*ys(iH2P)-ratraw(iR022)*ys(iHM)*ys(iH2P) &
          & +ratraw(iR023)*ys(iH2)*ys(ie)-ratraw(iR028)*ys(iHEP)*ys(iHM) &
          & -ratraw(iR029)*ys(iHE)*ys(iHM)+ratraw(iR052)*ys(iH)*ys(iDM) &
          & -ratraw(iR053)*ys(iD)*ys(iHM)-ratraw(iR054)*ys(iD)*ys(iHM) &
          & +ratraw(iR058)*ys(iHD)*ys(ie)-ratraw(iR061)*ys(iDP)*ys(iHM) &
          & -ratraw(iR066)*ys(iDP)*ys(iHM)-ratraw(iR071)*ys(iHDP)*ys(iHM) &
          & -ratraw(iR072)*ys(iHDP)*ys(iHM)-ratraw(iR075)*ys(iD2P)*ys(iHM) &
          & -ratraw(iR076)*ys(iD2P)*ys(iHM)-ratraw(iR116)*ys(iHM)

dydt(iH2P) =  ratraw(iR003)*ys(iH)*ys(iHP)-ratraw(iR004)*ys(iH)*ys(iH2P) &
           & -ratraw(iR006)*ys(iH2P)*ys(ie)+ratraw(iR007)*ys(iH2)*ys(iHP) &
           & +ratraw(iR016)*ys(iHP)*ys(iHM)-ratraw(iR021)*ys(iHM)*ys(iH2P) &
           & -ratraw(iR022)*ys(iHM)*ys(iH2P)+ratraw(iR025)*ys(iH2)*ys(iHEP) &
           & -ratraw(iR048)*ys(iH2P)*ys(iD)+ratraw(iR050)*ys(iHDP)*ys(iH) &
           & -ratraw(iR069)*ys(iH2P)*ys(iDM)-ratraw(iR070)*ys(iH2P)*ys(iDM) &
           & -ratraw(iR081)*ys(iD)*ys(iH2P)-ratraw(iR082)*ys(iH2P)*ys(iD) &
           & +ratraw(iR090)*ys(iH2)*ys(iDP)+ratraw(iR093)*ys(iHD)*ys(iHP) &
           & -ratraw(iR118)*ys(iH2P)

dydt(iH2) =  ratraw(iR002)*ys(iHM)*ys(iH)+ratraw(iR004)*ys(iH)*ys(iH2P) &
          & -ratraw(iR007)*ys(iH2)*ys(iHP)-ratraw(iR008)*ys(iH2)*ys(ie) &
          & -ratraw(iR009)*ys(iH2)*ys(iH)-ratraw(iR010)*ys(iH2)*ys(iH2) &
          & -ratraw(iR010)*ys(iH2)*ys(iH2)+ratraw(iR010)*ys(iH2)*ys(iH2) &
          & -ratraw(iR011)*ys(iH2)*ys(iHE)+ratraw(iR021)*ys(iHM)*ys(iH2P) &
          & -ratraw(iR023)*ys(iH2)*ys(ie)-ratraw(iR024)*ys(iH2)*ys(iHEP) &
          & -ratraw(iR025)*ys(iH2)*ys(iHEP)+ratraw(iR030)*ys(iH)*ys(iH)*ys(iH) &
          & -ratraw(iR031)*ys(iH)*ys(iH)*ys(iH2)+ratraw(iR031)*ys(iH)*ys(iH)*ys(iH2) &
          & +ratraw(iR031)*ys(iH)*ys(iH)*ys(iH2)+ratraw(iR032)*ys(iH)*ys(iH)*ys(iHE) &
          & -ratraw(iR037)*ys(iH2)*ys(iD)-ratraw(iR039)*ys(iH2)*ys(iDP) &
          & +ratraw(iR040)*ys(iHD)*ys(iH)+ratraw(iR041)*ys(iHD)*ys(iHP) &
          & +ratraw(iR069)*ys(iH2P)*ys(iDM)+ratraw(iR081)*ys(iD)*ys(iH2P) &
          & +ratraw(iR083)*ys(iHDP)*ys(iH)-ratraw(iR090)*ys(iH2)*ys(iDP) &
          & -ratraw(iR091)*ys(iH2)*ys(iDP)-ratraw(iR109)*ys(iHD)*ys(iH2) &
          & +ratraw(iR109)*ys(iHD)*ys(iH2)-ratraw(iR113)*ys(iD2)*ys(iH2) &
          & +ratraw(iR113)*ys(iD2)*ys(iH2)-ratraw(iR122)*ys(iH2)

dydt(iDP) =  -ratraw(iR033)*ys(iDP)*ys(ie)+ratraw(iR034)*ys(iD)*ys(iHP) &
           & -ratraw(iR035)*ys(iH)*ys(iDP)-ratraw(iR039)*ys(iH2)*ys(iDP) &
           & +ratraw(iR041)*ys(iHD)*ys(iHP)-ratraw(iR043)*ys(iH)*ys(iDP) &
           & +ratraw(iR045)*ys(iD)*ys(ie)+ratraw(iR046)*ys(iHEP)*ys(iD) &
           & -ratraw(iR047)*ys(iHE)*ys(iDP)+ratraw(iR049)*ys(iHDP)*ys(iD) &
           & -ratraw(iR061)*ys(iDP)*ys(iHM)-ratraw(iR062)*ys(iDP)*ys(iDM) &
           & -ratraw(iR066)*ys(iDP)*ys(iHM)-ratraw(iR068)*ys(iDP)*ys(iDM) &
           & -ratraw(iR080)*ys(iD)*ys(iDP)+ratraw(iR081)*ys(iD)*ys(iH2P) &
           & +ratraw(iR083)*ys(iHDP)*ys(iH)+ratraw(iR086)*ys(iD)*ys(iD2P) &
           & +ratraw(iR089)*ys(iD2P)*ys(iH)-ratraw(iR090)*ys(iH2)*ys(iDP) &
           & -ratraw(iR091)*ys(iH2)*ys(iDP)-ratraw(iR094)*ys(iHD)*ys(iDP) &
           & -ratraw(iR095)*ys(iHD)*ys(iDP)-ratraw(iR096)*ys(iHD)*ys(iDP) &
           & +ratraw(iR097)*ys(iD2)*ys(iHP)-ratraw(iR100)*ys(iD2)*ys(iDP) &
           & +ratraw(iR103)*ys(iHD)*ys(iHEP)+ratraw(iR105)*ys(iD2)*ys(iHEP) &
           & +ratraw(iR119)*ys(iHDP)+ratraw(iR121)*ys(iD2P)

dydt(iD) =  ratraw(iR033)*ys(iDP)*ys(ie)-ratraw(iR034)*ys(iD)*ys(iHP) &
         & +ratraw(iR035)*ys(iH)*ys(iDP)-ratraw(iR036)*ys(iH)*ys(iD) &
         & -ratraw(iR037)*ys(iH2)*ys(iD)+ratraw(iR040)*ys(iHD)*ys(iH) &
         & -ratraw(iR042)*ys(iD)*ys(iHP)+ratraw(iR044)*ys(iHDP)*ys(ie) &
         & -ratraw(iR045)*ys(iD)*ys(ie)-ratraw(iR046)*ys(iHEP)*ys(iD) &
         & +ratraw(iR047)*ys(iHE)*ys(iDP)-ratraw(iR048)*ys(iH2P)*ys(iD) &
         & -ratraw(iR049)*ys(iHDP)*ys(iD)+ratraw(iR050)*ys(iHDP)*ys(iH) &
         & -ratraw(iR051)*ys(iD)*ys(ie)+ratraw(iR052)*ys(iH)*ys(iDM) &
         & -ratraw(iR053)*ys(iD)*ys(iHM)-ratraw(iR054)*ys(iD)*ys(iHM) &
         & -ratraw(iR056)*ys(iD)*ys(iDM)+ratraw(iR058)*ys(iHD)*ys(ie) &
         & +ratraw(iR059)*ys(iD2)*ys(ie)+ratraw(iR063)*ys(iDM)*ys(ie) &
         & +ratraw(iR064)*ys(iDM)*ys(iH)+ratraw(iR065)*ys(iDM)*ys(iHE) &
         & +ratraw(iR066)*ys(iDP)*ys(iHM)+ratraw(iR067)*ys(iHP)*ys(iDM) &
         & +ratraw(iR068)*ys(iDP)*ys(iDM)+ratraw(iR068)*ys(iDP)*ys(iDM) &
         & +ratraw(iR069)*ys(iH2P)*ys(iDM)+ratraw(iR070)*ys(iH2P)*ys(iDM) &
         & +ratraw(iR072)*ys(iHDP)*ys(iHM)+ratraw(iR073)*ys(iHDP)*ys(iDM) &
         & +ratraw(iR074)*ys(iHDP)*ys(iDM)+ratraw(iR074)*ys(iHDP)*ys(iDM) &
         & +ratraw(iR076)*ys(iD2P)*ys(iHM)+ratraw(iR076)*ys(iD2P)*ys(iHM) &
         & +ratraw(iR077)*ys(iD2P)*ys(iDM)+ratraw(iR078)*ys(iD2P)*ys(iDM) &
         & +ratraw(iR078)*ys(iD2P)*ys(iDM)+ratraw(iR078)*ys(iD2P)*ys(iDM) &
         & +ratraw(iR079)*ys(iHEP)*ys(iDM)-ratraw(iR080)*ys(iD)*ys(iDP) &
         & -ratraw(iR081)*ys(iD)*ys(iH2P)-ratraw(iR082)*ys(iH2P)*ys(iD) &
         & -ratraw(iR084)*ys(iHDP)*ys(iD)-ratraw(iR085)*ys(iHDP)*ys(iD) &
         & -ratraw(iR086)*ys(iD)*ys(iD2P)+ratraw(iR088)*ys(iD2P)*ys(iH) &
         & +ratraw(iR090)*ys(iH2)*ys(iDP)+ratraw(iR093)*ys(iHD)*ys(iHP) &
         & +ratraw(iR094)*ys(iHD)*ys(iDP)+ratraw(iR098)*ys(iD2)*ys(iHP) &
         & +ratraw(iR100)*ys(iD2)*ys(iDP)+ratraw(iR102)*ys(iHD)*ys(iHEP) &
         & +ratraw(iR105)*ys(iD2)*ys(iHEP)-ratraw(iR106)*ys(iHD)*ys(iD) &
         & +ratraw(iR107)*ys(iD2)*ys(iH)+ratraw(iR108)*ys(iHD)*ys(iH) &
         & +ratraw(iR109)*ys(iHD)*ys(iH2)+ratraw(iR110)*ys(iHD)*ys(iHE) &
         & +ratraw(iR111)*ys(iHD)*ys(ie)+ratraw(iR112)*ys(iD2)*ys(iH) &
         & +ratraw(iR112)*ys(iD2)*ys(iH)+ratraw(iR113)*ys(iD2)*ys(iH2) &
         & +ratraw(iR113)*ys(iD2)*ys(iH2)+ratraw(iR114)*ys(iD2)*ys(iHE) &
         & +ratraw(iR114)*ys(iD2)*ys(iHE)+ratraw(iR115)*ys(iD2)*ys(ie) &
         & +ratraw(iR115)*ys(iD2)*ys(ie)+ratraw(iR117)*ys(iDM) &
         & +ratraw(iR120)*ys(iHDP)+ratraw(iR121)*ys(iD2P)+ratraw(iR123)*ys(iHD) &
         & +ratraw(iR124)*ys(iD2)+ratraw(iR124)*ys(iD2)

dydt(iDM) =  ratraw(iR051)*ys(iD)*ys(ie)-ratraw(iR052)*ys(iH)*ys(iDM) &
          & +ratraw(iR053)*ys(iD)*ys(iHM)-ratraw(iR055)*ys(iH)*ys(iDM) &
          & -ratraw(iR056)*ys(iD)*ys(iDM)+ratraw(iR057)*ys(iHD)*ys(ie) &
          & +ratraw(iR059)*ys(iD2)*ys(ie)-ratraw(iR060)*ys(iHP)*ys(iDM) &
          & -ratraw(iR062)*ys(iDP)*ys(iDM)-ratraw(iR063)*ys(iDM)*ys(ie) &
          & -ratraw(iR064)*ys(iDM)*ys(iH)-ratraw(iR065)*ys(iDM)*ys(iHE) &
          & -ratraw(iR067)*ys(iHP)*ys(iDM)-ratraw(iR068)*ys(iDP)*ys(iDM) &
          & -ratraw(iR069)*ys(iH2P)*ys(iDM)-ratraw(iR070)*ys(iH2P)*ys(iDM) &
          & -ratraw(iR073)*ys(iHDP)*ys(iDM)-ratraw(iR074)*ys(iHDP)*ys(iDM) &
          & -ratraw(iR077)*ys(iD2P)*ys(iDM)-ratraw(iR078)*ys(iD2P)*ys(iDM) &
          & -ratraw(iR079)*ys(iHEP)*ys(iDM)-ratraw(iR117)*ys(iDM)

dydt(iHDP) =  -ratraw(iR038)*ys(iHDP)*ys(iH)+ratraw(iR042)*ys(iD)*ys(iHP) &
            & +ratraw(iR043)*ys(iH)*ys(iDP)-ratraw(iR044)*ys(iHDP)*ys(ie) &
            & +ratraw(iR048)*ys(iH2P)*ys(iD)-ratraw(iR049)*ys(iHDP)*ys(iD) &
            & -ratraw(iR050)*ys(iHDP)*ys(iH)+ratraw(iR060)*ys(iHP)*ys(iDM) &
            & +ratraw(iR061)*ys(iDP)*ys(iHM)-ratraw(iR071)*ys(iHDP)*ys(iHM) &
            & -ratraw(iR072)*ys(iHDP)*ys(iHM)-ratraw(iR073)*ys(iHDP)*ys(iDM) &
            & -ratraw(iR074)*ys(iHDP)*ys(iDM)-ratraw(iR083)*ys(iHDP)*ys(iH) &
            & -ratraw(iR084)*ys(iHDP)*ys(iD)-ratraw(iR085)*ys(iHDP)*ys(iD) &
            & +ratraw(iR088)*ys(iD2P)*ys(iH)+ratraw(iR091)*ys(iH2)*ys(iDP) &
            & +ratraw(iR092)*ys(iHD)*ys(iHP)+ratraw(iR094)*ys(iHD)*ys(iDP) &
            & +ratraw(iR098)*ys(iD2)*ys(iHP)+ratraw(iR101)*ys(iHD)*ys(iHEP) &
            & -ratraw(iR119)*ys(iHDP)-ratraw(iR120)*ys(iHDP)

dydt(iHD) =  ratraw(iR036)*ys(iH)*ys(iD)+ratraw(iR037)*ys(iH2)*ys(iD) &
          & +ratraw(iR038)*ys(iHDP)*ys(iH)+ratraw(iR039)*ys(iH2)*ys(iDP) &
          & -ratraw(iR040)*ys(iHD)*ys(iH)-ratraw(iR041)*ys(iHD)*ys(iHP) &
          & +ratraw(iR049)*ys(iHDP)*ys(iD)+ratraw(iR054)*ys(iD)*ys(iHM) &
          & +ratraw(iR055)*ys(iH)*ys(iDM)-ratraw(iR057)*ys(iHD)*ys(ie) &
          & -ratraw(iR058)*ys(iHD)*ys(ie)+ratraw(iR071)*ys(iHDP)*ys(iHM) &
          & +ratraw(iR073)*ys(iHDP)*ys(iDM)+ratraw(iR082)*ys(iH2P)*ys(iD) &
          & +ratraw(iR089)*ys(iD2P)*ys(iH)-ratraw(iR092)*ys(iHD)*ys(iHP) &
          & -ratraw(iR093)*ys(iHD)*ys(iHP)-ratraw(iR094)*ys(iHD)*ys(iDP) &
          & -ratraw(iR095)*ys(iHD)*ys(iDP)-ratraw(iR096)*ys(iHD)*ys(iDP) &
          & +ratraw(iR097)*ys(iD2)*ys(iHP)-ratraw(iR101)*ys(iHD)*ys(iHEP) &
          & -ratraw(iR102)*ys(iHD)*ys(iHEP)-ratraw(iR103)*ys(iHD)*ys(iHEP) &
          & -ratraw(iR106)*ys(iHD)*ys(iD)+ratraw(iR107)*ys(iD2)*ys(iH) &
          & -ratraw(iR108)*ys(iHD)*ys(iH)-ratraw(iR109)*ys(iHD)*ys(iH2) &
          & -ratraw(iR110)*ys(iHD)*ys(iHE)-ratraw(iR111)*ys(iHD)*ys(ie) &
          & -ratraw(iR123)*ys(iHD)

dydt(iD2P) =  ratraw(iR062)*ys(iDP)*ys(iDM)-ratraw(iR075)*ys(iD2P)*ys(iHM) &
           & -ratraw(iR076)*ys(iD2P)*ys(iHM)-ratraw(iR077)*ys(iD2P)*ys(iDM) &
           & -ratraw(iR078)*ys(iD2P)*ys(iDM)+ratraw(iR080)*ys(iD)*ys(iDP) &
           & +ratraw(iR084)*ys(iHDP)*ys(iD)-ratraw(iR086)*ys(iD)*ys(iD2P) &
           & -ratraw(iR087)*ys(iH)*ys(iD2P)-ratraw(iR088)*ys(iD2P)*ys(iH) &
           & -ratraw(iR089)*ys(iD2P)*ys(iH)+ratraw(iR096)*ys(iHD)*ys(iDP) &
           & +ratraw(iR099)*ys(iD2)*ys(iHP)+ratraw(iR100)*ys(iD2)*ys(iDP) &
           & +ratraw(iR104)*ys(iD2)*ys(iHEP)-ratraw(iR121)*ys(iD2P)

dydt(iD2) =  ratraw(iR056)*ys(iD)*ys(iDM)-ratraw(iR059)*ys(iD2)*ys(ie) &
          & +ratraw(iR075)*ys(iD2P)*ys(iHM)+ratraw(iR077)*ys(iD2P)*ys(iDM) &
          & +ratraw(iR085)*ys(iHDP)*ys(iD)+ratraw(iR086)*ys(iD)*ys(iD2P) &
          & +ratraw(iR087)*ys(iH)*ys(iD2P)+ratraw(iR095)*ys(iHD)*ys(iDP) &
          & -ratraw(iR097)*ys(iD2)*ys(iHP)-ratraw(iR098)*ys(iD2)*ys(iHP) &
          & -ratraw(iR099)*ys(iD2)*ys(iHP)-ratraw(iR100)*ys(iD2)*ys(iDP) &
          & -ratraw(iR104)*ys(iD2)*ys(iHEP)-ratraw(iR105)*ys(iD2)*ys(iHEP) &
          & +ratraw(iR106)*ys(iHD)*ys(iD)-ratraw(iR107)*ys(iD2)*ys(iH) &
          & -ratraw(iR112)*ys(iD2)*ys(iH)-ratraw(iR113)*ys(iD2)*ys(iH2) &
          & -ratraw(iR114)*ys(iD2)*ys(iHE)-ratraw(iR115)*ys(iD2)*ys(ie) &
          & -ratraw(iR124)*ys(iD2)

dydt(iHEP) =   ratraw(iR017)*ys(iHE)*ys(ie)-ratraw(iR018)*ys(iHEP)*ys(ie) &
            & -ratraw(iR019)*ys(iHEP)*ys(ie)+ratraw(iR020)*ys(iHEPP)*ys(ie) &
            & -ratraw(iR024)*ys(iH2)*ys(iHEP)-ratraw(iR025)*ys(iH2)*ys(iHEP) &
            & -ratraw(iR026)*ys(iHEP)*ys(iH)+ratraw(iR027)*ys(iHE)*ys(iHP) &
            & -ratraw(iR028)*ys(iHEP)*ys(iHM)-ratraw(iR046)*ys(iHEP)*ys(iD) &
            & +ratraw(iR047)*ys(iHE)*ys(iDP)-ratraw(iR079)*ys(iHEP)*ys(iDM) &
            & -ratraw(iR101)*ys(iHD)*ys(iHEP)-ratraw(iR102)*ys(iHD)*ys(iHEP) &
            & -ratraw(iR103)*ys(iHD)*ys(iHEP)-ratraw(iR104)*ys(iD2)*ys(iHEP) &
            & -ratraw(iR105)*ys(iD2)*ys(iHEP)

dydt(iHE) =  -ratraw(iR011)*ys(iH2)*ys(iHE)+ratraw(iR011)*ys(iH2)*ys(iHE) &
           & -ratraw(iR017)*ys(iHE)*ys(ie)+ratraw(iR019)*ys(iHEP)*ys(ie) &
           & +ratraw(iR024)*ys(iH2)*ys(iHEP)+ratraw(iR025)*ys(iH2)*ys(iHEP) &
           & +ratraw(iR026)*ys(iHEP)*ys(iH)-ratraw(iR027)*ys(iHE)*ys(iHP) &
           & +ratraw(iR028)*ys(iHEP)*ys(iHM)-ratraw(iR029)*ys(iHE)*ys(iHM) &
           & +ratraw(iR029)*ys(iHE)*ys(iHM)-ratraw(iR032)*ys(iH)*ys(iH)*ys(iHE) &
           & +ratraw(iR032)*ys(iH)*ys(iH)*ys(iHE)+ratraw(iR046)*ys(iHEP)*ys(iD) &
           & -ratraw(iR047)*ys(iHE)*ys(iDP)-ratraw(iR065)*ys(iDM)*ys(iHE) &
           & +ratraw(iR065)*ys(iDM)*ys(iHE)+ratraw(iR079)*ys(iHEP)*ys(iDM) &
           & +ratraw(iR101)*ys(iHD)*ys(iHEP)+ratraw(iR102)*ys(iHD)*ys(iHEP) &
           & +ratraw(iR103)*ys(iHD)*ys(iHEP)+ratraw(iR104)*ys(iD2)*ys(iHEP) &
           & +ratraw(iR105)*ys(iD2)*ys(iHEP)-ratraw(iR110)*ys(iHD)*ys(iHE) &
           & +ratraw(iR110)*ys(iHD)*ys(iHE)-ratraw(iR114)*ys(iD2)*ys(iHE) &
           & +ratraw(iR114)*ys(iD2)*ys(iHE)

dydt(iHEPP) =  ratraw(iR018)*ys(iHEP)*ys(ie)-ratraw(iR020)*ys(iHEPP)*ys(ie)

dydt(iELEC) =  -ratraw(iR001)*ys(iH)*ys(ie)+ratraw(iR002)*ys(iHM)*ys(iH) &
             & -ratraw(iR006)*ys(iH2P)*ys(ie)-ratraw(iR008)*ys(iH2)*ys(ie) &
             & +ratraw(iR008)*ys(iH2)*ys(ie)-ratraw(iR012)*ys(iH)*ys(ie) &
             & +ratraw(iR012)*ys(iH)*ys(ie)+ratraw(iR012)*ys(iH)*ys(ie) &
             & -ratraw(iR013)*ys(iHP)*ys(ie)-ratraw(iR014)*ys(iHM)*ys(ie) &
             & +ratraw(iR014)*ys(iHM)*ys(ie)+ratraw(iR014)*ys(iHM)*ys(ie) &
             & +ratraw(iR015)*ys(iHM)*ys(iH)+ratraw(iR016)*ys(iHP)*ys(iHM) &
             & -ratraw(iR017)*ys(iHE)*ys(ie)+ratraw(iR017)*ys(iHE)*ys(ie) &
             & +ratraw(iR017)*ys(iHE)*ys(ie)-ratraw(iR018)*ys(iHEP)*ys(ie) &
             & +ratraw(iR018)*ys(iHEP)*ys(ie)+ratraw(iR018)*ys(iHEP)*ys(ie) &
             & -ratraw(iR019)*ys(iHEP)*ys(ie)-ratraw(iR020)*ys(iHEPP)*ys(ie) &
             & -ratraw(iR023)*ys(iH2)*ys(ie)+ratraw(iR029)*ys(iHE)*ys(iHM) &
             & -ratraw(iR033)*ys(iDP)*ys(ie)-ratraw(iR044)*ys(iHDP)*ys(ie) &
             & -ratraw(iR045)*ys(iD)*ys(ie)+ratraw(iR045)*ys(iD)*ys(ie) &
             & +ratraw(iR045)*ys(iD)*ys(ie)-ratraw(iR051)*ys(iD)*ys(ie) &
             & +ratraw(iR054)*ys(iD)*ys(iHM)+ratraw(iR055)*ys(iH)*ys(iDM) &
             & +ratraw(iR056)*ys(iD)*ys(iDM)-ratraw(iR057)*ys(iHD)*ys(ie) &
             & -ratraw(iR058)*ys(iHD)*ys(ie)-ratraw(iR059)*ys(iD2)*ys(ie) &
             & +ratraw(iR060)*ys(iHP)*ys(iDM)+ratraw(iR061)*ys(iDP)*ys(iHM) &
             & +ratraw(iR062)*ys(iDP)*ys(iDM)-ratraw(iR063)*ys(iDM)*ys(ie) &
             & +ratraw(iR063)*ys(iDM)*ys(ie)+ratraw(iR063)*ys(iDM)*ys(ie) &
             & +ratraw(iR064)*ys(iDM)*ys(iH)+ratraw(iR065)*ys(iDM)*ys(iHE) &
             & -ratraw(iR111)*ys(iHD)*ys(ie)+ratraw(iR111)*ys(iHD)*ys(ie) &
             & -ratraw(iR115)*ys(iD2)*ys(ie)+ratraw(iR115)*ys(iD2)*ys(ie) &
             & +ratraw(iR116)*ys(iHM)+ratraw(iR117)*ys(iDM)

   return
end subroutine pchem_network
                

