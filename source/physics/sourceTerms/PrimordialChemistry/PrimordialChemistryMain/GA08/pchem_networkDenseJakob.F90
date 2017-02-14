!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_networkDenseJakob
!!
!!
!! NAME
!!
!! pchem_networkDenseJakob
!!
!! SYNOPSIS
!!
!! call pchem_networkDenseJakok (real, intent(IN)    :: tt
!!                              real, intent(INOUT) :: y(:),
!!                              real, intent(OUT)   :: dfdy(nphys,nphys),
!!                           integer, intent(IN)    :: nlog,
!!                           integer, intent(IN)    :: nphys)
!!
!!
!! DESCRIPTION
!!
!!  routine networkDenseJakob sets up the dense PrimordialChemistry jacobian
!!
!!
!! ARGUMENTS
!!
!! tt: Not used
!! y:  Molar mass fractions
!! dfdy: Jacobian elements
!! nlog: logical size of dfdy
!! nphys: real size of dfdy
!!
!!***

subroutine pchem_networkDenseJakob(tt,ys,dfdy,nlog,nphys)

  use PrimordialChemistry_data
  use pchem_data
  use PrimordialChemistry_dataEOS

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)  :: nlog, nphys
  real, intent(IN)     :: tt
  real, intent(INOUT), dimension(*) :: ys
  real, intent(OUT), dimension(nphys,nphys) :: dfdy
!  real, intent(OUT), dimension(nlog,nlog) :: dfdy
  integer  i,j, ie
 
 

 !! zero the jacobian
 do j = 1,nlog
    do i = 1,nlog
       dfdy(i,j) = 0.0e0
    enddo
 enddo



!! positive definite mass fractions ??? What is this going to do?
 do i=1,NSPECIES
   ! print*, 'y(',i,')=',y(i)
    ys(i) = min(1.0e0, max(ys(i),1.0e-30))
   ! print *, 'y(',i,')=', y(i) 
enddo
  
!  y(iELEC) = y(iELEC)
!  print *, 'y(iELEC) = ', y(iELEC)

ie = iELEC



dfdy(iHP,iHP)=  -ratraw(iR003)*ys(iH)-ratraw(iR005)*ys(iHM)-ratraw(iR007)*ys(iH2) &
              & -ratraw(iR013)*ys(ie)-ratraw(iR016)*ys(iHM)-ratraw(iR027)*ys(iHE) &
              & -ratraw(iR034)*ys(iD)-ratraw(iR041)*ys(iHD)-ratraw(iR042)*ys(iD) &
              & -ratraw(iR060)*ys(iDM)-ratraw(iR067)*ys(iDM)-ratraw(iR092)*ys(iHD) &
              & -ratraw(iR093)*ys(iHD)-ratraw(iR097)*ys(iD2)-ratraw(iR098)*ys(iD2) &
              & -ratraw(iR099)*ys(iD2)

dfdy(iHP,iH)=  -ratraw(iR003)*ys(iHP)+ratraw(iR004)*ys(iH2P)+ratraw(iR012)*ys(ie) &
             & +ratraw(iR026)*ys(iHEP)+ratraw(iR035)*ys(iDP)+ratraw(iR038)*ys(iHDP) &
             & +ratraw(iR087)*ys(iD2P)

dfdy(iHP,iHM)=  -ratraw(iR005)*ys(iHP)-ratraw(iR016)*ys(iHP)

dfdy(iHP,iH2P)=  ratraw(iR004)*ys(iH)+ratraw(iR082)*ys(iD)+ratraw(iR118)

dfdy(iHP,iH2)=  -ratraw(iR007)*ys(iHP)+ratraw(iR024)*ys(iHEP)+ratraw(iR039)*ys(iDP)

dfdy(iHP,iDP)=  ratraw(iR035)*ys(iH)+ratraw(iR039)*ys(iH2)+ratraw(iR095)*ys(iHD)

dfdy(iHP,iD)=  -ratraw(iR034)*ys(iHP)-ratraw(iR042)*ys(iHP)+ratraw(iR082)*ys(iH2P) &
              & +ratraw(iR085)*ys(iHDP)

dfdy(iHP,iDM)=  -ratraw(iR060)*ys(iHP)-ratraw(iR067)*ys(iHP)

dfdy(iHP,iHDP)=  ratraw(iR038)*ys(iH)+ratraw(iR085)*ys(iD)+ratraw(iR120)

dfdy(iHP,iHD)=  -ratraw(iR041)*ys(iHP)-ratraw(iR092)*ys(iHP)-ratraw(iR093)*ys(iHP) &
              & +ratraw(iR095)*ys(iDP)+ratraw(iR102)*ys(iHEP)

dfdy(iHP,iD2P)=  ratraw(iR087)*ys(iH)

dfdy(iHP,iD2)=  -ratraw(iR097)*ys(iHP)-ratraw(iR098)*ys(iHP)-ratraw(iR099)*ys(iHP)

dfdy(iHP,iHEP)=  ratraw(iR024)*ys(iH2)+ratraw(iR026)*ys(iH)+ratraw(iR102)*ys(iHD)

dfdy(iHP,iHE)=  -ratraw(iR027)*ys(iHP)

dfdy(iHP,iHEPP)= 0.0e0 

dfdy(iHP,iELEC)=  ratraw(iR012)*ys(iH)-ratraw(iR013)*ys(iHP)



dfdy(iH,iHP)=  -ratraw(iR003)*ys(iH)+ratraw(iR005)*ys(iHM)+ratraw(iR005)*ys(iHM) &
             & +ratraw(iR007)*ys(iH2)+ratraw(iR013)*ys(ie)+ratraw(iR027)*ys(iHE) &
             & +ratraw(iR034)*ys(iD)+ratraw(iR067)*ys(iDM)+ratraw(iR092)*ys(iHD) &
             & +ratraw(iR099)*ys(iD2)

dfdy(iH,iH)=  -ratraw(iR001)*ys(ie)-ratraw(iR002)*ys(iHM)-ratraw(iR003)*ys(iHP) &
            & -ratraw(iR004)*ys(iH2P)-ratraw(iR009)*ys(iH2)+ratraw(iR009)*ys(iH2) &
            & +ratraw(iR009)*ys(iH2)+ratraw(iR009)*ys(iH2)-ratraw(iR012)*ys(ie) &
            & -ratraw(iR015)*ys(iHM)+ratraw(iR015)*ys(iHM)+ratraw(iR015)*ys(iHM) &
            & -ratraw(iR026)*ys(iHEP)-ratraw(iR030)*ys(iH)*ys(iH)-ratraw(iR030)*ys(iH)*ys(iH) &
            & -ratraw(iR030)*ys(iH)*ys(iH)-ratraw(iR030)*ys(iH)*ys(iH) &
            & -ratraw(iR030)*ys(iH)*ys(iH)-ratraw(iR030)*ys(iH)*ys(iH) &
            & -ratraw(iR030)*ys(iH)*ys(iH)-ratraw(iR030)*ys(iH)*ys(iH) &
            & -ratraw(iR030)*ys(iH)*ys(iH)+ratraw(iR030)*ys(iH)*ys(iH) &
            & +ratraw(iR030)*ys(iH)*ys(iH)+ratraw(iR030)*ys(iH)*ys(iH) &
            & -ratraw(iR031)*ys(iH)*ys(iH2)-ratraw(iR031)*ys(iH)*ys(iH2) &
            & -ratraw(iR031)*ys(iH)*ys(iH2)-ratraw(iR031)*ys(iH)*ys(iH2) &
            & -ratraw(iR032)*ys(iH)*ys(iHE)-ratraw(iR032)*ys(iH)*ys(iHE) &
            & -ratraw(iR032)*ys(iH)*ys(iHE)-ratraw(iR032)*ys(iH)*ys(iHE) &
            & -ratraw(iR035)*ys(iDP)-ratraw(iR036)*ys(iD)-ratraw(iR038)*ys(iHDP) &
            & -ratraw(iR040)*ys(iHD)-ratraw(iR043)*ys(iDP)-ratraw(iR050)*ys(iHDP) &
            & -ratraw(iR052)*ys(iDM)-ratraw(iR055)*ys(iDM)-ratraw(iR064)*ys(iDM) &
            & +ratraw(iR064)*ys(iDM)-ratraw(iR083)*ys(iHDP)-ratraw(iR087)*ys(iD2P) &
            & -ratraw(iR088)*ys(iD2P)-ratraw(iR089)*ys(iD2P)-ratraw(iR107)*ys(iD2) &
            & -ratraw(iR108)*ys(iHD)+ratraw(iR108)*ys(iHD)+ratraw(iR108)*ys(iHD) &
            & -ratraw(iR112)*ys(iD2)+ratraw(iR112)*ys(iD2)

dfdy(iH,iHM)=  -ratraw(iR002)*ys(iH)+ratraw(iR005)*ys(iHP)+ratraw(iR005)*ys(iHP) &
             & +ratraw(iR014)*ys(ie)-ratraw(iR015)*ys(iH)+ratraw(iR015)*ys(iH) &
             & +ratraw(iR015)*ys(iH)+ratraw(iR021)*ys(iH2P)+ratraw(iR022)*ys(iH2P) &
             & +ratraw(iR022)*ys(iH2P)+ratraw(iR022)*ys(iH2P)+ratraw(iR028)*ys(iHEP) &
             & +ratraw(iR029)*ys(iHE)+ratraw(iR053)*ys(iD)+ratraw(iR066)*ys(iDP) &
             & +ratraw(iR071)*ys(iHDP)+ratraw(iR072)*ys(iHDP)+ratraw(iR072)*ys(iHDP) &
             & +ratraw(iR075)*ys(iD2P)+ratraw(iR076)*ys(iD2P)+ratraw(iR116)

dfdy(iH,iH2P)=  -ratraw(iR004)*ys(iH)+ratraw(iR006)*ys(ie)+ratraw(iR006)*ys(ie) &
              & +ratraw(iR021)*ys(iHM)+ratraw(iR022)*ys(iHM)+ratraw(iR022)*ys(iHM) &
              & +ratraw(iR022)*ys(iHM)+ratraw(iR048)*ys(iD)+ratraw(iR070)*ys(iDM) &
              & +ratraw(iR070)*ys(iDM)+ratraw(iR118)

dfdy(iH,iH2)=  ratraw(iR007)*ys(iHP)+ratraw(iR008)*ys(ie)+ratraw(iR008)*ys(ie) &
            & -ratraw(iR009)*ys(iH)+ratraw(iR009)*ys(iH)+ratraw(iR009)*ys(iH) &
            & +ratraw(iR009)*ys(iH)+ratraw(iR010)*ys(iH2)+ratraw(iR010)*ys(iH2) &
            & +ratraw(iR010)*ys(iH2)+ratraw(iR010)*ys(iH2)+ratraw(iR011)*ys(iHE) &
            & +ratraw(iR011)*ys(iHE)+ratraw(iR023)*ys(ie)+ratraw(iR024)*ys(iHEP) &
            & -ratraw(iR031)*ys(iH)*ys(iH)-ratraw(iR031)*ys(iH)*ys(iH) &
            & +ratraw(iR037)*ys(iD)+ratraw(iR091)*ys(iDP)+ratraw(iR109)*ys(iHD) &
            & +ratraw(iR122)+ratraw(iR122)

dfdy(iH,iDP)=  -ratraw(iR035)*ys(iH)-ratraw(iR043)*ys(iH)+ratraw(iR066)*ys(iHM) &
             & +ratraw(iR091)*ys(iH2)+ratraw(iR096)*ys(iHD)

dfdy(iH,iD)=  ratraw(iR034)*ys(iHP)-ratraw(iR036)*ys(iH)+ratraw(iR037)*ys(iH2) &
           & +ratraw(iR048)*ys(iH2P)+ratraw(iR053)*ys(iHM)+ratraw(iR084)*ys(iHDP) &
           & +ratraw(iR106)*ys(iHD)

dfdy(iH,iDM)=  -ratraw(iR052)*ys(iH)-ratraw(iR055)*ys(iH)-ratraw(iR064)*ys(iH) &
             & +ratraw(iR064)*ys(iH)+ratraw(iR067)*ys(iHP)+ratraw(iR070)*ys(iH2P) &
             & +ratraw(iR070)*ys(iH2P)+ratraw(iR074)*ys(iHDP)

dfdy(iH,iHDP)=  -ratraw(iR038)*ys(iH)+ratraw(iR044)*ys(ie)-ratraw(iR050)*ys(iH) &
              & +ratraw(iR071)*ys(iHM)+ratraw(iR072)*ys(iHM)+ratraw(iR072)*ys(iHM) &
              & +ratraw(iR074)*ys(iDM)-ratraw(iR083)*ys(iH)+ratraw(iR084)*ys(iD)+ratraw(iR119)

dfdy(iH,iHD)=  -ratraw(iR040)*ys(iH)+ratraw(iR057)*ys(ie)+ratraw(iR092)*ys(iHP) &
             & +ratraw(iR096)*ys(iDP)+ratraw(iR103)*ys(iHEP)+ratraw(iR106)*ys(iD) &
             & -ratraw(iR108)*ys(iH)+ratraw(iR108)*ys(iH)+ratraw(iR108)*ys(iH) &
             & +ratraw(iR109)*ys(iH2)+ratraw(iR110)*ys(iHE)+ratraw(iR111)*ys(ie)+ratraw(iR123)

dfdy(iH,iD2P)=  ratraw(iR075)*ys(iHM)+ratraw(iR076)*ys(iHM)-ratraw(iR087)*ys(iH) &
             & -ratraw(iR088)*ys(iH)-ratraw(iR089)*ys(iH)

dfdy(iH,iD2)=  ratraw(iR099)*ys(iHP)-ratraw(iR107)*ys(iH)-ratraw(iR112)*ys(iH)+ratraw(iR112)*ys(iH)

dfdy(iH,iHEP)=  ratraw(iR024)*ys(iH2)-ratraw(iR026)*ys(iH)+ratraw(iR028)*ys(iHM)+ratraw(iR103)*ys(iHD)

dfdy(iH,iHE)=  ratraw(iR011)*ys(iH2)+ratraw(iR011)*ys(iH2)+ratraw(iR027)*ys(iHP) &
            & +ratraw(iR029)*ys(iHM)-ratraw(iR032)*ys(iH)*ys(iH)-ratraw(iR032)*ys(iH)*ys(iH) &
            & +ratraw(iR110)*ys(iHD)

dfdy(iH,iHEPP)= 0.0e0  

dfdy(iH,iELEC)=  -ratraw(iR001)*ys(iH)+ratraw(iR006)*ys(iH2P)+ratraw(iR006)*ys(iH2P) &
               & +ratraw(iR008)*ys(iH2)+ratraw(iR008)*ys(iH2)-ratraw(iR012)*ys(iH) &
               & +ratraw(iR013)*ys(iHP)+ratraw(iR014)*ys(iHM)+ratraw(iR023)*ys(iH2) &
               & +ratraw(iR044)*ys(iHDP)+ratraw(iR057)*ys(iHD)+ratraw(iR111)*ys(iHD)



dfdy(iHM,iHP)=  -ratraw(iR005)*ys(iHM)-ratraw(iR016)*ys(iHM)

dfdy(iHM,iH)=  ratraw(iR001)*ys(ie)-ratraw(iR002)*ys(iHM)-ratraw(iR015)*ys(iHM) &
            & +ratraw(iR052)*ys(iDM)

dfdy(iHM,iHM)=  -ratraw(iR002)*ys(iH)-ratraw(iR005)*ys(iHP)-ratraw(iR014)*ys(ie) &
              & -ratraw(iR015)*ys(iH)-ratraw(iR016)*ys(iHP)-ratraw(iR021)*ys(iH2P) &
              & -ratraw(iR022)*ys(iH2P)-ratraw(iR028)*ys(iHEP)-ratraw(iR029)*ys(iHE) &
              & -ratraw(iR053)*ys(iD)-ratraw(iR054)*ys(iD)-ratraw(iR061)*ys(iDP) &
              & -ratraw(iR066)*ys(iDP)-ratraw(iR071)*ys(iHDP)-ratraw(iR072)*ys(iHDP) &
              & -ratraw(iR075)*ys(iD2P)-ratraw(iR076)*ys(iD2P)-ratraw(iR116)

dfdy(iHM,iH2P)=  -ratraw(iR021)*ys(iHM)-ratraw(iR022)*ys(iHM)

dfdy(iHM,iH2)=  ratraw(iR023)*ys(ie)

dfdy(iHM,iDP)=  -ratraw(iR061)*ys(iHM)-ratraw(iR066)*ys(iHM)

dfdy(iHM,iD)=  -ratraw(iR053)*ys(iHM)-ratraw(iR054)*ys(iHM)

dfdy(iHM,iDM)=  ratraw(iR052)*ys(iH)

dfdy(iHM,iHDP)=  -ratraw(iR071)*ys(iHM)-ratraw(iR072)*ys(iHM)

dfdy(iHM,iHD)=  ratraw(iR058)*ys(ie)

dfdy(iHM,iD2P)=  -ratraw(iR075)*ys(iHM)-ratraw(iR076)*ys(iHM)

dfdy(iHM,iD2)=  0.0e0

dfdy(iHM,iHEP)=  -ratraw(iR028)*ys(iHM)

dfdy(iHM,iHE)=  -ratraw(iR029)*ys(iHM)

dfdy(iHM,iHEPP)=  0.0e0

dfdy(iHM,iELEC)=  ratraw(iR001)*ys(iH)-ratraw(iR014)*ys(iHM) &
               & +ratraw(iR023)*ys(iH2)+ratraw(iR058)*ys(iHD)



dfdy(iH2P,iHP)=  ratraw(iR003)*ys(iH)+ratraw(iR007)*ys(iH2)+ratraw(iR016)*ys(iHM) &
              & +ratraw(iR093)*ys(iHD)

dfdy(iH2P,iH)=  ratraw(iR003)*ys(iHP)-ratraw(iR004)*ys(iH2P)+ratraw(iR050)*ys(iHDP)

dfdy(iH2P,iHM)=  ratraw(iR016)*ys(iHP)-ratraw(iR021)*ys(iH2P)-ratraw(iR022)*ys(iH2P)

dfdy(iH2P,iH2P)=  -ratraw(iR004)*ys(iH)-ratraw(iR006)*ys(ie)-ratraw(iR021)*ys(iHM) &
                & -ratraw(iR022)*ys(iHM)-ratraw(iR048)*ys(iD)-ratraw(iR069)*ys(iDM) &
                & -ratraw(iR070)*ys(iDM)-ratraw(iR081)*ys(iD)-ratraw(iR082)*ys(iD) &
                & -ratraw(iR118)

dfdy(iH2P,iH2)=  ratraw(iR007)*ys(iHP)+ratraw(iR025)*ys(iHEP)+ratraw(iR090)*ys(iDP)

dfdy(iH2P,iDP)=  ratraw(iR090)*ys(iH2)

dfdy(iH2P,iD)=  -ratraw(iR048)*ys(iH2P)-ratraw(iR081)*ys(iH2P)-ratraw(iR082)*ys(iH2P)

dfdy(iH2P,iDM)=  -ratraw(iR069)*ys(iH2P)-ratraw(iR070)*ys(iH2P)

dfdy(iH2P,iHDP)=  ratraw(iR050)*ys(iH)

dfdy(iH2P,iHD)=  ratraw(iR093)*ys(iHP)

dfdy(iH2P,iD2P)= 0.0e0 

dfdy(iH2P,iD2)=  0.0e0

dfdy(iH2P,iHEP)=  ratraw(iR025)*ys(iH2)

dfdy(iH2P,iHE)=  0.0e0

dfdy(iH2P,iHEPP)=  0.0e0

dfdy(iH2P,iELEC)=  -ratraw(iR006)*ys(iH2P)



dfdy(iH2,iHP)=  -ratraw(iR007)*ys(iH2)+ratraw(iR041)*ys(iHD)

dfdy(iH2,iH)=  ratraw(iR002)*ys(iHM)+ratraw(iR004)*ys(iH2P)-ratraw(iR009)*ys(iH2) &
            & +ratraw(iR030)*ys(iH)*ys(iH)+ratraw(iR030)*ys(iH)*ys(iH) &
            & +ratraw(iR030)*ys(iH)*ys(iH)-ratraw(iR031)*ys(iH)*ys(iH2) &
            & -ratraw(iR031)*ys(iH)*ys(iH2)+ratraw(iR031)*ys(iH)*ys(iH2) &
            & +ratraw(iR031)*ys(iH)*ys(iH2)+ratraw(iR031)*ys(iH)*ys(iH2) &
            & +ratraw(iR031)*ys(iH)*ys(iH2)+ratraw(iR032)*ys(iH)*ys(iHE) &
            & +ratraw(iR032)*ys(iH)*ys(iHE)+ratraw(iR040)*ys(iHD)+ratraw(iR083)*ys(iHDP)

dfdy(iH2,iHM)=  ratraw(iR002)*ys(iH)+ratraw(iR021)*ys(iH2P)

dfdy(iH2,iH2P)=  ratraw(iR004)*ys(iH)+ratraw(iR021)*ys(iHM) &
              & +ratraw(iR069)*ys(iDM)+ratraw(iR081)*ys(iD)

dfdy(iH2,iH2)=  -ratraw(iR007)*ys(iHP)-ratraw(iR008)*ys(ie)-ratraw(iR009)*ys(iH) &
              & -ratraw(iR010)*ys(iH2)-ratraw(iR010)*ys(iH2)-ratraw(iR010)*ys(iH2) &
              & -ratraw(iR010)*ys(iH2)+ratraw(iR010)*ys(iH2)+ratraw(iR010)*ys(iH2) &
              & -ratraw(iR011)*ys(iHE)-ratraw(iR023)*ys(ie)-ratraw(iR024)*ys(iHEP) &
              & -ratraw(iR025)*ys(iHEP)-ratraw(iR031)*ys(iH)*ys(iH)+ratraw(iR031)*ys(iH)*ys(iH) &
              & +ratraw(iR031)*ys(iH)*ys(iH)-ratraw(iR037)*ys(iD)-ratraw(iR039)*ys(iDP) &
              & -ratraw(iR090)*ys(iDP)-ratraw(iR091)*ys(iDP)-ratraw(iR109)*ys(iHD) &
              & +ratraw(iR109)*ys(iHD)-ratraw(iR113)*ys(iD2)+ratraw(iR113)*ys(iD2)-ratraw(iR122)

dfdy(iH2,iDP)=  -ratraw(iR039)*ys(iH2)-ratraw(iR090)*ys(iH2)-ratraw(iR091)*ys(iH2)

dfdy(iH2,iD)=  -ratraw(iR037)*ys(iH2)+ratraw(iR081)*ys(iH2P)

dfdy(iH2,iDM)=  ratraw(iR069)*ys(iH2P)

dfdy(iH2,iHDP)=  ratraw(iR083)*ys(iH)

dfdy(iH2,iHD)=  ratraw(iR040)*ys(iH)+ratraw(iR041)*ys(iHP)-ratraw(iR109)*ys(iH2)+ratraw(iR109)*ys(iH2)

dfdy(iH2,iD2P)= 0.0e0 

dfdy(iH2,iD2)=  -ratraw(iR113)*ys(iH2)+ratraw(iR113)*ys(iH2)

!print *, 'Next line should be ZERO!!!'
!print *, 'dfdy(ih2,id2): ', dfdy(iH2,iD2)

dfdy(iH2,iHEP)=  -ratraw(iR024)*ys(iH2)-ratraw(iR025)*ys(iH2)

dfdy(iH2,iHE)=  -ratraw(iR011)*ys(iH2)+ratraw(iR032)*ys(iH)*ys(iH)

dfdy(iH2,iHEPP)=  0.0e0

dfdy(iH2,iELEC)=  -ratraw(iR008)*ys(iH2)-ratraw(iR023)*ys(iH2)



dfdy(iDP,iHP)=  ratraw(iR034)*ys(iD)+ratraw(iR041)*ys(iHD)+ratraw(iR097)*ys(iD2)

dfdy(iDP,iH)=  -ratraw(iR035)*ys(iDP)-ratraw(iR043)*ys(iDP)+ratraw(iR083)*ys(iHDP)+ratraw(iR089)*ys(iD2P)

dfdy(iDP,iHM)=  -ratraw(iR061)*ys(iDP)-ratraw(iR066)*ys(iDP)

dfdy(iDP,iH2P)=  ratraw(iR081)*ys(iD)

dfdy(iDP,iH2)=  -ratraw(iR039)*ys(iDP)-ratraw(iR090)*ys(iDP)-ratraw(iR091)*ys(iDP)

dfdy(iDP,iDP)=  -ratraw(iR033)*ys(ie)-ratraw(iR035)*ys(iH)-ratraw(iR039)*ys(iH2) &
              & -ratraw(iR043)*ys(iH)-ratraw(iR047)*ys(iHE)-ratraw(iR061)*ys(iHM) &
              & -ratraw(iR062)*ys(iDM)-ratraw(iR066)*ys(iHM)-ratraw(iR068)*ys(iDM) &
              & -ratraw(iR080)*ys(iD)-ratraw(iR090)*ys(iH2)-ratraw(iR091)*ys(iH2) &
              & -ratraw(iR094)*ys(iHD)-ratraw(iR095)*ys(iHD)-ratraw(iR096)*ys(iHD) &
              & -ratraw(iR100)*ys(iD2)

dfdy(iDP,iD)=  ratraw(iR034)*ys(iHP)+ratraw(iR045)*ys(ie)+ratraw(iR046)*ys(iHEP) &
            & +ratraw(iR049)*ys(iHDP)-ratraw(iR080)*ys(iDP)+ratraw(iR081)*ys(iH2P) &
            & +ratraw(iR086)*ys(iD2P)

dfdy(iDP,iDM)=  -ratraw(iR062)*ys(iDP)-ratraw(iR068)*ys(iDP)

dfdy(iDP,iHDP)=  ratraw(iR049)*ys(iD)+ratraw(iR083)*ys(iH)+ratraw(iR119)

dfdy(iDP,iHD)=  ratraw(iR041)*ys(iHP)-ratraw(iR094)*ys(iDP)-ratraw(iR095)*ys(iDP) &
             & -ratraw(iR096)*ys(iDP)+ratraw(iR103)*ys(iHEP)

dfdy(iDP,iD2P)=  ratraw(iR086)*ys(iD)+ratraw(iR089)*ys(iH)+ratraw(iR121)

dfdy(iDP,iD2)=  ratraw(iR097)*ys(iHP)-ratraw(iR100)*ys(iDP)+ratraw(iR105)*ys(iHEP)

dfdy(iDP,iHEP)=  ratraw(iR046)*ys(iD)+ratraw(iR103)*ys(iHD)+ratraw(iR105)*ys(iD2)

dfdy(iDP,iHE)=  -ratraw(iR047)*ys(iDP)

dfdy(iDP,iHEPP)= 0.0e0 

dfdy(iDP,iELEC)=  -ratraw(iR033)*ys(iDP)+ratraw(iR045)*ys(iD)



dfdy(iD,iHP)=  -ratraw(iR034)*ys(iD)-ratraw(iR042)*ys(iD)+ratraw(iR067)*ys(iDM) &
             & +ratraw(iR093)*ys(iHD)+ratraw(iR098)*ys(iD2)

dfdy(iD,iH)=  ratraw(iR035)*ys(iDP)-ratraw(iR036)*ys(iD)+ratraw(iR040)*ys(iHD) &
           & +ratraw(iR050)*ys(iHDP)+ratraw(iR052)*ys(iDM)+ratraw(iR064)*ys(iDM) &
           & +ratraw(iR088)*ys(iD2P)+ratraw(iR107)*ys(iD2)+ratraw(iR108)*ys(iHD) &
           & +ratraw(iR112)*ys(iD2)+ratraw(iR112)*ys(iD2)

dfdy(iD,iHM)=  -ratraw(iR053)*ys(iD)-ratraw(iR054)*ys(iD)+ratraw(iR066)*ys(iDP) &
             & +ratraw(iR072)*ys(iHDP)+ratraw(iR076)*ys(iD2P)+ratraw(iR076)*ys(iD2P)

dfdy(iD,iH2P)=  -ratraw(iR048)*ys(iD)+ratraw(iR069)*ys(iDM)+ratraw(iR070)*ys(iDM) &
              & -ratraw(iR081)*ys(iD)-ratraw(iR082)*ys(iD)

dfdy(iD,iH2)=  -ratraw(iR037)*ys(iD)+ratraw(iR090)*ys(iDP)+ratraw(iR109)*ys(iHD) &
             & +ratraw(iR113)*ys(iD2)+ratraw(iR113)*ys(iD2)

dfdy(iD,iDP)=  ratraw(iR033)*ys(ie)+ratraw(iR035)*ys(iH)+ratraw(iR047)*ys(iHE) &
            & +ratraw(iR066)*ys(iHM)+ratraw(iR068)*ys(iDM)+ratraw(iR068)*ys(iDM) &
            & -ratraw(iR080)*ys(iD)+ratraw(iR090)*ys(iH2)+ratraw(iR094)*ys(iHD) &
            & +ratraw(iR100)*ys(iD2)

dfdy(iD,iD)=  -ratraw(iR034)*ys(iHP)-ratraw(iR036)*ys(iH)-ratraw(iR037)*ys(iH2) &
            & -ratraw(iR042)*ys(iHP)-ratraw(iR045)*ys(ie)-ratraw(iR046)*ys(iHEP) &
            & -ratraw(iR048)*ys(iH2P)-ratraw(iR049)*ys(iHDP)-ratraw(iR051)*ys(ie) &
            & -ratraw(iR053)*ys(iHM)-ratraw(iR054)*ys(iHM)-ratraw(iR056)*ys(iDM) &
            & -ratraw(iR080)*ys(iDP)-ratraw(iR081)*ys(iH2P)-ratraw(iR082)*ys(iH2P) &
            & -ratraw(iR084)*ys(iHDP)-ratraw(iR085)*ys(iHDP)-ratraw(iR086)*ys(iD2P) &
            & -ratraw(iR106)*ys(iHD)

dfdy(iD,iDM)=  ratraw(iR052)*ys(iH)-ratraw(iR056)*ys(iD)+ratraw(iR063)*ys(ie) &
            & +ratraw(iR064)*ys(iH)+ratraw(iR065)*ys(iHE)+ratraw(iR067)*ys(iHP) &
            & +ratraw(iR068)*ys(iDP)+ratraw(iR068)*ys(iDP)+ratraw(iR069)*ys(iH2P) &
            & +ratraw(iR070)*ys(iH2P)+ratraw(iR073)*ys(iHDP)+ratraw(iR074)*ys(iHDP) &
            & +ratraw(iR074)*ys(iHDP)+ratraw(iR077)*ys(iD2P)+ratraw(iR078)*ys(iD2P) &
            & +ratraw(iR078)*ys(iD2P)+ratraw(iR078)*ys(iD2P)+ratraw(iR079)*ys(iHEP) &
            & +ratraw(iR117)

dfdy(iD,iHDP)=  ratraw(iR044)*ys(ie)-ratraw(iR049)*ys(iD)+ratraw(iR050)*ys(iH) &
             & +ratraw(iR072)*ys(iHM)+ratraw(iR073)*ys(iDM)+ratraw(iR074)*ys(iDM) &
             & +ratraw(iR074)*ys(iDM)-ratraw(iR084)*ys(iD)-ratraw(iR085)*ys(iD)+ratraw(iR120)

dfdy(iD,iHD)=  ratraw(iR040)*ys(iH)+ratraw(iR058)*ys(ie)+ratraw(iR093)*ys(iHP) &
            & +ratraw(iR094)*ys(iDP)+ratraw(iR102)*ys(iHEP)-ratraw(iR106)*ys(iD) &
            & +ratraw(iR108)*ys(iH)+ratraw(iR109)*ys(iH2)+ratraw(iR110)*ys(iHE) &
            & +ratraw(iR111)*ys(ie)+ratraw(iR123)

dfdy(iD,iD2P)=  ratraw(iR076)*ys(iHM)+ratraw(iR076)*ys(iHM)+ratraw(iR077)*ys(iDM) &
             & +ratraw(iR078)*ys(iDM)+ratraw(iR078)*ys(iDM)+ratraw(iR078)*ys(iDM) &
             & -ratraw(iR086)*ys(iD)+ratraw(iR088)*ys(iH)+ratraw(iR121)

dfdy(iD,iD2)=  ratraw(iR059)*ys(ie)+ratraw(iR098)*ys(iHP)+ratraw(iR100)*ys(iDP) &
            & +ratraw(iR105)*ys(iHEP)+ratraw(iR107)*ys(iH)+ratraw(iR112)*ys(iH) &
            & +ratraw(iR112)*ys(iH)+ratraw(iR113)*ys(iH2)+ratraw(iR113)*ys(iH2) &
            & +ratraw(iR114)*ys(iHE)+ratraw(iR114)*ys(iHE)+ratraw(iR115)*ys(ie) &
            & +ratraw(iR115)*ys(ie)+ratraw(iR124)+ratraw(iR124)

dfdy(iD,iHEP)=  -ratraw(iR046)*ys(iD)+ratraw(iR079)*ys(iDM) &
              & +ratraw(iR102)*ys(iHD)+ratraw(iR105)*ys(iD2)

dfdy(iD,iHE)=  ratraw(iR047)*ys(iDP)+ratraw(iR065)*ys(iDM)+ratraw(iR110)*ys(iHD) &
            & +ratraw(iR114)*ys(iD2)+ratraw(iR114)*ys(iD2)

dfdy(iD,iHEPP)= 0.0e0  

dfdy(iD,iELEC)=  ratraw(iR033)*ys(iDP)+ratraw(iR044)*ys(iHDP)-ratraw(iR045)*ys(iD) &
              & -ratraw(iR051)*ys(iD)+ratraw(iR058)*ys(iHD)+ratraw(iR059)*ys(iD2) &
              & +ratraw(iR063)*ys(iDM)+ratraw(iR111)*ys(iHD)+ratraw(iR115)*ys(iD2) &
              & +ratraw(iR115)*ys(iD2)



dfdy(iDM,iHP)=  -ratraw(iR060)*ys(iDM)-ratraw(iR067)*ys(iDM)

dfdy(iDM,iH)=  -ratraw(iR052)*ys(iDM)-ratraw(iR055)*ys(iDM)-ratraw(iR064)*ys(iDM)

dfdy(iDM,iHM)=  ratraw(iR053)*ys(iD)

dfdy(iDM,iH2P)=  -ratraw(iR069)*ys(iDM)-ratraw(iR070)*ys(iDM)

dfdy(iDM,iH2)=  0.0e0

dfdy(iDM,iDP)=  -ratraw(iR062)*ys(iDM)-ratraw(iR068)*ys(iDM)

dfdy(iDM,iD)=  ratraw(iR051)*ys(ie)+ratraw(iR053)*ys(iHM)-ratraw(iR056)*ys(iDM)

dfdy(iDM,iDM)=  -ratraw(iR052)*ys(iH)-ratraw(iR055)*ys(iH)-ratraw(iR056)*ys(iD) &
              & -ratraw(iR060)*ys(iHP)-ratraw(iR062)*ys(iDP)-ratraw(iR063)*ys(ie) &
              & -ratraw(iR064)*ys(iH)-ratraw(iR065)*ys(iHE)-ratraw(iR067)*ys(iHP) &
              & -ratraw(iR068)*ys(iDP)-ratraw(iR069)*ys(iH2P)-ratraw(iR070)*ys(iH2P) &
              & -ratraw(iR073)*ys(iHDP)-ratraw(iR074)*ys(iHDP)-ratraw(iR077)*ys(iD2P) &
              & -ratraw(iR078)*ys(iD2P)-ratraw(iR079)*ys(iHEP)-ratraw(iR117)

dfdy(iDM,iHDP)=  -ratraw(iR073)*ys(iDM)-ratraw(iR074)*ys(iDM)

dfdy(iDM,iHD)=  ratraw(iR057)*ys(ie)

dfdy(iDM,iD2P)=  -ratraw(iR077)*ys(iDM)-ratraw(iR078)*ys(iDM)

dfdy(iDM,iD2)=  ratraw(iR059)*ys(ie)

dfdy(iDM,iHEP)=  -ratraw(iR079)*ys(iDM)

dfdy(iDM,iHE)=  -ratraw(iR065)*ys(iDM)

dfdy(iDM,iHEPP)= 0.0e0  

dfdy(iDM,iELEC)=  ratraw(iR051)*ys(iD)+ratraw(iR057)*ys(iHD)+ratraw(iR059)*ys(iD2)-ratraw(iR063)*ys(iDM)



dfdy(iHDP,iHP)=  ratraw(iR042)*ys(iD)+ratraw(iR060)*ys(iDM)+ratraw(iR092)*ys(iHD)+ratraw(iR098)*ys(iD2)

dfdy(iHDP,iH)=  -ratraw(iR038)*ys(iHDP)+ratraw(iR043)*ys(iDP)-ratraw(iR050)*ys(iHDP) &
              & -ratraw(iR083)*ys(iHDP)+ratraw(iR088)*ys(iD2P)

dfdy(iHDP,iHM)=  ratraw(iR061)*ys(iDP)-ratraw(iR071)*ys(iHDP)-ratraw(iR072)*ys(iHDP)

dfdy(iHDP,iH2P)=  ratraw(iR048)*ys(iD)

dfdy(iHDP,iH2)=  ratraw(iR091)*ys(iDP)

dfdy(iHDP,iDP)=  ratraw(iR043)*ys(iH)+ratraw(iR061)*ys(iHM)+ratraw(iR091)*ys(iH2)+ratraw(iR094)*ys(iHD)

dfdy(iHDP,iD)=  ratraw(iR042)*ys(iHP)+ratraw(iR048)*ys(iH2P)-ratraw(iR049)*ys(iHDP) &
             & -ratraw(iR084)*ys(iHDP)-ratraw(iR085)*ys(iHDP)

dfdy(iHDP,iDM)=  ratraw(iR060)*ys(iHP)-ratraw(iR073)*ys(iHDP)-ratraw(iR074)*ys(iHDP)

dfdy(iHDP,iHDP)=  -ratraw(iR038)*ys(iH)-ratraw(iR044)*ys(ie)-ratraw(iR049)*ys(iD) &
                & -ratraw(iR050)*ys(iH)-ratraw(iR071)*ys(iHM)-ratraw(iR072)*ys(iHM) &
                & -ratraw(iR073)*ys(iDM)-ratraw(iR074)*ys(iDM)-ratraw(iR083)*ys(iH) &
                & -ratraw(iR084)*ys(iD)-ratraw(iR085)*ys(iD)-ratraw(iR119)-ratraw(iR120)

dfdy(iHDP,iHD)=  ratraw(iR092)*ys(iHP)+ratraw(iR094)*ys(iDP)+ratraw(iR101)*ys(iHEP)

dfdy(iHDP,iD2P)=  ratraw(iR088)*ys(iH)

dfdy(iHDP,iD2)=  ratraw(iR098)*ys(iHP)

dfdy(iHDP,iHEP)=  ratraw(iR101)*ys(iHD)

dfdy(iHDP,iHE)=  0.0e0

dfdy(iHDP,iHEPP)=  0.0e0

dfdy(iHDP,iELEC)=  -ratraw(iR044)*ys(iHDP)



dfdy(iHD,iHP)=  -ratraw(iR041)*ys(iHD)-ratraw(iR092)*ys(iHD)-ratraw(iR093)*ys(iHD)+ratraw(iR097)*ys(iD2)

dfdy(iHD,iH)=  ratraw(iR036)*ys(iD)+ratraw(iR038)*ys(iHDP)-ratraw(iR040)*ys(iHD) &
            & +ratraw(iR055)*ys(iDM)+ratraw(iR089)*ys(iD2P)+ratraw(iR107)*ys(iD2) &
            & -ratraw(iR108)*ys(iHD)

dfdy(iHD,iHM)=  ratraw(iR054)*ys(iD)+ratraw(iR071)*ys(iHDP)

dfdy(iHD,iH2P)=  ratraw(iR082)*ys(iD)

dfdy(iHD,iH2)=  ratraw(iR037)*ys(iD)+ratraw(iR039)*ys(iDP)-ratraw(iR109)*ys(iHD)

dfdy(iHD,iDP)=  ratraw(iR039)*ys(iH2)-ratraw(iR094)*ys(iHD)-ratraw(iR095)*ys(iHD)-ratraw(iR096)*ys(iHD)

dfdy(iHD,iD)=  ratraw(iR036)*ys(iH)+ratraw(iR037)*ys(iH2)+ratraw(iR049)*ys(iHDP) &
            & +ratraw(iR054)*ys(iHM)+ratraw(iR082)*ys(iH2P)-ratraw(iR106)*ys(iHD)

dfdy(iHD,iDM)=  ratraw(iR055)*ys(iH)+ratraw(iR073)*ys(iHDP)

dfdy(iHD,iHDP)=  ratraw(iR038)*ys(iH)+ratraw(iR049)*ys(iD)+ratraw(iR071)*ys(iHM)+ratraw(iR073)*ys(iDM)

dfdy(iHD,iHD)=  -ratraw(iR040)*ys(iH)-ratraw(iR041)*ys(iHP)-ratraw(iR057)*ys(ie) &
              & -ratraw(iR058)*ys(ie)-ratraw(iR092)*ys(iHP)-ratraw(iR093)*ys(iHP) &
              & -ratraw(iR094)*ys(iDP)-ratraw(iR095)*ys(iDP)-ratraw(iR096)*ys(iDP) &
              & -ratraw(iR101)*ys(iHEP)-ratraw(iR102)*ys(iHEP)-ratraw(iR103)*ys(iHEP) &
              & -ratraw(iR106)*ys(iD)-ratraw(iR108)*ys(iH)-ratraw(iR109)*ys(iH2) &
              & -ratraw(iR110)*ys(iHE)-ratraw(iR111)*ys(ie)-ratraw(iR123)

dfdy(iHD,iD2P)=  ratraw(iR089)*ys(iH)

dfdy(iHD,iD2)=  ratraw(iR097)*ys(iHP)+ratraw(iR107)*ys(iH)

dfdy(iHD,iHEP)=  -ratraw(iR101)*ys(iHD)-ratraw(iR102)*ys(iHD)-ratraw(iR103)*ys(iHD)

dfdy(iHD,iHE)=  -ratraw(iR110)*ys(iHD)

dfdy(iHD,iHEPP)= 0.0e0 

dfdy(iHD,iELEC)=  -ratraw(iR057)*ys(iHD)-ratraw(iR058)*ys(iHD)-ratraw(iR111)*ys(iHD)



dfdy(iD2P,iHP)=  ratraw(iR099)*ys(iD2)

dfdy(iD2P,iH)=  -ratraw(iR087)*ys(iD2P)-ratraw(iR088)*ys(iD2P)-ratraw(iR089)*ys(iD2P)

dfdy(iD2P,iHM)=  -ratraw(iR075)*ys(iD2P)-ratraw(iR076)*ys(iD2P)

dfdy(iD2P,iH2P)= 0.0e0 

dfdy(iD2P,iH2)=  0.0e0

dfdy(iD2P,iDP)=  ratraw(iR062)*ys(iDM)+ratraw(iR080)*ys(iD)+ratraw(iR096)*ys(iHD)+ratraw(iR100)*ys(iD2)

dfdy(iD2P,iD)=  ratraw(iR080)*ys(iDP)+ratraw(iR084)*ys(iHDP)-ratraw(iR086)*ys(iD2P)

dfdy(iD2P,iDM)=  ratraw(iR062)*ys(iDP)-ratraw(iR077)*ys(iD2P)-ratraw(iR078)*ys(iD2P)

dfdy(iD2P,iHDP)=  ratraw(iR084)*ys(iD)

dfdy(iD2P,iHD)=  ratraw(iR096)*ys(iDP)

dfdy(iD2P,iD2P)=  -ratraw(iR075)*ys(iHM)-ratraw(iR076)*ys(iHM)-ratraw(iR077)*ys(iDM) &
                & -ratraw(iR078)*ys(iDM)-ratraw(iR086)*ys(iD)-ratraw(iR087)*ys(iH) &
                & -ratraw(iR088)*ys(iH)-ratraw(iR089)*ys(iH)-ratraw(iR121)

dfdy(iD2P,iD2)=  ratraw(iR099)*ys(iHP)+ratraw(iR100)*ys(iDP)+ratraw(iR104)*ys(iHEP)

dfdy(iD2P,iHEP)=  ratraw(iR104)*ys(iD2)

dfdy(iD2P,iHE)=  0.0e0

dfdy(iD2P,iHEPP)=  0.0e0

dfdy(iD2P,iELEC)=  0.0e0



dfdy(iD2,iHP)=  -ratraw(iR097)*ys(iD2)-ratraw(iR098)*ys(iD2)-ratraw(iR099)*ys(iD2)

dfdy(iD2,iH)=  ratraw(iR087)*ys(iD2P)-ratraw(iR107)*ys(iD2)-ratraw(iR112)*ys(iD2)

dfdy(iD2,iHM)=  ratraw(iR075)*ys(iD2P)

dfdy(iD2,iH2P)= 0.0e0 

dfdy(iD2,iH2)=  -ratraw(iR113)*ys(iD2)

dfdy(iD2,iDP)=  ratraw(iR095)*ys(iHD)-ratraw(iR100)*ys(iD2)

dfdy(iD2,iD)=  ratraw(iR056)*ys(iDM)+ratraw(iR085)*ys(iHDP)+ratraw(iR086)*ys(iD2P)+ratraw(iR106)*ys(iHD)

dfdy(iD2,iDM)=  ratraw(iR056)*ys(iD)+ratraw(iR077)*ys(iD2P)

dfdy(iD2,iHDP)=  ratraw(iR085)*ys(iD)

dfdy(iD2,iHD)=  ratraw(iR095)*ys(iDP)+ratraw(iR106)*ys(iD)

dfdy(iD2,iD2P)=  ratraw(iR075)*ys(iHM)+ratraw(iR077)*ys(iDM)+ratraw(iR086)*ys(iD)+ratraw(iR087)*ys(iH)

dfdy(iD2,iD2)=  -ratraw(iR059)*ys(ie)-ratraw(iR097)*ys(iHP)-ratraw(iR098)*ys(iHP) &
              & -ratraw(iR099)*ys(iHP)-ratraw(iR100)*ys(iDP)-ratraw(iR104)*ys(iHEP) &
              & -ratraw(iR105)*ys(iHEP)-ratraw(iR107)*ys(iH)-ratraw(iR112)*ys(iH) &
              & -ratraw(iR113)*ys(iH2)-ratraw(iR114)*ys(iHE)-ratraw(iR115)*ys(ie)-ratraw(iR124)

dfdy(iD2,iHEP)=  -ratraw(iR104)*ys(iD2)-ratraw(iR105)*ys(iD2)

dfdy(iD2,iHE)=  -ratraw(iR114)*ys(iD2)

dfdy(iD2,iHEPP)= 0.0e0 

dfdy(iD2,iELEC)=  -ratraw(iR059)*ys(iD2)-ratraw(iR115)*ys(iD2)



dfdy(iHEP,iHP)=  ratraw(iR027)*ys(iHE)

dfdy(iHEP,iH)=  -ratraw(iR026)*ys(iHEP)

dfdy(iHEP,iHM)=  -ratraw(iR028)*ys(iHEP)

dfdy(iHEP,iH2P)= 0.0e0 

dfdy(iHEP,iH2)=  -ratraw(iR024)*ys(iHEP)-ratraw(iR025)*ys(iHEP)

dfdy(iHEP,iDP)=  ratraw(iR047)*ys(iHE)

dfdy(iHEP,iD)=  -ratraw(iR046)*ys(iHEP)

dfdy(iHEP,iDM)=  -ratraw(iR079)*ys(iHEP)

dfdy(iHEP,iHDP)= 0.0e0 

dfdy(iHEP,iHD)=  -ratraw(iR101)*ys(iHEP)-ratraw(iR102)*ys(iHEP)-ratraw(iR103)*ys(iHEP)

dfdy(iHEP,iD2P)=  0.0e0

dfdy(iHEP,iD2)=  -ratraw(iR104)*ys(iHEP)-ratraw(iR105)*ys(iHEP)

dfdy(iHEP,iHEP)=  -ratraw(iR018)*ys(ie)-ratraw(iR019)*ys(ie)-ratraw(iR024)*ys(iH2) &
                & -ratraw(iR025)*ys(iH2)-ratraw(iR026)*ys(iH)-ratraw(iR028)*ys(iHM) &
                & -ratraw(iR046)*ys(iD)-ratraw(iR079)*ys(iDM)-ratraw(iR101)*ys(iHD) &
                & -ratraw(iR102)*ys(iHD)-ratraw(iR103)*ys(iHD)-ratraw(iR104)*ys(iD2) &
                & -ratraw(iR105)*ys(iD2)

dfdy(iHEP,iHE)=  ratraw(iR017)*ys(ie)+ratraw(iR027)*ys(iHP)+ratraw(iR047)*ys(iDP)

dfdy(iHEP,iHEPP)=  ratraw(iR020)*ys(ie)

dfdy(iHEP,iELEC)=  ratraw(iR017)*ys(iHE)-ratraw(iR018)*ys(iHEP)-ratraw(iR019)*ys(iHEP) &
                & +ratraw(iR020)*ys(iHEPP)



dfdy(iHE,iHP)=  -ratraw(iR027)*ys(iHE)

dfdy(iHE,iH)=  ratraw(iR026)*ys(iHEP)-ratraw(iR032)*ys(iH)*ys(iHE)-ratraw(iR032)*ys(iH)*ys(iHE) &
            & +ratraw(iR032)*ys(iH)*ys(iHE)+ratraw(iR032)*ys(iH)*ys(iHE)

dfdy(iHE,iHM)=  ratraw(iR028)*ys(iHEP)-ratraw(iR029)*ys(iHE)+ratraw(iR029)*ys(iHE)

dfdy(iHE,iH2P)= 0.0e0 

dfdy(iHE,iH2)=  -ratraw(iR011)*ys(iHE)+ratraw(iR011)*ys(iHE)+ratraw(iR024)*ys(iHEP)+ratraw(iR025)*ys(iHEP)

dfdy(iHE,iDP)=  -ratraw(iR047)*ys(iHE)

dfdy(iHE,iD)=  ratraw(iR046)*ys(iHEP)

dfdy(iHE,iDM)=  -ratraw(iR065)*ys(iHE)+ratraw(iR065)*ys(iHE)+ratraw(iR079)*ys(iHEP)

dfdy(iHE,iHDP)= 0.0e0  

dfdy(iHE,iHD)=  ratraw(iR101)*ys(iHEP)+ratraw(iR102)*ys(iHEP)+ratraw(iR103)*ys(iHEP) &
             & -ratraw(iR110)*ys(iHE)+ratraw(iR110)*ys(iHE)

dfdy(iHE,iD2P)= 0.0e0 

dfdy(iHE,iD2)=  ratraw(iR104)*ys(iHEP)+ratraw(iR105)*ys(iHEP)-ratraw(iR114)*ys(iHE)+ratraw(iR114)*ys(iHE)

dfdy(iHE,iHEP)=  ratraw(iR019)*ys(ie)+ratraw(iR024)*ys(iH2)+ratraw(iR025)*ys(iH2) &
              & +ratraw(iR026)*ys(iH)+ratraw(iR028)*ys(iHM)+ratraw(iR046)*ys(iD) &
              & +ratraw(iR079)*ys(iDM)+ratraw(iR101)*ys(iHD)+ratraw(iR102)*ys(iHD) &
              & +ratraw(iR103)*ys(iHD)+ratraw(iR104)*ys(iD2)+ratraw(iR105)*ys(iD2)

dfdy(iHE,iHE)=  -ratraw(iR011)*ys(iH2)+ratraw(iR011)*ys(iH2)-ratraw(iR017)*ys(ie) &
              & -ratraw(iR027)*ys(iHP)-ratraw(iR029)*ys(iHM)+ratraw(iR029)*ys(iHM) &
              & -ratraw(iR032)*ys(iH)*ys(iH)+ratraw(iR032)*ys(iH)*ys(iH) &
              & -ratraw(iR047)*ys(iDP)-ratraw(iR065)*ys(iDM)+ratraw(iR065)*ys(iDM) &
              & -ratraw(iR110)*ys(iHD)+ratraw(iR110)*ys(iHD)-ratraw(iR114)*ys(iD2) &
              & +ratraw(iR114)*ys(iD2)

dfdy(iHE,iHEPP)= 0.0e0 

dfdy(iHE,iELEC)=  -ratraw(iR017)*ys(iHE)+ratraw(iR019)*ys(iHEP)



dfdy(iHEPP,iHP)=  0.0e0

dfdy(iHEPP,iH)=  0.0e0

dfdy(iHEPP,iHM)=  0.0e0

dfdy(iHEPP,iH2P)=  0.0e0

dfdy(iHEPP,iH2)=  0.0e0

dfdy(iHEPP,iDP)=  0.0e0

dfdy(iHEPP,iD)=  0.0e0

dfdy(iHEPP,iDM)=  0.0e0

dfdy(iHEPP,iHDP)=  0.0e0

dfdy(iHEPP,iHD)=  0.0e0

dfdy(iHEPP,iD2P)=  0.0e0

dfdy(iHEPP,iD2)=  0.0e0

dfdy(iHEPP,iHEP)=  ratraw(iR018)*ys(ie)

dfdy(iHEPP,iHE)=  0.0e0

dfdy(iHEPP,iHEPP)=  -ratraw(iR020)*ys(ie)

dfdy(iHEPP,iELEC)=  ratraw(iR018)*ys(iHEP)-ratraw(iR020)*ys(iHEPP)



dfdy(iELEC,iHP)=  -ratraw(iR013)*ys(ie)+ratraw(iR016)*ys(iHM)+ratraw(iR060)*ys(iDM)

dfdy(iELEC,iH)=  -ratraw(iR001)*ys(ie)+ratraw(iR002)*ys(iHM)-ratraw(iR012)*ys(ie) &
               & +ratraw(iR012)*ys(ie)+ratraw(iR012)*ys(ie)+ratraw(iR015)*ys(iHM) &
               & +ratraw(iR055)*ys(iDM)+ratraw(iR064)*ys(iDM)

dfdy(iELEC,iHM)=  ratraw(iR002)*ys(iH)-ratraw(iR014)*ys(ie)+ratraw(iR014)*ys(ie) &
               & +ratraw(iR014)*ys(ie)+ratraw(iR015)*ys(iH)+ratraw(iR016)*ys(iHP) &
               & +ratraw(iR029)*ys(iHE)+ratraw(iR054)*ys(iD)+ratraw(iR061)*ys(iDP) &
               & +ratraw(iR116)

dfdy(iELEC,iH2P)=  -ratraw(iR006)*ys(ie)

dfdy(iELEC,iH2)=  -ratraw(iR008)*ys(ie)+ratraw(iR008)*ys(ie)-ratraw(iR023)*ys(ie)

dfdy(iELEC,iDP)=  -ratraw(iR033)*ys(ie)+ratraw(iR061)*ys(iHM)+ratraw(iR062)*ys(iDM)

dfdy(iELEC,iD)=  -ratraw(iR045)*ys(ie)+ratraw(iR045)*ys(ie)+ratraw(iR045)*ys(ie) &
               & -ratraw(iR051)*ys(ie)+ratraw(iR054)*ys(iHM)+ratraw(iR056)*ys(iDM)

dfdy(iELEC,iDM)=  ratraw(iR055)*ys(iH)+ratraw(iR056)*ys(iD)+ratraw(iR060)*ys(iHP) &
               & +ratraw(iR062)*ys(iDP)-ratraw(iR063)*ys(ie)+ratraw(iR063)*ys(ie) &
               & +ratraw(iR063)*ys(ie)+ratraw(iR064)*ys(iH)+ratraw(iR065)*ys(iHE) &
               & +ratraw(iR117)

dfdy(iELEC,iHDP)=  -ratraw(iR044)*ys(ie)

dfdy(iELEC,iHD)=  -ratraw(iR057)*ys(ie)-ratraw(iR058)*ys(ie)-ratraw(iR111)*ys(ie)+ratraw(iR111)*ys(ie)

dfdy(iELEC,iD2P)= 0.0e0 

dfdy(iELEC,iD2)=  -ratraw(iR059)*ys(ie)-ratraw(iR115)*ys(ie)+ratraw(iR115)*ys(ie)

dfdy(iELEC,iHEP)=  -ratraw(iR018)*ys(ie)+ratraw(iR018)*ys(ie)+ratraw(iR018)*ys(ie)-ratraw(iR019)*ys(ie)

dfdy(iELEC,iHE)=  -ratraw(iR017)*ys(ie)+ratraw(iR017)*ys(ie)+ratraw(iR017)*ys(ie) &
                & +ratraw(iR029)*ys(iHM)+ratraw(iR065)*ys(iDM)

dfdy(iELEC,iHEPP)=  -ratraw(iR020)*ys(ie)

dfdy(iELEC,iELEC)=  -ratraw(iR001)*ys(iH)-ratraw(iR006)*ys(iH2P)-ratraw(iR008)*ys(iH2) &
                  & +ratraw(iR008)*ys(iH2)-ratraw(iR012)*ys(iH)+ratraw(iR012)*ys(iH) &
                  & +ratraw(iR012)*ys(iH)-ratraw(iR013)*ys(iHP)-ratraw(iR014)*ys(iHM) &
                  & +ratraw(iR014)*ys(iHM)+ratraw(iR014)*ys(iHM)-ratraw(iR017)*ys(iHE) &
                  & +ratraw(iR017)*ys(iHE)+ratraw(iR017)*ys(iHE)-ratraw(iR018)*ys(iHEP) &
                  & +ratraw(iR018)*ys(iHEP)+ratraw(iR018)*ys(iHEP)-ratraw(iR019)*ys(iHEP) &
                  & -ratraw(iR020)*ys(iHEPP)-ratraw(iR023)*ys(iH2)-ratraw(iR033)*ys(iDP) &
                  & -ratraw(iR044)*ys(iHDP)-ratraw(iR045)*ys(iD)+ratraw(iR045)*ys(iD) &
                  & +ratraw(iR045)*ys(iD)-ratraw(iR051)*ys(iD)-ratraw(iR057)*ys(iHD) &
                  & -ratraw(iR058)*ys(iHD)-ratraw(iR059)*ys(iD2)-ratraw(iR063)*ys(iDM) &
                  & +ratraw(iR063)*ys(iDM)+ratraw(iR063)*ys(iDM)-ratraw(iR111)*ys(iHD) &
                  & +ratraw(iR111)*ys(iHD)-ratraw(iR115)*ys(iD2)+ratraw(iR115)*ys(iD2)



!   do i = 1,nlog
!    do j = 1,nlog
!       dfdy(i,j) = i+j
!       print *, 'dfdy(',i,j,')= ', dfdy(i,j)
!    enddo
!   enddo   
 !  print *, 'nlog is =' , nlog
 !  print *, 'nphys is = ', nphys


 !  print *, 'IN pchem_networkDenseJakob'
 !  do i=1,nlog
 !     do j=1,nlog
 !       print *, 'dfdy(',i,j,')= ' , dfdy(i,j)
 !      enddo
 !   enddo


  return
  end subroutine pchem_networkDenseJakob
