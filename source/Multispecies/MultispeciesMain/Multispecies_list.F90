!!****if* source/Multispecies/MultispeciesMain/Multispecies_list
!!
!! NAME
!!  Multispecies_list
!!
!! SYNOPSIS
!!
!!  Multispecies_list(integer(in) :: fileUnit)            
!!
!! DESCRIPTION
!!
!!  Writes the Multispecies to a already open file
!! 
!!
!! ARGUMENTS
!!
!!     fileUnit - file number to write
!!
!! NOTES
!!
!!***            

subroutine Multispecies_list(fileUnit)
  
  use Multispecies_data !, ONLY   : ms_array
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(in)                :: fileUnit
  integer                            :: i, n
  character(len=4)            :: speciesName
  

  !  List the multispecies

  write(fileUnit,900)
  do i=1, NSPECIES
     !! Since species are defined in the Config file, they are actually NUMBERS
     call Simulation_mapIntToStr(ms_Array(i)%name,speciesName,MAPBLOCK_UNK)

     write(fileUnit,901)speciesName, &
          ms_Array(i)%name, &
          ms_Array(i)%numTotal, &
          ms_Array(i)%numPositive, &
          ms_Array(i)%numNeutral, &
          ms_Array(i)%numNegative, &
          ms_Array(i)%bindingEnergy, &
          ms_Array(i)%adiabaticIndex, &
          ms_Array(i)%eosType
  enddo

  ! Print an additional block of more esoteric fields currently mostly used
  ! with multiTemp setups, only if a multiTemp Eos is used:
#if N_EOS_TEMP > 0
  write(fileUnit,920)
  do i=1, NSPECIES
     !! Since species are defined in the Config file, they are actually NUMBERS
     call Simulation_mapIntToStr(ms_Array(i)%name,speciesName,MAPBLOCK_UNK)

     write(fileUnit,921)speciesName, &
          ms_Array(i)%name, &
          ms_Array(i)%opacityLowTemp, &
          ms_Array(i)%zMin, &
          ms_Array(i)%eosSubtype, &
          ms_Array(i)%eosZfreeTableFile, &
          ms_Array(i)%eosEnerTableFile, &
          ms_Array(i)%eosPresTableFile
  enddo
#endif

  ! Print an additional block of even more esoteric fields if useful.
  if (ANY(ms_Array(:)%eosGroupName(1:6) .NE. "-none-")) then
  
     write(fileUnit,930)
     do i=1, NSPECIES
        !! Since species are defined in the Config file, they are actually NUMBERS
        call Simulation_mapIntToStr(ms_Array(i)%name,speciesName,MAPBLOCK_UNK)

        write(fileUnit,931)speciesName, &
          ms_Array(i)%name, &
          ms_Array(i)%eosGroupName
     enddo
  end if

  write(fileUnit, '(/a)') "Species Constituents"

  do i=1, NSPECIES
     !! Since species are defined in the Config file, they are actually NUMBERS
     call Simulation_mapIntToStr(ms_Array(i)%name,speciesName,MAPBLOCK_UNK)

     write(fileUnit, '(a,i6,a)') speciesName, ms_Array(i)%numElems, " constituents"

     if(ms_Array(i)%numElems > 0) then
        
        write(fileUnit, '(a3,a3,a11,a11)') "#","Z","A","Frac"
        
        do n = 1, ms_Array(i)%numElems
           write(fileUnit, '(i3,i3,e11.3,e11.3)') &
                n, &
                ms_Array(i)%zElems(n), &
                ms_Array(i)%aElems(n), &
                ms_Array(i)%fractions(n)
        end do
        
        write(fileUnit, *)

     end if
  enddo
     



  return           
  !------------------------------------------------------------------------           
900 format("Initially defined values of species:",/,                        &
         &         "Name",T10,"Index",T25,                                           &
         &         "Total",T33,"Positive",T43,"Neutral",T53,"Negative",T63,         &
         &         "bind Ener",T73,"Gamma",T81,"eosType")
  
901 format(A4,T10,I5,6X,6(ES9.2,1x), &
         I5)

920 format(/,                        &
         &         "Name",T10,"Index",T17,                                      &
         &         "OpacityLowTemp",T33,"Zmin",T41,"Subtype",T49, &
         &         "ZFreeTableFile",T82,"EnerTableFile",T115,"PresTableFile")
  
921 format(A4,T10,I5,5X, &
         2(1x,ES9.2), &
         1x,I7, &
         3(1X,A32))

930 format(/,                        &
         &         "Name",T10,"Index",T21,                                      &
         &         "GroupName")
  
931 format(A4,T10,I5,T21, A128)


end subroutine Multispecies_list
