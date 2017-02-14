!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_extractRaysData
!!
!! NAME
!!
!!  ed_extractRaysData
!!
!! SYNOPSIS
!!
!!  use ed_extractRaysData
!!
!! DESCRIPTION
!!
!!  This module allows the external user to extract all data from all the rays.
!!  The character keyword necessary to identify the type of data wanted must
!!  be the same as the name of the data variable under which the ray data
!!  is stored. The main function is overloaded with specific functions for each
!!  data type. If any of the ray number ID, the ray tag ID or the character keyword
!!  for the data does not match which what is currently stored for the processor,
!!  the routine returns the value of .false. in the dataFound keyword.
!!
!! ARGUMENTS
!!
!!  rayID      : the ID number of the ray                       ! optional
!!  rayTag     : the global Tag number of the ray               ! optional
!!  entryField : a character string identifying the data
!!  dataValue  : the returned value of the data requested
!!  dataFound  : logical keyword indicating if the data was found
!!
!! NOTES
!!
!!  No overloading is currently needed, as all data entries of the rays are real.
!!  The module character of this operation is for future changes, if the ray field
!!  will contain mixed types of data.
!!
!!***

Module ed_extractRaysData

implicit none

interface ed_extractRayData
   module procedure ed_extractRayDataReal
end interface

contains
!
!
!     ...The real version.
!
!
subroutine ed_extractRayDataReal (rayID,                  &
                                  rayTag,                 &
                                  entryField,             &
                                               dataValue, &
                                               dataFound  )

  use Driver_interface,       ONLY : Driver_abortFlash

  use EnergyDeposition_data,  ONLY : ed_maxRayCount, &
                                     ed_rays

  implicit none

#include "EnergyDeposition.h"

  integer, optional, intent (in)  :: rayID
  integer, optional, intent (in)  :: rayTag
  integer,           intent (in)  :: entryField
  real,              intent (out) :: dataValue
  logical,           intent (out) :: dataFound

  logical :: useRayID
  logical :: useRayTag

  integer :: n
  integer :: ray
  integer :: tag

  dataFound = .true.                ! default value
  useRayID  = present (rayID)
  useRayTag = present (rayTag)

  if (useRayID .and. useRayTag) then
      dataValue = 0.0
      dataFound = .false.
      return
  end if

  if (useRayTag) then
      ray = 0
      do n = 1,ed_maxRayCount
         tag = int (ed_rays (RAY_TAGS,n))
         if (tag == rayTag) then
             ray = n
             exit
         end if
      end do
  end if

  if (useRayID) then
      ray = rayID
  end if

  if ((ray < 1) .or. (ray > ed_maxRayCount) ) then
       dataValue = 0.0
       dataFound = .false.
       return
  end if

  select case (entryField)

  case (RAY_POSX)
         dataValue = ed_rays (RAY_POSX,ray)
  case (RAY_POSY)
         dataValue = ed_rays (RAY_POSY,ray)
  case (RAY_POSZ)
         dataValue = ed_rays (RAY_POSZ,ray)
  case (RAY_VELX)
         dataValue = ed_rays (RAY_VELX,ray)
  case (RAY_VELY)
         dataValue = ed_rays (RAY_VELY,ray)
  case (RAY_VELZ)
         dataValue = ed_rays (RAY_VELZ,ray)
  case (RAY_POWR)
         dataValue = ed_rays (RAY_POWR,ray)
  case (RAY_DENC)
         dataValue = ed_rays (RAY_DENC,ray)
  case (RAY_BLCK)
         dataValue = ed_rays (RAY_BLCK,ray)
  case (RAY_PROC)
         dataValue = ed_rays (RAY_PROC,ray)
  case (RAY_TAGS)
         dataValue = ed_rays (RAY_TAGS,ray)
  case default
         dataFound = .false.
  end select

  return
end subroutine ed_extractRayDataReal

end Module ed_extractRaysData
