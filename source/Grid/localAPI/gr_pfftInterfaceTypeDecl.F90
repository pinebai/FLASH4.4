!!****ih* source/Grid/localAPI/gr_pfftInterfaceTypeDecl
!!
!! NAME
!!
!!  gr_pfftInterfaceTypeDecl
!!
!! SYNOPSIS
!!
!!  use gr_pfftInterfaceTypeDecl
!!
!! DESCRIPTION
!!
!!  Contains derived data type declarations used by gr_pfftInterface.
!!  We need to place the derived data types in a module in order that 
!!  different subroutines have access to the same type declaration.
!!
!!***

module gr_pfftInterfaceTypeDecl
  implicit none
  type PossibleGrid_t
     integer :: jProcs
     integer :: kProcs
  end type PossibleGrid_t
end module gr_pfftInterfaceTypeDecl
