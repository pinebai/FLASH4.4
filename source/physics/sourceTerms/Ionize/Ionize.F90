!!****f* source/physics/sourceTerms/Ionize/Ionize
!!
!! NAME
!!  Ionize
!!
!! SYNOPSIS
!!
!!  call Ionize(integer(IN) :: blockCount,
!!              integer(IN) :: blockList(blockCount),
!!              real(IN)    :: dt,
!!              real(IN)    :: time)
!!
!!
!!
!! DESCRIPTION
!! Apply the ionization operator
!! on the list of blocks provided as input
!!
!! ARGUMENTS
!! blockCount : The number of blocks in the list
!! blockList(:) : The list of blocks on which to apply the stirring operator
!! dt : the current timestep
!! time : the current time
!!
!! NOTES
!!          In the present version there is no evaluation of the
!!          energetic contribution of the ionization and recombination
!!          to the energy equation (heating and cooling the plasma). We
!!          are going to implement this computation soon.
!!
!!          The source term we considered is adequate to solve the
!!          problem for optically thin plasma in the ``coronal''
!!          approximation. This means that we are considering
!!          collisional ionization, auto-ionization, radiative
!!          recombination and dielectronic recombination. The
!!          possibility to change the source term is foreseen to treat
!!          other problems (e.g. the nebular model). We used the
!!          ionization and recombination coefficients computed by
!!          Summers.
!!
!!***

!======================================================================
subroutine Ionize(blockCount,blockList,dt,time)
  
  implicit none
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  return
 end subroutine Ionize


