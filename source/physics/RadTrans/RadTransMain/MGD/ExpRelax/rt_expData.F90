!!****if* source/physics/RadTrans/RadTransMain/MGD/ExpRelax/rt_expData
!!
!!  NAME 
!!    rt_expData
!!
!!  SYNOPSIS
!!    use rt_expData
!!
!!  DESCRIPTION 
!!    Stores additional data for ExpRelax MGD
!!
!!***

module rt_expData
  implicit none
  
  ! Maximum number of iterations in the ExpRelax outer loop.
  integer, save :: rt_expRelaxMaxIter

end module rt_expData
