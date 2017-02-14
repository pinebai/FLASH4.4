!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_finalize
!!
!! NAME
!!
!!  Diffuse_finalize
!!
!! SYNOPSIS
!!
!!  Diffuse_finalize()
!!
!! DESCRIPTION
!!
!!  Finalizes local data for Unit Diffuse 
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!***

subroutine Diffuse_finalize() 
  
  implicit none  
  
  call diff_saFinalize ()  
  
  return
end subroutine Diffuse_finalize
