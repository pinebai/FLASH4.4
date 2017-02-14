!!****if* source/Grid/GridMain/paramesh/paramesh4/amr_restrict_fc_ins
!!
!! NAME
!!
!!  amr_restrict_fc_ins
!!
!! SYNOPSIS
!!
!!  call amr_restrict_fc_ins(Real, Intent(in)  :: recv,
!!                           Real, Intent(inout)  :: temp,
!!                           Integer, Intent(in)  :: icoord,
!!                           Integer, Intent(in)  :: order,
!!                           Integer, Intent(in)  :: ivar)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   recv : 
!!
!!   temp : 
!!
!!   icoord : 
!!
!!   order : 
!!
!!   ivar : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****if* source/amr_restrict_fc_ins.F90
!! NAME
!!   amr_restrict_fc_ins
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_fc_ins(recv,temp,icoord,order,ivar)
!!   Call amr_restrict_fc_ins(real array,real array,integer,integer,integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: recv(:,:,:,:)  Face-centered data to restrict             
!!   Real,    Intent(inout) :: temp(:,:,:,:)  Restricted face-centered data              
!!   Integer, Intent(in)    :: icoord         selects which edge to operate on           
!!                                            iccord = 1, selects x-face                 
!!                                            iccord = 2, selects y-face                 
!!                                            iccord = 3, selects z-face                 
!!   Integer, Intent(in)    :: order          order of Lagrange polynomial to use        
!!   Integer, Intent(in)    :: ivar           which variable in unk_e_x(y,z) to          
!!                                            operate on.                                
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! DESCRIPTION
!!
!!   Stub for a routine that performs interpolation for the restriction operation on
!!   on face-centered data stored in 'facevar?', adapted for use by incompressible 
!!   Navier-Stokes solver.
!!
!! AUTHORS
!!
!!   Named _ins        Klaus Weide      Sept  2016
!!***

Subroutine amr_restrict_fc_ins(recv,temp,icoord,order,ivar)

    use Driver_interface, ONLY: Driver_abortFlash

    Implicit None

!-----Input/Output arguments
    Real,    Intent(in)    :: recv(:,:,:,:)
    Real,    Intent(inout) :: temp(:,:,:,:)
    Integer, Intent(in)    :: icoord, order, ivar


    call Driver_abortFlash(&
           "Attempting to use amr_restrict_fc_ins, not comfigured in!")


End Subroutine amr_restrict_fc_ins
