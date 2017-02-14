!!****if* source/Simulation/SimulationMain/WindTunnel/user_bnd
!!
!!  NAME
!!    tot_bnd
!!
!!  SYNOPSIS
!!         call user_bnd(solnData,      block_no,neigh  , bndcnd)
!!         call user_bnd(real(:,:,:,:), integer, integer, integer)
!! 
!! DESCRIPTION
!!              Called for a block if the boundary condition is bnd_user.
!!              The routine must fill the appropriate guardcells (as determined
!!              by `neigh', which is 1 for -X, 2 for +X, 3 for -Y, 4 for +Y,
!!              5 for -Z, 6 for +Z) for the given block number (block_no).
!!              For completeness, the boundary condition is passed as well;
!!              at the moment, we expect this to be bnd_user.
!!
!! NOTES
!!              This is intended to be a template.   Please describe here
!!              what any actual functional versions do.
!!
!!***

      subroutine user_bnd( solnData, block_no,neigh  , bndcnd)

!===============================================================================


      use runtime_parameters, only : get_parm_from_context, GLOBAL_PARM_CONTEXT

      use dBase, ONLY: & 
     &     nxb, nyb, nzb, k2d, k3d, nguard, mfaces, ionmax, & 
     &     iLo_gc, jLo_gc, kLo_gc, &
     &     iHi_gc, jHi_gc, kHi_gc, nvar, &
     &     bnd_reflect, bnd_outflow, bnd_user, bnd_hydrostat, & 
     &     dBaseLocalBlockCount, & 
     &     dBaseKeyNumber, dBaseSpecies, & 
     &     dBaseGetDataPtrSingleBlock, & 
     &     dBaseNeighborBlockList,     &
     &     dBaseGetCoords
      implicit none

      real,intent(INOUT), &
           DIMENSION(nvar,iLo_gc:iHi_gc,jLo_gc:jHi_gc,kLo_gc:kHi_gc) &
           :: solnData 

      integer, intent(IN) :: block_no, neigh, bndcnd
      logical, save :: first_call = .true.

      integer :: n, i, j, k

      integer :: ic, kk

      integer, save :: ivelx, ively, ivelz, idens, ipres, itemp, iener, &
           igamc, igame
      integer, save :: inuc_begin

      real, save :: gamma, p_ambient, rho_ambient, smallp, wind_vel
      real :: kine

!===============================================================================

      if (first_call) then

         ivelx = dBaseKeyNumber("velx")
         ively = dBaseKeyNumber("vely")
         ivelz = dBaseKeyNumber("velz")

         idens = dBaseKeyNumber("dens")

         ipres = dBaseKeyNumber("pres")
         itemp = dBaseKeyNumber("temp")
         iener = dBaseKeyNumber("ener")

         igame = dBaseKeyNumber("game")
         igamc = dBaseKeyNumber("gamc")

         inuc_begin = dBaseSpecies(1)

! get the runtime parameters needed for the windtunnel BC
         call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'gamma', gamma)

         call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'p_ambient', p_ambient)
         call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'rho_ambient', rho_ambient)
         call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'smallp', smallp)
         call get_parm_from_context(GLOBAL_PARM_CONTEXT, 'wind_vel', wind_vel)
         


         first_call = .FALSE.
         
      endif


      if (bndcnd .ne. bnd_user) then
         call abort_flash('user_bnd: dont know about non-bnd_user bcs!')
      endif

      

!-------------------------------------------------------------------------------

!               -X user-defined boundary
         
         if (neigh == 1 .and. bndcnd == bnd_user) then

            do k = 1,2*nguard*k3d+nzb
               do j = 1,2*nguard*k2d+nyb
                  do i = 1,nguard
                     
                     solnData(igame,i,j,k) = gamma
                     solnData(igamc,i,j,k) = gamma
                     
                     solnData(idens,i,j,k) = rho_ambient

                     solnData(ivelx,i,j,k) = wind_vel
                     solnData(ively,i,j,k) = 0.
                     solnData(ivelz,i,j,k) = 0.

                     solnData(ipres,i,j,k) = p_ambient

                     kine = 0.5 * (solnData(ivelx,i,j,k)**2)

                     solnData(iener,i,j,k) =  & 
                          solnData(ipres,i,j,k) / & 
                          (solnData(igame,i,j,k)-1.)

                     solnData(iener,i,j,k) =  & 
                          solnData(iener,i,j,k) / & 
                          solnData(idens,i,j,k)

                     solnData(iener,i,j,k) =  & 
                          solnData(iener,i,j,k)  & 
                          + kine

                     solnData(iener,i,j,k) =  & 
                          max(solnData(iener,i,j,k),smallp)

                  end do
               end do
            end do


         endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!               +X user-defined boundary
         
         if (neigh == 2 .and. bndcnd == bnd_user) then

           call abort_flash &
     &       ("tot_bnd:  no +X user boundary condition defined!")

!          do k = 1, 2*nguard*k3d+nzb
!            do j = 1, 2*nguard*k2d+nyb
!              do i = nguard+nxb+1, 2*nguard+nxb
!                solnData(:,i,j,k) = user-defined values
!              end do
!            end do
!          end do

         end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!               -Y user-defined boundary
         
         if (neigh == 3 .and. bndcnd == bnd_user) then

           call abort_flash &
     &       ("tot_bnd:  no -Y user boundary condition defined!")

!          do k = 1, 2*nguard*k3d+nzb
!            do i = 1, 2*nguard+nxb
!              do j = 1, nguard
!                solnData(:,i,j,k) = user-defined values
!              end do
!            end do
!          end do

         end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!               +Y user-defined boundary

         if (neigh == 4 .and. bndcnd == bnd_user) then

           call abort_flash &
     &       ("tot_bnd:  no +Y user boundary condition defined!")

!          do k = 1, 2*nguard*k3d+nzb
!            do i = 1, 2*nguard+nxb
!              do j = 2*nguard+nyb, nguard+nyb+1, -1
!                solnData(:,i,j,k) = user-defined values
!              end do
!            end do
!          end do

         end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!               -Z user-defined boundary

         if (neigh == 5 .and. bndcnd == bnd_user) then

           call abort_flash &
     &       ("tot_bnd:  no -Z user boundary condition defined!")

!          do j = 1, 2*nguard+nyb
!            do i = 1, 2*nguard+nxb
!              do k = 1, nguard
!                solnData(:,i,j,k) = user-defined values
!              end do
!            end do
!          end do

         end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!               +Z user-defined boundary

         if (neigh == 6 .and.  bndcnd == bnd_user) then

           call abort_flash &
     &       ("tot_bnd:  no +Z user boundary condition defined!")

           
!          do j = 1, 2*nguard+nyb
!            do i = 1, 2*nguard+nxb
!              do k = 2*nguard+nzb, nguard+nzb+1, -1
!                solnData(:,i,j,k) = user-defined values
!              end do
!            end do
!          end do

         end if

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      return
      end 

