!!****if* source/physics/Eos/EosMain/Helmholtz/Ye/eos_initHelmholtz
!!
!! NAME
!!
!!  eos_initHelmholtz
!!
!! SYNOPSIS
!!
!!  call eos_initHelmholtz()
!!
!! DESCRIPTION
!!
!!  Initialize the Helmholtz EOS.  The table data is read in on processor 0
!!  and broadcast to all other processors at the start of FLASH execution.
!!  This routine first checks for a binary copy of the table (helm_table.bdat), 
!!  and then for the ASCII version (helm_table.dat).  If the binary table 
!!  file was not found, it is created by this routine for subsequent use.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!   eos_coulombMult[Real, default 1.0] -- Coulomb correction multiplier. Set
!!               to zero to ignore the Coloumb correction.
!!   eos_tolerance[Real, 1.0e-8] -- Convergence tolerance for the Newton-Rhapson
!!               iterations
!!   eos_maxNewton[Integer, 50] -- Maximum number of Newton-Raphson iterations
!!   eos_forceConstantInput     -- This switch forces the Eos implementation
!!                                 to never modify EINT or PRES in MODE_DENS_EI
!!                                 and MODE_DENS_PRES. If this is .false. (the
!!                                 default), calls to Eos may slightly modify
!!                                 these input variables in order to preserve
!!                                 thermodynamic equilibrium.
!!
!!  NOTES
!!
!!  Helmholtz law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Flash.h
!!
!!***



subroutine eos_initHelmholtz()

  use Eos_data, ONLY : eos_type, eos_meshMe, eos_smallt
  use eos_helmData 
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  ! vector_eos.fh computes the vector length from nxb, nyb, nzb, so 
  ! this information must be provided

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
  include 'Flash_mpi.h'
  integer:: unitEos =2
  integer :: i, j
  real :: tstp, dstp
  integer :: istat, ierr


  !get the runtime parameters

  call RuntimeParameters_get('smallt', eos_smallt)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
  call RuntimeParameters_get('eos_coulombMult', eos_coulombMult)
  call RuntimeParameters_get('eos_coulombAbort', eos_coulombAbort)
  call RuntimeParameters_get("eos_forceConstantInput",eos_forceConstantInput)

#ifndef USE_EOS_YE
  write(*,*)"WARNING!! When using Eos/Helmholtz/Ye, generally USE_EOS_YE should be defined!")
  call Logfile_stampMessage( &
       "WARNING! When using Eos/Helmholtz/Ye, generally USE_EOS_YE should be defined!")
#endif


  if (eos_meshMe==MASTER_PE) then
     open(unit=unitEos,file='helm_table.bdat',status='old',iostat=istat)
     close(unit=unitEos)
     print*,'about to open file'
     istat=1
     if (istat.ne.0) then
        write(*,*) '[Eos_init] Cannot open helm_table.bdat!'
        write(*,*) '[Eos_init] Trying old helm_table.dat!'
        
        open(unit=unitEos,file='helm_table.dat',status='old',iostat=istat)

        if (istat .ne. 0) then
           write(*,*) '[Eos_init] ERROR: opening helm_table.dat!'
           call Driver_abortFlash("[Eos_init] ERROR: opening helm_table.dat")
        endif

        !..read the helmholtz free energy table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_f(i,j),eos_fd(i,j),eos_ft(i,j),&
                   eos_fdd(i,j),eos_ftt(i,j), & 
                   eos_fdt(i,j), eos_fddt(i,j),eos_fdtt(i,j),eos_fddtt(i,j)
           enddo
        enddo

        !..read the pressure derivative with density table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_dpdf(i,j),eos_dpdfd(i,j),&
                   eos_dpdft(i,j),eos_dpdfdt(i,j)
           enddo
        enddo

        !..read the electron chemical potential table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_ef(i,j),eos_efd(i,j),&
                   eos_eft(i,j),eos_efdt(i,j)
           enddo
        enddo

        !..read the number density table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_xf(i,j),eos_xfd(i,j),&
                   eos_xft(i,j),eos_xfdt(i,j)
           enddo
        enddo

        !..close up the data file and write a message
        close(unitEos)

        !..dump binary version of table for later use
        istat = EOSIMAX*EOSJMAX
        call eos_writeHfet(istat, & 
             eos_f,eos_fd,eos_ft,eos_fdd,&
             eos_ftt,eos_fdt,eos_fddt,eos_fdtt,eos_fddtt, & 
             eos_dpdf,eos_dpdfd,eos_dpdft,eos_dpdfdt, & 
             eos_ef,eos_efd,eos_eft,eos_efdt, & 
             eos_xf,eos_xfd,eos_xft,eos_xfdt)

        !..read binary version of table
     else
        istat = EOSIMAX*EOSJMAX
        call eos_readHfet(istat, & 
             eos_f,eos_fd,eos_ft,eos_fdd,&
             eos_ftt,eos_fdt,eos_fddt,eos_fdtt,eos_fddtt, & 
             eos_dpdf,eos_dpdfd,eos_dpdft,eos_dpdfdt, & 
             eos_ef,eos_efd,eos_eft,eos_efdt, & 
             eos_xf,eos_xfd,eos_xft,eos_xfdt)
     endif
  endif

  !..broadcast to rest of processors
  istat = EOSIMAX*EOSJMAX
  call MPI_BCAST(eos_f,      istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_fd,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_ft,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_fdd,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_ftt,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_fdt,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_fddt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_fdtt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_fddtt,  istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_dpdf,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_dpdfd,  istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_dpdft,  istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_dpdfdt, istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_ef,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_efd,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_eft,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_efdt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_xf,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_xfd,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_xft,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_xfdt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)


  eos_tlo   = 4.0e0
  tstp  = (11.0e0 - eos_tlo)/float(EOSJMAX-1)
  eos_tstpi = 1.0e0/tstp
  eos_dlo   = -10.0e0
  dstp  = (11.0e0 - eos_dlo)/float(EOSIMAX-1)
  eos_dstpi = 1.0e0/dstp
  do j=1,EOSJMAX
     eos_t(j) = 10.0e0**(eos_tlo + (j-1)*tstp)
     do i=1,EOSIMAX
        eos_d(i) = 10.0e0**(eos_dlo + (i-1)*dstp)
     enddo
  enddo


  !..store the temperature and density differences and their inverses 
  do j=1,EOSJMAX-1
     eos_dt(j)   = eos_t(j+1) - eos_t(j)
     eos_dtSqr(j)  = eos_dt(j)*eos_dt(j)
     eos_dtInv(j)  = 1.0e0/eos_dt(j)
     eos_dtSqrInv(j) = 1.0e0/eos_dtSqr(j)
  enddo
  do i=1,EOSIMAX-1
     eos_dd(i)   = eos_d(i+1) - eos_d(i)
     eos_ddSqr(i)  = eos_dd(i)*eos_dd(i)
     eos_ddInv(i)  = 1.0e0/eos_dd(i)
     eos_ddSqrInv(i) = 1.0e0/eos_ddSqr(i)
  enddo
  eos_type=EOS_HLM

  return
end subroutine eos_initHelmholtz
