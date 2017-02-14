!!****if* source/physics/Eos/EosNuclear/kernel/readtable
!!
!! NAME
!!
!!  readtable
!!
!! SYNOPSIS
!!
!!  call readtable(:: eos_filename)
!!
!! DESCRIPTION
!!
!! Reads the Eos table  
!!
!! ARGUMENTS
!!
!!   eos_filename : the file name of the table
!!
!!
!!
!!***

subroutine readtable(eos_filename)
! This routine reads the table and initializes
! all variables in the module. 

  use eosmodule
  use hdf5 

  implicit none

  character(*) eos_filename

  character(len=100) message

! HDF5 vars
  integer(HID_T) file_id,dset_id,dspace_id
  integer(HSIZE_T) dims1(1), dims3(3)
  integer error,rank,accerr
  integer i,j,k, ind

  real*8 amu_cgs_andi
  real*8 buffer1,buffer2,buffer3,buffer4

  real, allocatable :: tempArray(:,:,:)
  real, allocatable :: logtemp0(:)
  accerr=0

!  write(*,*) "Reading Ott EOS Table"

  call h5open_f(error)

  call h5fopen_f (trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, file_id, error)

!  write(6,*) trim(adjustl(eos_filename))

! read scalars
  dims1(1)=1
  call h5dopen_f(file_id, "pointsrho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
  call h5dclose_f(dset_id,error)

  if(error.ne.0) then
     stop "Could not read EOS table file"
  endif

  dims1(1)=1
  call h5dopen_f(file_id, "pointstemp", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
  call h5dclose_f(dset_id,error)

  if(error.ne.0) then
     stop "Could not read EOS table file"
  endif

  dims1(1)=1
  call h5dopen_f(file_id, "pointsye", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
  call h5dclose_f(dset_id,error)

  if(error.ne.0) then
     stop "Could not read EOS table file"
  endif

!  write(message,"(a25,i5,i5,i5)") "We have nrho ntemp nye: ", nrho,ntemp,nye
!  write(*,*) message

  allocate(logrho(nrho))
  dims1(1)=nrho
  call h5dopen_f(file_id, "logrho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  allocate(logtemp0(ntemp))
!  write(*,*) 'have allocated'
  dims1(1)=ntemp
  call h5dopen_f(file_id, "logtemp", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp0, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
!  write(*,*) 'have read'
  allocate(ye(nye))
  dims1(1)=nye
  call h5dopen_f(file_id, "ye", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

!!$  do i = 1, ntemp
!!$     if (logtemp0(i) <= log10(eos_table_tmax)) ind = i
!!$  enddo
!  write(*,*) 'old, new ntemp', ntemp, ind
  ind=ntemp

  allocate(tempArray(nrho,ntemp,nye))
  allocate(logtemp(ind))
  logtemp = logtemp0(1:ind)

!  allocate(alltables(nrho,ntemp,nye,nvars))
  allocate(alltables(nrho,ind,nye,nvars))

  ! index variable mapping:
  !  1 -> logpress
  !  2 -> logenergy
  !  3 -> entropy
  !  4 -> munu
  !  5 -> cs2
  !  6 -> dedT
  !  7 -> dpdrhoe
  !  8 -> dpderho
  !  9 -> muhat
  ! 10 -> mu_e
  ! 11 -> mu_p
  ! 12 -> mu_n
  ! 13 -> xa
  ! 14 -> xh
  ! 15 -> xn
  ! 16 -> xp
  ! 17 -> abar
  ! 18 -> zbar


  dims3(1)=nrho
  dims3(2)=ntemp
  dims3(3)=nye
  call h5dopen_f(file_id, "logpress", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,1) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "logenergy", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,2) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "entropy", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,3) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "munu", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,4) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "cs2", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,5) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "dedt", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,6) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "dpdrhoe", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,7) = tempArray(:,1:ind,:)
  call h5dopen_f(file_id, "dpderho", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,8) = tempArray(:,1:ind,:)

! chemical potentials
  call h5dopen_f(file_id, "muhat", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,9) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "mu_e", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,10) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "mu_p", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,11) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "mu_n", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,12) = tempArray(:,1:ind,:)

! compositions
  call h5dopen_f(file_id, "Xa", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,13) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "Xh", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,14) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "Xn", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,15) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "Xp", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,16) = tempArray(:,1:ind,:)


! average nucleus
  call h5dopen_f(file_id, "Abar", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,17) = tempArray(:,1:ind,:)

  call h5dopen_f(file_id, "Zbar", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,18) = tempArray(:,1:ind,:)

! Gamma
  call h5dopen_f(file_id, "gamma", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempArray(:,:,:), dims3, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error
  alltables(:,:,:,19) = tempArray(:,1:ind,:)


  call h5dopen_f(file_id, "energy_shift", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, energy_shift, dims1, error)
  call h5dclose_f(dset_id,error)
  accerr=accerr+error

  if(accerr.ne.0) then
    stop "Problem reading EOS table file"
  endif


  call h5fclose_f (file_id,error)

  call h5close_f (error)

  ntemp = ind

  ! set min-max values:

  eos_rhomin = 10.0d0**logrho(1)
  eos_rhomax = 10.0d0**logrho(nrho)

  eos_yemin = ye(1)
  eos_yemax = ye(nye)

  eos_tempmin = 10.0d0**logtemp(1)
  eos_tempmax = 10.0d0**logtemp(ntemp)

  deallocate(tempArray)

!  write(6,*) "Done reading eos tables"


end subroutine readtable


