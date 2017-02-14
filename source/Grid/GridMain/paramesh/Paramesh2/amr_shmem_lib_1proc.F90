!
! This file contains a replacement library for all the shmem routines
! called by the package, which will enable the user to compile and
! run their code on a single processor machine which does not have
! the true shmem library or mpi. This may be useful for debugging
! purposes.
!


! $RCSfile: amr_shmem_lib_1proc.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $



        subroutine shmem_real_get(target,source,len,pe)

        implicit none
        real    target(len),source(len)
        integer len,pe

        target(:) = source(:)

        return
        end



        subroutine shmem_integer_get(target,source,len,pe)

        implicit none
        integer target(len),source(len)
        integer len,pe

        target(:) = source(:)

        return
        end



        subroutine shmem_real_put(target,source,len,pe)

        implicit none
        real    target(len),source(len)
        integer len,pe

        target(:) = source(:)

        return
        end



        subroutine shmem_integer_put(target,source,len,pe)

        implicit none
        integer target(len),source(len)
        integer len,pe

        target(:) = source(:)

        return
        end



        subroutine shmem_udcflush()

        implicit none

        return
        end



        subroutine barrier()

        implicit none

        return
        end



        subroutine shmem_barrier_all()

        implicit none

        return
        end



        subroutine shmem_real8_min_to_all(target,source,nred,pestart, & 
     &          pestride,pesize,pwrk,ipsync)

        implicit none
        real    target,source,pwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end



        subroutine shmem_real8_max_to_all(target,source,nred,pestart, & 
     &          pestride,pesize,pwrk,ipsync)

        implicit none
        real    target,source,pwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end



        subroutine shmem_real8_sum_to_all(target,source,nred,pestart, & 
     &          pestride,pesize,pwrk,ipsync)

        implicit none
        real    target,source,pwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end



        subroutine shmem_int8_min_to_all(target,source,nred,pestart, & 
     &          pestride,pesize,ipwrk,ipsync)

        implicit none
        integer         target,source,ipwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end



        subroutine shmem_int8_max_to_all(target,source,nred,pestart, & 
     &          pestride,pesize,ipwrk,ipsync)

        implicit none
        integer         target,source,ipwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end



        subroutine shmem_int8_sum_to_all(target,source,nred,pestart, & 
     &          pestride,pesize,ipwrk,ipsync)

        implicit none
        integer         target,source,ipwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end



        integer function shmem_my_pe()

        implicit none
        shmem_my_pe = 0
        return
        end



        integer function shmem_n_pes()

        implicit none
        shmem_n_pes = 1
        return
        end
