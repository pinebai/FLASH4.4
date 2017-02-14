!       PARAMETER ( mvx_m   = max(nxb,nyb,nzb) + 2*nguard )
!       PARAMETER ( mui_m   = nvar )

! The following should work in most cases.

        PARAMETER ( mvx_m   = 128 + 2 * 4 )
        PARAMETER ( mui_m   = 35 )
        PARAMETER ( mvm_m   = mvx_m )
        PARAMETER ( mvxu_m  = mvx_m * mui_m )
        PARAMETER ( mvxmu_m = mvx_m * mvxu_m )
