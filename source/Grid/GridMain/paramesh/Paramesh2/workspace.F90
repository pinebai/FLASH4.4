!!****ih* source/Grid/GridMain/paramesh/Paramesh2/workspace
!!
!! NAME
!!
!!    workspace
!!
!!
!! SYNOPSIS
!!
!!   workspace   
!!
!! DESCRIPTION
!!
!!     Defines work space for use with paramesh 2.0
!!  
!!***


module workspace
  use physicaldata
#include "Flash.h"	
  integer len_wblock
  integer ilw,iuw,jlw,juw,klw,kuw
  integer nguard_work,ngw2
  parameter(nguard_work=NGUARD)
  parameter(ngw2=2*nguard_work)
  parameter(len_wblock=(NXB+ngw2)*(NYB+ngw2)*(NZB+ngw2*k3d))
  parameter(ilw=1,iuw=NXB+ngw2)
  parameter(jlw=1,juw=NYB+ngw2)
  parameter(klw=1,kuw=NZB+ngw2*k3d)
  real ::                                                                 &
       &     work(ilw:iuw,jlw:juw,klw:kuw,MAXBLOCKS,1),                    &
       &     recv1(ilw:iuw,jlw:juw,klw:kuw),                              &
       &     send1(ilw:iuw,jlw:juw,klw:kuw),                              &
       &     temp1(ilw:iuw,jlw:juw,klw:kuw)
end module workspace
