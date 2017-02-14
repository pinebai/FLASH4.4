!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_dataAprox13
!!
!! NAME
!!  
!!  bn_dataAprox13
!!
!!
!! SYNOPSIS
!! 
!!  use bn_dataAprox13
!!
!! DESCRIPTION
!!
!!  Contains variables required by the Aprox13 nuclear network
!!
!! NOTES
!!
!!   In Flash2, this data was contained in common blocks netc12 and netc13
!!
!!***

module bn_dataAprox13

!..for easy aprox13 rate identification:
!! these were in common block /netc12/
      integer, save ::          ir3a,   irg3a,  ircag,  ir1212, ir1216, iroga,   &   
     &                 iroag,  irnega, irneag, irmgga, irmgag, irsiga,  &
     &                 irmgap, iralpa, iralpg, irsigp, irsiag, irsga,   &
     &                 irsiap, irppa,  irppg,  irsgp,  irsag,  irarga,  &
     &                 irsap,  irclpa, irclpg, irargp, irarag, ircaga,  &
     &                 irarap, irkpa,  irkpg,  ircagp, ircaag, irtiga,  &
     &                 ircaap, irscpa, irscpg, irtigp, irtiag, ircrga,  &
     &                 irtiap, irvpa,  irvpg,  ircrgp, ircrag, irfega,  &
     &                 ircrap, irmnpa, irmnpg, irfegp, irfeag, irniga,  &
     &                 ir1616

!..add these rates for the aprox19 network
!! note that these are also needed in Aprox13!
!! these were in common block /netc13/
      integer, save ::          irpp,   ir33,   ir34,   ircpg,  irnpg,  iropg,   &     
     &                 irnag,  irfeap, ircopa, ircopg, irnigp, irfepg,  &
     &                 ircogp, ir52ng, ir53gn, ir53ng, ir54gn, irheng,  &
     &                 irhegn, irhng,  irdgn,  irdpg,  irhegp, irpen,   &
     &                 ispen,  irnep,  isnep,  irn56ec,isn56ec,ifa,     &
     &                 ifg

end module bn_dataAprox13
