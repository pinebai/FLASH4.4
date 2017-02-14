; Convenient astrophysical constants. Note that
;   units are all in cgs. 

function constants

const = {				$
; mathematical constants
	 pi 	: 3.141592654,		$
	 G  	: 6.6726e-8, 		$
	 c    	: 2.99792458e10,	$
	 h    	: 6.6260755e-27,	$
;	 hbar 	: h / (2 * pi),		$
	 e 	: 4.8032068e-10,	$
	 me  	: 9.1093898e-28,	$
 	 mp	: 1.672631e-24,		$
 	 mn	: 1.6749286e-24,	$
	 Na 	: 6.0221367e23,		$	
	 kboltz : 1.380658e-16,		$
	 Rgas 	: 8.314511e7,		$
	 sigma  : 5.670e-5,		$
	 a      : 7.56523e-15,		$
	 evtoerg: 1.60217733e-12,	$
	 amu    : 1.6605402e-24,	$

; Astrophysical constants

	Lsun    : 3.847e33,		$
	Rsun 	: 6.96e10,		$
	Msun	: 1.9891e33,		$
	Teffsun : 5780,			$
	pc	: 3.086e18,		$
	au 	: 1.496e13,		$
	Mearth 	: 5.976e27,		$
	Rearth  : 6.37e8,		$
	Mjupiter: 1.899e30,		$
	Rjupiter: 7.140e9,		$

; time constants

	hour	: 3600,			$
	day 	: 86400,		$
	month 	: 2592000,		$
	year	: 31536000}		

return, const

end
