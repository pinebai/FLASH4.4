function amr_vecmag, v1, v2, v3

if n_elements(v2) eq 0 then nvec=1 $
else if n_elements(v3) eq 0 then nvec=2 $
else nvec=3

if nvec eq 1 then begin
	return, amr_abs(v1)
endif else if nvec ge 2 then begin
	v1sqr=amr_multiply(v1,v1)
	v2sqr=amr_multiply(v2,v2)
	vtmp=amr_add(v1sqr,v2sqr)
	vmag=amr_pow(vtmp,0.5)
	amr_free, v1sqr
	amr_free, v2sqr
	amr_free, vtmp
	return, vmag
endif else begin
	v1sqr=amr_multiply(v1,v1)
	v2sqr=amr_multiply(v2,v2)
	v3sqr=amr_multiply(v3,v3)
	vtmp1=amr_add(v1sqr,v2sqr)
	vtmp2=amr_add(vtmp1,v3sqr)
	vmag=amr_pow(vtmp2,0.5)
	amr_free, v1sqr
	amr_free, v2sqr
	amr_free, v3sqr
	amr_free, vtmp1
	amr_free, vtmp2
	return, vmag
endelse

return, -1
end

