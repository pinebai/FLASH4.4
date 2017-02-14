pro get_j_amr, px, py, pz, jx, jy, jz, comp=comp, pos=pos

; Procedure to compute the cellwise angular momentum of amr data about
; an arbitrary point. The keyword comp specifies the component of j to
; calculate. Setting comp=0 gives jx, comp=1 gives jy, comp=2 gives
; jz. The component is always returned in its proper slot -- the other
; slots are unaltered if comp is set. If comp is not set, all three
; components are calculated. The position to compute the angular
; momentum about is specified by pos. If it is unset, the origin
; is used.

; Read keywords
if not keyword_set(pos) then pos=fltarr(px.ndim)
if n_elements(comp) eq 0 then comp=-1 else begin
    if ((comp lt 0) or (comp ge amr.ndim)) then begin
	print, 'Error: comp must be in the range 0 to ', amr.ndim-1
	return
    endif
endelse

; Create x, y, and z vectors
if (comp ne 0) then xvec=make_rad_amr(px, comp=0, pos=pos)
if (comp ne 1) then yvec=make_rad_amr(px, comp=1, pos=pos)
if (comp ne 2) then zvec=make_rad_amr(px, comp=2, pos=pos)

; Compute the components
if ((comp eq -1) or (comp eq 0)) then begin
	rypz = amr_multiply(yvec, pz)
	rzpy = amr_multiply(zvec, py)
	jx = amr_subtract(rypz, rzpy)
	amr_free, rypz
	amr_free, rzpy
endif
if ((comp eq -1) or (comp eq 1)) then begin
	rzpx = amr_multiply(zvec, px)
	rxpz = amr_multiply(xvec, pz)
	jy = amr_subtract(rzpx, rxpz)
	amr_free, rzpx
	amr_free, rxpz
endif
if ((comp eq -1) or (comp eq 2)) then begin
	rxpy = amr_multiply(xvec, py)
	rypx = amr_multiply(yvec, px)
	jz = amr_subtract(rxpy, rypx)
	amr_free, rxpy
	amr_free, rypx
endif

; Free memory
if (comp ne 0) then amr_free, xvec
if (comp ne 1) then amr_free, yvec
if (comp ne 2) then amr_free, zvec

return
end