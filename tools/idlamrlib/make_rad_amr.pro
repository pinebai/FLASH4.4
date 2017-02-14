function make_rad_amr, amr, pos=pos, comp=comp

; Function to make an amr object with the same structure as amr, with
; each cell storing a radial vector. The keyword pos is used to specify
; the location of the origin for the radial vector (default is [0,0,0]).
; The keyword comp specifies to compute only the x (comp=0), y (comp=1),
; or z (comp=2) component of the vector. The default is the compute the
; magnitude of the radial vector.

; Read keywords
if not keyword_set(pos) then pos=fltarr(amr.ndim)
if n_elements(comp) eq 0 then comp=-1 else begin
    if ((comp lt 0) or (comp ge amr.ndim)) then begin
	print, 'Error: comp must be in the range 0 to ', amr.ndim-1
	return, -1
    endif
endelse

; Create the new amr object to hold the result. Null out the pointers
; in it.
amr_new = amr
for n=0, amr_new.maxlevel do begin
	for m=0, amr_new.levels[n].nfab do $
		amr_new.levels[n].fabptr = ptr_new()
endfor

; change component name field in new object
if (comp eq -1) then amr_new.componentname = 'Radius, r=0 at (' $
else if (comp eq 0) then amr_new.componentname = 'x, x=0 at (' $
else if (comp eq 1) then amr_new.componentname = 'y, y=0 at (' $
else if (comp eq 2) then amr_new.componentname = 'z, z=0 at ('
if (amr_new.ndim eq 1) then begin
    amr_new.componentname = amr_new.componentname + $
	strstrim(string(pos[0]),2)+')'
endif else if (amr_new.ndim eq 2) then begin
    amr_new.componentname = amr_new.componentname + $
	strtrim(string(pos[0]),2)+', '+strtrim(string(pos[1]),2)+')'
endif else if (amr_new.ndim eq 3) then begin
    amr_new.componentname = amr_new.componentname + $
	strtrim(string(pos[0]),2)+', '+strtrim(string(pos[1]),2) $
	    +', '+strtrim(string(pos[2]),2)+')'
endif

; go through levels
for n=0, amr_new.maxlevel do begin

    ; Useful shorthand
    dx = amr_new.gridspacing[*,n]

    ; create a new fab pointer for amr_new
    amr_new.levels[n].fabptr = ptr_new(*amr.levels[n].fabptr)

    ; for safety, null out all the data pointers within this fab pointer
    for m=0, amr_new.levels[n].nfab-1 do $
	(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new()

    ; now go through fabs
    for m=0, amr_new.levels[n].nfab-1 do begin

	; Create shorthands
	idxlo = (*amr_new.levels[n].fabptr)[m].idxlo
	idxhi = (*amr_new.levels[n].fabptr)[m].idxhi
	xlo = (*amr_new.levels[n].fabptr)[m].xlo
	xhi = (*amr_new.levels[n].fabptr)[m].xhi

	; create arrays for x, y, and z
	if ((comp eq -1) or (comp eq 0)) then begin
	    xvec = xlo[0] + (dindgen(idxhi[0]-idxlo[0]+1)+0.5)*dx[0] - pos[0]
	    if amr_new.ndim eq 1 then begin
		xarr = xvec
	    endif else if amr_new.ndim eq 2 then begin
		xarr = dblarr(idxhi[0]-idxlo[0]+1, idxhi[1]-idxlo[1]+1)
		for j=idxlo[1],idxhi[1] do xarr[*,j-idxlo[1]]=xvec
	    endif else begin
		xarr = dblarr(idxhi[0]-idxlo[0]+1, idxhi[1]-idxlo[1]+1, $
			      idxhi[2]-idxlo[2]+1)
		for j=idxlo[1],idxhi[1] do for k=idxlo[2],idxhi[2] do $
		    xarr[*,j-idxlo[1],k-idxlo[2]]=xvec
	    endelse
	endif
	if ((comp eq -1) or (comp eq 1)) then begin
	    if amr.ndim ge 2 then $
		yvec = xlo[1] + (dindgen(idxhi[1]-idxlo[1]+1)+0.5)*dx[1] $
		    - pos[1]
	    if amr_new.ndim eq 1 then begin
		yarr = dblarr(idxhi[0]-idxlo[0]+1)
	    endif else if amr_new.ndim eq 2 then begin
		yarr = dblarr(idxhi[0]-idxlo[0]+1, idxhi[1]-idxlo[1]+1)
		for i=idxlo[0],idxhi[0] do yarr[i-idxlo[1],*]=yvec
	    endif else begin
		yarr = dblarr(idxhi[0]-idxlo[0]+1, idxhi[1]-idxlo[1]+1, $
			      idxhi[2]-idxlo[2]+1)
		for i=idxlo[0],idxhi[0] do for k=idxlo[2],idxhi[2] do $
		    yarr[i-idxlo[0],*,k-idxlo[2]]=yvec
	    endelse
	endif
	if ((comp eq -1) or (comp eq 2)) then begin
	    if amr_new.ndim ge 3 then $
		zvec = xlo[2] + (dindgen(idxhi[2]-idxlo[2]+1)+0.5)*dx[2] $
		    - pos[2]
	    if amr_new.ndim eq 1 then begin
		zarr = dblarr(idxhi[0]-idxlo[0]+1)
	    endif else if amr_new.ndim eq 2 then begin
		zarr = dblarr(idxhi[0]-idxlo[0]+1, idxhi[1]-idxlo[1]+1)
	    endif else begin
		zarr = dblarr(idxhi[0]-idxlo[0]+1, idxhi[1]-idxlo[1]+1, $
			      idxhi[2]-idxlo[2]+1)
		for i=idxlo[0],idxhi[0] do for j=idxlo[1],idxhi[1] do $
		    zarr[i-idxlo[0],j-idxlo[1],*]=zvec
	    endelse
	endif

	; store data
	if (comp eq 0) then $
	    (*amr_new.levels[n].fabptr)[m].dataptr = ptr_new(xarr) $
	else if (comp eq 1) then $
	    (*amr_new.levels[n].fabptr)[m].dataptr = ptr_new(yarr) $
	else if (comp eq 2) then $
	    (*amr_new.levels[n].fabptr)[m].dataptr = ptr_new(zarr) $
	else if (comp eq -1) then $
	    (*amr_new.levels[n].fabptr)[m].dataptr = $
		ptr_new(sqrt(xarr^2+yarr^2+zarr^2))

    endfor
endfor

; return the new amr object
return, amr_new
end

