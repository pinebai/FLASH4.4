function amr_particle, amr, npart, min=min, max=max, bdbox=bdbox

; returns an array of particles that trace the structure of an
; amr object. min sets the minimum value for which particles will
; be placed, while max sets the threshhold, above which no extra
; particles will be placed. Note that this function places particles
; probabilistically, so the number of particles is not guaranteed to
; be exactly npart, though it should be close for large npart.

; read keywords
if not keyword_set(min) then min=min_amr(amr)
if not keyword_set(max) then max=max_amr(amr)
if keyword_set(bdbox) then begin
    amrbd=amr_add(amr, 0)
    for n=0,amrbd.maxlevel do begin
	for m=0, amrbd.levels[n].nfab-1 do begin
	    fabidxmin=(*amrbd.levels[n].fabptr)[m].idxlo
	    fabidxmax=(*amrbd.levels[n].fabptr)[m].idxhi
	    xlo=(*amrbd.levels[n].fabptr)[m].xlo
	    xhi=(*amrbd.levels[n].fabptr)[m].xhi
	    data=*(*amrbd.levels[n].fabptr)[m].dataptr
	    xpos=xlo[0]+amrbd.gridspacing[0,n]* $
		(findgen(fabidxmax[0]-fabidxmin[0]+1)+0.5)
	    xmask=(xpos ge bdbox[0,0]) and (xpos le bdbox[0,1])
	    loc=where(xmask ne 1)
	    if (loc[0] eq -1) then continue
	    for i=0,fabidxmax[1]-fabidxmin[1] do begin
		for j=0,fabidxmax[2]-fabidxmin[2] do begin
		    (*(*amrbd.levels[n].fabptr)[m].dataptr)[*,i,j] = $
			data[*,i,j]*xmask
		endfor
	    endfor
	    ypos=xlo[1]+amrbd.gridspacing[1,n]* $
		(findgen(fabidxmax[1]-fabidxmin[1]+1)+0.5)
	    ymask=(ypos ge bdbox[1,0]) and (ypos le bdbox[1,1])
	    loc=where(ymask ne 1)
	    if (loc[0] eq -1) then continue
	    for i=0,fabidxmax[0]-fabidxmin[0] do begin
		for j=0,fabidxmax[2]-fabidxmin[2] do begin
		    (*(*amrbd.levels[n].fabptr)[m].dataptr)[i,*,j] = $
			data[i,*,j]*ymask
		endfor
	    endfor
	    zpos=xlo[2]+amrbd.gridspacing[2,n]* $
		(findgen(fabidxmax[2]-fabidxmin[2]+1)+0.5)
	    zmask=(zpos ge bdbox[2,0]) and (zpos le bdbox[2,1])
	    loc=where(zmask ne 1)
	    if (loc[0] eq -1) then continue
	    for i=0,fabidxmax[0]-fabidxmin[0] do begin
		for j=0,fabidxmax[1]-fabidxmin[1] do begin
		    (*(*amrbd.levels[n].fabptr)[m].dataptr)[i,j,*] = $
			data[i,j,*]*zmask
		endfor
	    endfor
	endfor
    endfor
endif else amrbd=amr_add(amr, 0)

; set cells below min to zero
mask=amr_ge(amrbd, min)
amrtmp1=amr_multiply(amr, mask)
amr_free, mask
amr_free, amrbd

;mask=amr_lt(amrbd, min)
;amrdiv=amr_divide(amrbd, min)
;amrdiv1=amr_multiply(amrdiv, mask)
;amrexp=amr_exp(amrdiv1, 1.1)
;amrtmp1=amr_multiply(amrbd, amrexp)
;amr_free, mask
;amr_free, amrdiv
;amr_free, amrdiv1
;amr_free, amrexp

; set cells above max to max
maskle=amr_le(amr, max)
amrtmp2=amr_multiply(amrtmp1, maskle)
maskgt=amr_gt(amr, max)
maskrepl=amr_multiply(maskgt, max)
amrtmp=amr_add(maskrepl, amrtmp2)
amr_free, maskle
amr_free, maskgt
amr_free, maskrepl
amr_free, amrtmp1
amr_free, amrtmp2

; figure out how much weight to assign to each particle
wgttot=amr_sum(amrtmp)
wgtpart=wgttot/npart

; initialize particle list, being conservative about how many particles
; we will need, since extending the list later is slow
partlist=dblarr(2*npart,amr.ndim)
partptr=0L

; loop through grid structure
for n=amrtmp.maxlevel,0,-1 do begin

    ; set the refinement ratio between this level and the next finest level
    if (n ne amrtmp.maxlevel) then refratio = amr.refratio[n]

    ; grab cell volume
    cellvol=1.0d0
    for i=0,amrtmp.ndim-1 do cellvol=cellvol*amrtmp.gridspacing[i,n]

    ; now go through fabs
    for m=0, amrtmp.levels[n].nfab-1 do begin

	; grab the data
	fabidxmin=(*amrtmp.levels[n].fabptr)[m].idxlo
	fabidxmax=(*amrtmp.levels[n].fabptr)[m].idxhi
	xlo=(*amrtmp.levels[n].fabptr)[m].xlo
	xhi=(*amrtmp.levels[n].fabptr)[m].xhi
	data=*(*amrtmp.levels[n].fabptr)[m].dataptr

	; We want to get the list of all fine fabs that overlay
	; this current coarse fab. We will store the result in
	; overlay_list. Don't do this if we're on maxlevel,
	; though.
	if n lt amrtmp.maxlevel then begin
	    overlay_list = lonarr(amr.levels[n+1].nfab+1) - 1
	    overlay_list_ptr = 0
	    for i=0, amr.levels[n+1].nfab-1 do begin

		; get limits of the possibly overlaying fab,
		; coarsened to this level
		overlaymin = (*amr.levels[n+1].fabptr)[i].idxlo $
		    / refratio
		overlaymax = ((*amr.levels[n+1].fabptr)[i].idxhi+1) $
		    / refratio - 1

		; check if this fine fab overlaps our current fab
		if total(fabidxmin gt overlaymax) ne 0 then continue
		if total(fabidxmax lt overlaymin) ne 0 then continue

		; if we're here, this is an overlaying fab, so
		; record its number
		overlay_list[overlay_list_ptr] = i
		overlay_list_ptr = overlay_list_ptr + 1
	    endfor
	endif

	; If there isn't an overlay, we can skip this next part.
	; If there is, we construct a mask to block out cells that
	; are overlayed by finer data.
	if n lt amrtmp.maxlevel then begin

	    ; initialize the mask
	    mask = data * 0

	    ; loop through the overlaying fabs
	    overlay_list_ptr=0
	    while overlay_list[overlay_list_ptr] ne -1 do begin

		; get limits of the possibly overlaying fab,
		; coarsened to this level
		overlaymin = (*amr.levels[n+1].fabptr)$
		   [overlay_list[overlay_list_ptr]].idxlo $
		   / refratio
		overlaymax = ((*amr.levels[n+1].fabptr)$
		   [overlay_list[overlay_list_ptr]].idxhi+1) $
		   / refratio - 1

		; create an object to record the intersection limits
		intersectmin = lonarr(amr.ndim)
		intersectmax = lonarr(amr.ndim)

		; loop through dimensions
		for i=0, amr.ndim-1 do begin

		    ; figure out the limits of the intersection region
		    ; in this dimesion
		    intersectmin[i] = max([fabidxmin[i], overlaymin[i]])
		    intersectmax[i] = min([fabidxmax[i], overlaymax[i]])

		endfor

		; convert to mask / data indices
		maskmin = intersectmin - fabidxmin
		maskmax = intersectmax - fabidxmin

		; add 1 to every mask cell for each dimension where
		; that cell is within the intersection limits
		if amr.ndim eq 1 then begin
		    mask[maskmin[0]:maskmax[0]] = $
			mask[maskmin[0]:maskmax[0]] + 1
		endif
		if amr.ndim eq 2 then begin
		    mask[maskmin[0]:maskmax[0],*] = $
			mask[maskmin[0]:maskmax[0],*] + 1
		    mask[*,maskmin[1]:maskmax[1]] = $
			mask[*,maskmin[1]:maskmax[1]] + 1
		endif
		if amr.ndim eq 3 then begin
		    mask[maskmin[0]:maskmax[0],*,*] = $
			mask[maskmin[0]:maskmax[0],*,*] + 1
		    mask[*,maskmin[1]:maskmax[1],*] = $
			mask[*,maskmin[1]:maskmax[1],*] + 1
		    mask[*,*,maskmin[2]:maskmax[2]] = $
			mask[*,*,maskmin[2]:maskmax[2]] + 1
		endif

		; put a 1 in mask cells that are inside the overlap
		; region in every dimension, a 0 otherwise
		mask = (mask eq amr.ndim)

		; apply the mask to the region
		data = (1 - mask) * data

		; increment the pointer
		overlay_list_ptr = overlay_list_ptr + 1

	    endwhile
	endif

	; divide each cell by weight and multiply by volume
	datawgt=data*cellvol/wgtpart

	; assign integer numbers of particles
	datanpart=floor(datawgt)

	; find remainders
	dataremainder = datawgt-datanpart
	;dataremainder=0

	; get size of data
	sz=size(datawgt)

	; generate random numbers
	if (amrtmp.ndim eq 1) then $
	    datarand = randomu(dummy, sz[1]) $
	else if (amrtmp.ndim eq 2) then $
	    datarand = randomu(dummy, sz[1], sz[2]) $
	else if (amrtmp.ndim eq 3) then $
	    datarand = randomu(dummy, sz[1], sz[2], sz[3])

	; compare random numbers to remainders, and increase particle
	; number appropriately
	addpart = datarand lt dataremainder
	datanpart = datanpart + addpart
	newpart = total(datanpart)
	newptr = 0L

	; find cells that have particles and loop through them
	partloc = where(datanpart ne 0)
	if (partloc[0] eq -1) then continue
	partidx=lonarr(amrtmp.ndim)
	for l=0, n_elements(partloc)-1 do begin

	    ; get x, y, and z index of this cell in the fab
	    fabidx=partloc[l]
	    for i=long(sz[0]-1),0,-1 do begin
		div=1
		for j=1,i do div=div*sz[j]
		partidx[i]=fabidx/div
		fabidx=fabidx mod div
	    endfor

	    ; get physical locations
	    partpos = xlo + partidx*amrtmp.gridspacing[*,n]

	    ; generate random offsets for each particle and add them in
	    offset = randomu(dummy, datanpart[partloc[l]], amrtmp.ndim)
	    ;offset = fltarr(datanpart[partloc[l]], amrtmp.ndim)
	    newpartlist = dblarr(datanpart[partloc[l]], amrtmp.ndim)
	    for i=0, amrtmp.ndim-1 do $
		newpartlist[*,i]=partpos[i]+offset[*,i]*amrtmp.gridspacing[i,n]

	    ; make sure not to overflow particle list
	    if (partptr+datanpart[partloc[l]] ge n_elements(partlist[*,0])) $
		then begin
		parttmp = partlist
		partlist = dblarr(n_elements(parttmp[*,0]), amrtmp.ndim)
		partlist[0:partptr-1,*] = parttmp
	    endif

	    ; add new particles to list
	    partlist[partptr:partptr-1+datanpart[partloc[l]],*] = newpartlist
	    partptr = partptr + datanpart[partloc[l]]

	endfor
    endfor
endfor

amr_free, amrtmp

partlist = partlist[0:partptr-1,*]
return, partlist

end