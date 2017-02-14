function get_ppv_amr, rho, v, dir, xrange=xrange, yrange=yrange,$
	verbose=verbose, maxlevel=maxlevel, nvelbin=nvelbin,	$
	velbinlim=velbinlim, velres=velres, rhocut=rhocut,	$
	velbinctr=velbinctr, smearsize=smearsize,		$
	smeardist=smeardist, csound=csound, zrange=zrange

; Copyright Mark Krumholz (2002)
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

; This routine observes an amr object along the direction specified by
; dir and returns ab observation in the form of a ppv (position-position
; -velocity) cube. The user can specify the number of velocity bins (via
; the nvelbin keyword), can explicitly specify a series of velocity bins
; (by setting velbinlim equal to an array containing the bin limits), and
; can specify the minimum density rho for which the observation will
; see anything (corresponding to an excitation threshhold). All
; observations are assumed to be in the optically thin limit. If velbinlim
; is not set, it returns the velocity bin limits used. The keyword
; velbinctr returns the centers of the velocity bins.

; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=rho.maxlevel
if not keyword_set(csound) then csound=0.0

; figure out which direction we're observing along
if rho.ndim lt 3 then begin
	print, 'Error: rho object must have 3 dimensions.'
	return, -1
endif
if dir eq 0 then begin
	dim1=1
	dim2=2
	dim3=0
endif else if dir eq 1 then begin
	dim1=0
	dim2=2
	dim3=1
endif else if dir eq 2 then begin
	dim1=0
	dim2=1
	dim3=2
endif else begin
	print, 'Error: direction must be 0, 1, or 2'
	return, -1
endelse

; set up the velocity bins
if n_elements(velbinlim) lt 2 then begin
	velmax = max_amr(v, velmin)
	velmax = velmax + 0.001*(velmax-velmin) + 3*csound
	velmin = velmin - 0.001*(velmax-velmin) - 3*csound
	if not keyword_set(velres) then begin
		if not keyword_set(nvelbin) then nvelbin = 20
		velbinlim = velmin + dindgen(nvelbin+1)*(velmax-velmin)/nvelbin
	endif else begin
		nvelbin = ceil((velmax-velmin)/velres)
		velbinlim = velmin + dindgen(nvelbin+1)*velres
	endelse
endif else nvelbin=n_elements(velbinlim)-1
velbinctr=0.5*(velbinlim[0:nvelbin-1]+velbinlim[1:nvelbin])

; set range we're going to extract
if not keyword_set(xrange) then xlim=[rho.boxmin[dim1], rho.boxmax[dim1]] $
else xlim=xrange
if not keyword_set(yrange) then ylim=[rho.boxmin[dim2], rho.boxmax[dim2]] $
else ylim=yrange
if not keyword_set(zrange) then zlim=[rho.boxmin[dim3], rho.boxmax[dim3]] $
else zlim=zrange

; if asked to do so, snap range to grid on finest level
if keyword_set(xrange) and (n_elements(snaptogrid) ne 0) then begin
	xlim[0] = floor( (xlim[0] - rho.boxmin[dim1]) / $
			 rho.gridspacing[dim1,maxlevel] ) * $
		  rho.gridspacing[dim1, maxlevel] + $
		  rho.boxmin[dim1]
	xlim[1] = ceil( (xlim[1] - rho.boxmin[dim1]) / $
			 rho.gridspacing[dim1,maxlevel] ) * $
		  rho.gridspacing[dim1, maxlevel] + $
		  rho.boxmin[dim1]
	if keyword_set(verbose) then print, 'New x range = ', xlim
endif
if keyword_set(yrange) and (n_elements(snaptogrid) ne 0) then begin
	ylim[0] = floor( (ylim[0] - rho.boxmin[dim2]) / $
			 rho.gridspacing[dim2,maxlevel] ) * $
		  rho.gridspacing[dim2, maxlevel] + $
		  rho.boxmin[dim2]
	ylim[1] = ceil( (ylim[1] - rho.boxmin[dim2]) / $
			 rho.gridspacing[dim2,maxlevel] ) * $
		  rho.gridspacing[dim2, maxlevel] + $
		  rho.boxmin[dim2]
	if keyword_set(verbose) then print, 'New y range: ', ylim
endif
if keyword_set(zrange) and (n_elements(snaptogrid) ne 0) then begin
	zlim[0] = floor( (zlim[0] - rho.boxmin[dim3]) / $
			 rho.gridspacing[dim3,maxlevel] ) * $
		  rho.gridspacing[dim3, maxlevel] + $
		  rho.boxmin[dim3]
	zlim[1] = ceil( (zlim[1] - rho.boxmin[dim3]) / $
			 rho.gridspacing[dim3,maxlevel] ) * $
		  rho.gridspacing[dim3, maxlevel] + $
		  rho.boxmin[dim3]
	if keyword_set(verbose) then print, 'New z range: ', zlim
endif

; set column image size equal to number of cells on maxlevel
; covered by the requested range
collo = [floor((xlim[0]-rho.boxmin[dim1]) / $
		rho.gridspacing[dim1,maxlevel]), $
	 floor((ylim[0]-rho.boxmin[dim2]) / $
		rho.gridspacing[dim2,maxlevel])]
colhi = [ceil((xlim[1]-rho.boxmin[dim1]) / $
		rho.gridspacing[dim1,maxlevel]), $
	 ceil((ylim[1]-rho.boxmin[dim2]) / $
		rho.gridspacing[dim2,maxlevel])] - 1
ppv=fltarr(colhi[0]-collo[0]+1, colhi[1]-collo[1]+1, nvelbin)

; go through amr structure, filling appropriate values into
; column image box
boxmin = [rho.boxmin[dim1], rho.boxmin[dim2], rho.boxmin[dim3]]
boxmax = [rho.boxmax[dim1], rho.boxmax[dim2], rho.boxmax[dim3]]
for n=maxlevel, 0, -1 do begin

	; set the refinement ratio between this level and the finest level
	if n ne 0 then refratio = rho.refratio[n-1]

	; create a shorthand for gridspacing
	dx = [rho.gridspacing[dim1,n], rho.gridspacing[dim2,n], $
	      rho.gridspacing[dim3,n]]

	; figure out the indices on this level corresponding the the
	; physical limits given
	xidxlim = lonarr(2)
	yidxlim = lonarr(2)
	zidxlim = lonarr(2)
	xidxlim[0] = floor((xlim[0] - boxmin[0]) / dx[0])
	xidxlim[1] = ceil((xlim[1] - boxmin[0]) / dx[0]) - 1
	yidxlim[0] = floor((ylim[0] - boxmin[1]) / dx[1])
	yidxlim[1] = ceil((ylim[1] - boxmin[1]) / dx[1]) - 1
	zidxlim[0] = floor((zlim[0] - boxmin[2]) / dx[2])
	zidxlim[1] = ceil((zlim[1] - boxmin[2]) / dx[2]) - 1

	; loop through fabs
	for m=0, rho.levels[n].nfab-1 do begin

		if keyword_set(verbose) then print, 'Level ', $
			strtrim(string(n),2), ', fab ', strtrim(string(m),2)

		; set up some shorthands
		fabidxmin = (*rho.levels[n].fabptr)[m].idxlo
		fabidxmax = (*rho.levels[n].fabptr)[m].idxhi

		; find the overlap between our target index range and the
		; index range stored in this fab
		overlapmin = lonarr(3)
		overlapmax = lonarr(3)
		overlapmin[0] = (fabidxmin[dim1] gt xidxlim[0]) * $
			fabidxmin[dim1] + $
			(fabidxmin[dim1] le xidxlim[0]) * xidxlim[0]
		overlapmin[1] = (fabidxmin[dim2] gt yidxlim[0]) * $
			fabidxmin[dim2] + $
			(fabidxmin[dim2] le yidxlim[0]) * yidxlim[0]
		overlapmin[2] = (fabidxmin[dim3] gt zidxlim[0]) * $
			fabidxmin[dim3] + $
			(fabidxmin[dim3] le zidxlim[0]) * zidxlim[0]
		overlapmax[0] = (fabidxmax[dim1] le xidxlim[1]) * $
			fabidxmax[dim1] + $
			(fabidxmax[dim1] gt xidxlim[1]) * xidxlim[1]
		overlapmax[1] = (fabidxmax[dim2] le yidxlim[1]) * $
			fabidxmax[dim2] + $
			(fabidxmax[dim2] gt yidxlim[1]) * yidxlim[1]
		overlapmax[2] = (fabidxmax[dim3] le zidxlim[1]) * $
			fabidxmax[dim3] + $
			(fabidxmax[dim3] gt zidxlim[1]) * zidxlim[1]

		; if there is no overlap between this fab and the
		; image box, then move to the next fab
		if (overlapmin[0] gt overlapmax[0]) or $
		   (overlapmin[1] gt overlapmax[1]) or $
		   (overlapmin[2] gt overlapmax[2]) then continue

		; We want to get the list of all fine fabs that overlay
		; this current coarse fab. We will store the result in
		; overlay_list. Don't do this if we're on maxlevel,
		; though.
		if n lt maxlevel then begin
		   overlay_list = lonarr(rho.levels[n+1].nfab+1) - 1
		   overlay_list_ptr = 0
		   for i=0, rho.levels[n+1].nfab-1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*rho.levels[n+1].fabptr)[i].idxlo $
			   / refratio
			overlaymax = ((*rho.levels[n+1].fabptr)[i].idxhi+1) $
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

		; grab the data for this region
		rhofab = *(*rho.levels[n].fabptr)[m].dataptr
		vfab = *(*v.levels[n].fabptr)[m].dataptr

		; If there isn't an overlay, we can skip this next part.
		; If there is, we construct a mask to block out cells that
		; are overlayed by finer data.
		if n lt maxlevel then begin

		   ; initialize the mask
		   mask = rhofab * 0

		   ; loop through the overlaying fabs
		   overlay_list_ptr=0
		   while overlay_list[overlay_list_ptr] ne -1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*rho.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxlo $
			   / refratio
			overlaymax = ((*rho.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxhi+1) $
			   / refratio - 1

			; create an object to record the intersection limits
			intersectmin = lonarr(rho.ndim)
			intersectmax = lonarr(rho.ndim)

			; loop through dimensions
			for i=0, rho.ndim-1 do begin

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
			if rho.ndim eq 1 then begin
			   mask[maskmin[0]:maskmax[0]] = $
				mask[maskmin[0]:maskmax[0]] + 1
			endif
			if rho.ndim eq 2 then begin
			   mask[maskmin[0]:maskmax[0],*] = $
				mask[maskmin[0]:maskmax[0],*] + 1
			   mask[*,maskmin[1]:maskmax[1]] = $
				mask[*,maskmin[1]:maskmax[1]] + 1
			endif
			if rho.ndim eq 3 then begin
			   mask[maskmin[0]:maskmax[0],*,*] = $
				mask[maskmin[0]:maskmax[0],*,*] + 1
			   mask[*,maskmin[1]:maskmax[1],*] = $
				mask[*,maskmin[1]:maskmax[1],*] + 1
			   mask[*,*,maskmin[2]:maskmax[2]] = $
				mask[*,*,maskmin[2]:maskmax[2]] + 1
			endif

			; put a 1 in mask cells that are inside the overlap
			; region in every dimension, a 0 otherwise
			mask = (mask eq rho.ndim)

			; apply the mask to the region
			rhofab = (1 - mask) * rhofab

			; increment the pointer
			overlay_list_ptr = overlay_list_ptr + 1

		   endwhile
		endif

		; extract the portion of the data that intersects the
		; requested image box
		xidxlist = lindgen(overlapmax[0]-overlapmin[0]+1) $
			+ overlapmin[0]
		yidxlist = lindgen(overlapmax[1]-overlapmin[1]+1) $
			+ overlapmin[1]
		zidxlist = lindgen(overlapmax[2]-overlapmin[2]+1) $
			+ overlapmin[2]
		datasz = size(rhofab)
		rhofab = reform(rhofab, 1, datasz[1], datasz[2], datasz[3], $
			      /overwrite)
		vfab = reform(vfab, 1, datasz[1], datasz[2], datasz[3], $
			      /overwrite)
		if dir eq 0 then begin
			rhofab = reform(rhofab[0, zidxlist-fabidxmin[dim3], $
				xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(zidxlist), $
				n_elements(xidxlist), $
				n_elements(yidxlist))
			vfab = reform(vfab[0, zidxlist-fabidxmin[dim3], $
				xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(zidxlist), $
				n_elements(xidxlist), $
				n_elements(yidxlist))
		endif else if dir eq 1 then begin
			rhofab = reform(rhofab[0, xidxlist-fabidxmin[dim1], $
				zidxlist-fabidxmin[dim3], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(xidxlist), $
				n_elements(zidxlist), $
				n_elements(yidxlist))
			vfab = reform(vfab[0, xidxlist-fabidxmin[dim1], $
				zidxlist-fabidxmin[dim3], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(xidxlist), $
				n_elements(zidxlist), $
				n_elements(yidxlist))
		endif else begin
			rhofab = reform(rhofab[0, xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2], $
				zidxlist-fabidxmin[dim3]], $
				n_elements(xidxlist), $
				n_elements(yidxlist), $
				n_elements(zidxlist))
			vfab = reform(vfab[0, xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2], $
				zidxlist-fabidxmin[dim3]], $
				n_elements(xidxlist), $
				n_elements(yidxlist), $
				n_elements(zidxlist))
		endelse

		; set densities below the cutoff to zero
		if keyword_set(rhocut) then $
			rhofab = rhofab * (rhofab gt rhocut)

		; refine the overlap index range to maxlevel
		overlapminref = overlapmin
		overlapmaxref = overlapmax + 1
		for l=maxlevel-1,n,-1 do begin
		   overlapminref = overlapminref * rho.refratio[l]
		   overlapmaxref = overlapmaxref * rho.refratio[l]
		endfor
		overlapmaxref = overlapmaxref - 1

		imgidxlo = long(overlapminref[0:1]-collo)
		imgidxhi = long(overlapmaxref[0:1]-collo)

		; crop the image box index range to fit what we
		; have available
		imgidxlo = (imgidxlo gt 0) * imgidxlo
		imgidxhi = (imgidxhi le (colhi-collo)) * imgidxhi + $
			   (imgidxhi gt (colhi-collo)) * $
				(colhi-collo)

		; loop through velocity bins
		for i=0, nvelbin-1 do begin

			; create a mask for this bin
			vsize=size(vfab)
			if not keyword_set(csound) then begin
				; just bin assuming no maxwellian distrib.
				vmask = (vfab ge velbinlim[i]) and $
					(vfab lt velbinlim[i+1])
			endif else begin
				; bin by maxwellian distribution
				vmat1 = (velbinlim[i]-vfab) / $
					(sqrt(2.)*csound)
				vmat2 = (velbinlim[i+1]-vfab) / $
					(sqrt(2.)*csound)
				vmask = 0.5 * (errorf(vmat2)-errorf(vmat1))
			endelse
			vmask=reform(vmask, vsize[1], vsize[2], vsize[3])

			; now sum in the appropriate direction
			rhofabtmp=reform(rhofab*vmask, vsize[1], vsize[2], $
					 vsize[3])
			column = reform(rho.gridspacing[dim3,n] * $
					total(rhofabtmp, dim3 + 1, $
					/double), $
				        n_elements(xidxlist), $
					n_elements(yidxlist))

			; add data to ppv array
			ppv[imgidxlo[0]:imgidxhi[0], $
		    	    imgidxlo[1]:imgidxhi[1], i] = $
		 	  ppv[imgidxlo[0]:imgidxhi[0], $
		    	      imgidxlo[1]:imgidxhi[1], i] + $
				congrid(column, imgidxhi[0]-imgidxlo[0]+1, $
					imgidxhi[1]-imgidxlo[1]+1)

		endfor

	endfor
endfor

; smear the data with a beam if requested
if keyword_set(smearsize) then begin
	for n=0,nvelbin-1 do begin
		ppv[*,*,n] = beamsmear(ppv[*,*,n], $
				       rho.gridspacing[0,maxlevel], $
				       smeardist, smearsize)
	endfor
endif

; return the ppv data
return, ppv

end
