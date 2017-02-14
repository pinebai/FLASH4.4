function get_column_amr, amr, plane, xrange=xrange, yrange=yrange, $
	verbose=verbose, maxlevel=maxlevel

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

; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel

; figure out which direction we're summing along
if amr.ndim lt 3 then begin
	print, 'Error: amr object must have 3 dimensions.'
	return, -1
endif
if plane eq 0 then begin
	dim1=1
	dim2=2
	dim3=0
endif else if plane eq 1 then begin
	dim1=0
	dim2=2
	dim3=1
endif else if plane eq 2 then begin
	dim1=0
	dim2=1
	dim3=2
endif else begin
	print, 'Error: plane must be 0, 1, or 2'
	return, -1
endelse

; set range of plot range
if not keyword_set(xrange) then xlim=[amr.boxmin[dim1], amr.boxmax[dim1]] $
else xlim=xrange
if not keyword_set(yrange) then ylim=[amr.boxmin[dim2], amr.boxmax[dim2]] $
else ylim=yrange

; if asked to do so, snap range to grid on finest level
if keyword_set(xrange) and (n_elements(snaptogrid) ne 0) then begin
	xlim[0] = floor( (xlim[0] - amr.boxmin[dim1]) / $
			 amr.gridspacing[dim1,maxlevel] ) * $
		  amr.gridspacing[dim1, maxlevel] + $
		  amr.boxmin[dim1]
	xlim[1] = ceil( (xlim[1] - amr.boxmin[dim1]) / $
			 amr.gridspacing[dim1,maxlevel] ) * $
		  amr.gridspacing[dim1, maxlevel] + $
		  amr.boxmin[dim1]
	if keyword_set(verbose) then print, 'New x range = ', xlim
endif
if keyword_set(yrange) and (n_elements(snaptogrid) ne 0) then begin
	ylim[0] = floor( (ylim[0] - amr.boxmin[dim2]) / $
			 amr.gridspacing[dim2,maxlevel] ) * $
		  amr.gridspacing[dim2, maxlevel] + $
		  amr.boxmin[dim2]
	ylim[1] = ceil( (ylim[1] - amr.boxmin[dim2]) / $
			 amr.gridspacing[dim2,maxlevel] ) * $
		  amr.gridspacing[dim2, maxlevel] + $
		  amr.boxmin[dim2]
	if keyword_set(verbose) then print, 'New y range: ', ylim
endif

; set column image size equal to number of cells on maxlevel
; covered by the requested range
imglo = [floor((xlim[0]-amr.boxmin[dim1]) / $
		amr.gridspacing[dim1,maxlevel]), $
	 floor((ylim[0]-amr.boxmin[dim2]) / $
		amr.gridspacing[dim2,maxlevel])]
imghi = [ceil((xlim[1]-amr.boxmin[dim1]) / $
		amr.gridspacing[dim1,maxlevel]), $
	 ceil((ylim[1]-amr.boxmin[dim2]) / $
		amr.gridspacing[dim2,maxlevel])] - 1
img=fltarr(imghi[0]-imglo[0]+1, imghi[1]-imglo[1]+1)

; go through amr structure, filling appropriate values into
; column image box
boxmin = [amr.boxmin[dim1], amr.boxmin[dim2]]
boxmax = [amr.boxmax[dim1], amr.boxmax[dim2]]
for n=maxlevel, 0, -1 do begin

	; set the refinement ratio between this level and the finest level
	if n ne 0 then refratio = amr.refratio[n-1]

	; create a shorthand for gridspacing
	dx = [amr.gridspacing[dim1,n], amr.gridspacing[dim2,n]]

	; figure out the indices on this level corresponding the the
	; physical limits given
	xidxlim = lonarr(2)
	yidxlim = lonarr(2)
	xidxlim[0] = floor((xlim[0] - boxmin[0]) / dx[0])
	xidxlim[1] = ceil((xlim[1] - boxmin[0]) / dx[0]) - 1
	yidxlim[0] = floor((ylim[0] - boxmin[1]) / dx[1])
	yidxlim[1] = ceil((ylim[1] - boxmin[1]) / dx[1]) - 1

	; loop through fabs
	for m=0, amr.levels[n].nfab-1 do begin

		if keyword_set(verbose) then print, 'Level ', $
			strtrim(string(n),2), ', fab ', strtrim(string(m),2)

		; set up some shorthands
		fabidxmin = (*amr.levels[n].fabptr)[m].idxlo
		fabidxmax = (*amr.levels[n].fabptr)[m].idxhi

		; find the overlap between our target index range and the
		; index range stored in this fab
		overlapmin = lonarr(2)
		overlapmax = lonarr(2)
		overlapmin[0] = (fabidxmin[dim1] gt xidxlim[0]) * $
			fabidxmin[dim1] + $
			(fabidxmin[dim1] le xidxlim[0]) * xidxlim[0]
		overlapmin[1] = (fabidxmin[dim2] gt yidxlim[0]) * $
			fabidxmin[dim2] + $
			(fabidxmin[dim2] le yidxlim[0]) * yidxlim[0]
		overlapmax[0] = (fabidxmax[dim1] le xidxlim[1]) * $
			fabidxmax[dim1] + $
			(fabidxmax[dim1] gt xidxlim[1]) * xidxlim[1]
		overlapmax[1] = (fabidxmax[dim2] le yidxlim[1]) * $
			fabidxmax[dim2] + $
			(fabidxmax[dim2] gt yidxlim[1]) * yidxlim[1]

		; if there is no overlap between this fab and the
		; image box, then move to the next fab
		if (overlapmin[0] gt overlapmax[0]) or $
		   (overlapmin[1] gt overlapmax[1]) then continue

		; We want to get the list of all fine fabs that overlay
		; this current coarse fab. We will store the result in
		; overlay_list. Don't do this if we're on maxlevel,
		; though.
		if n lt maxlevel then begin
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

		; grab the data for this region
		data = *(*amr.levels[n].fabptr)[m].dataptr

		; If there isn't an overlay, we can skip this next part.
		; If there is, we construct a mask to block out cells that
		; are overlayed by finer data.
		if n lt maxlevel then begin

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

		; extract the portion of the data that intersects the
		; requested image box
		xidxlist = lindgen(overlapmax[0]-overlapmin[0]+1) $
			+ overlapmin[0]
		yidxlist = lindgen(overlapmax[1]-overlapmin[1]+1) $
			+ overlapmin[1]
		zidxlist = lindgen(fabidxmax[dim3]-fabidxmin[dim3]+1) $
			+ fabidxmin[dim3]
		datasz = size(data)
		data = reform(data, 1, datasz[1], datasz[2], datasz[3], $
			      /overwrite)
		if plane eq 0 then begin
			data =  reform(data[0, zidxlist-fabidxmin[dim3], $
				xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(zidxlist), $
				n_elements(xidxlist), $
				n_elements(yidxlist))
		endif else if plane eq 1 then begin
			data =  reform(data[0, xidxlist-fabidxmin[dim1], $
				zidxlist-fabidxmin[dim3], $
				yidxlist-fabidxmin[dim2]], $
				n_elements(xidxlist), $
				n_elements(zidxlist), $
				n_elements(yidxlist))
		endif else begin
			data =  reform(data[0, xidxlist-fabidxmin[dim1], $
				yidxlist-fabidxmin[dim2], $
				zidxlist-fabidxmin[dim3]], $
				n_elements(xidxlist), $
				n_elements(yidxlist), $
				n_elements(zidxlist))
		endelse

		; now sum in the appropriate direction
		column = reform(amr.gridspacing[dim3,n] * $
				total(data, dim3 + 1, /double), $
			        n_elements(xidxlist), n_elements(yidxlist))

		; refine the overlap index range to maxlevel
		overlapminref = overlapmin
		overlapmaxref = overlapmax + 1
		for l=maxlevel-1,n,-1 do begin
		   overlapminref = overlapminref * amr.refratio[l]
		   overlapmaxref = overlapmaxref * amr.refratio[l]
		endfor
		overlapmaxref = overlapmaxref - 1

		imgidxlo = long(overlapminref-imglo)
		imgidxhi = long(overlapmaxref-imglo)

		; crop the image box index range to fit what we
		; have available
		imgidxlo = (imgidxlo gt 0) * imgidxlo
		imgidxhi = (imgidxhi le (imghi-imglo)) * imgidxhi + $
			   (imgidxhi gt (imghi-imglo)) * $
				(imghi-imglo)

		; add data to image array
		img[imgidxlo[0]:imgidxhi[0], $
		    imgidxlo[1]:imgidxhi[1]] = $
		   img[imgidxlo[0]:imgidxhi[0], $
		       imgidxlo[1]:imgidxhi[1]] + $
			congrid(column, imgidxhi[0]-imgidxlo[0]+1, $
				imgidxhi[1]-imgidxlo[1]+1)

	endfor
endfor

; return the column image box
return, img

end
