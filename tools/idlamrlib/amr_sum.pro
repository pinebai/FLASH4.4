function amr_sum, amr, maxlevel=maxlevel, xrange=xrange, yrange=yrange, zrange=zrange

; Copyright Mark Krumholz (2001);
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

; This function does a volume-weighted sum of an amr object. If
; maxlevel is set, the sum only uses data on levels up to maxlevel.


; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel
if not keyword_set(xrange) then xrange=[amr.boxmin[0], amr.boxmax[0]]
if amr.ndim gt 1 then $
	if not keyword_set(yrange) then yrange=[amr.boxmin[1], amr.boxmax[1]]
if amr.ndim gt 2 then $
	if not keyword_set(zrange) then zrange=[amr.boxmin[2], amr.boxmax[2]]

; initialize the data sum
amrsum = 0.0d

; loop through the levels, starting at the finest
for n=maxlevel, 0, -1 do begin

	; set grid spacing, cell volume, and refinement ratio
	dx = reform(amr.gridspacing[*,n])
	cellvol = double(dx[0])
	for i=1, amr.ndim-1 do cellvol = cellvol * dx[i]
	if n ne 0 then refratio = amr.refratio[n-1]

	; figure out the index limits on this level corresponding to physical
	; limits given
        xidxlim = lonarr(2)
        xidxlim[0] = floor((xrange[0] - amr.boxmin[0]) / dx[0])
        xidxlim[1] = ceil((xrange[1] - amr.boxmin[0]) / dx[0]) - 1
	if amr.ndim gt 1 then begin
		yidxlim = lonarr(2)
	        yidxlim[0] = floor((yrange[0] - amr.boxmin[1]) / dx[1])
        	yidxlim[1] = ceil((yrange[1] - amr.boxmin[1]) / dx[1]) - 1
	endif
	if amr.ndim gt 2 then begin
		zidxlim = lonarr(2)
	        zidxlim[0] = floor((zrange[0] - amr.boxmin[2]) / dx[1])
        	zidxlim[1] = ceil((zrange[1] - amr.boxmin[2]) / dx[1]) - 1
	endif

	; loop through fabs
	for m=0, amr.levels[n].nfab-1 do begin

		; create shortcut for fab limits
		fabidxmin = (*amr.levels[n].fabptr)[m].idxlo
		fabidxmax = (*amr.levels[n].fabptr)[m].idxhi

		; see how this fab overlaps with the desired range
                overlapmin = lonarr(amr.ndim)
                overlapmax = lonarr(amr.ndim)
                overlapmin[0] = (fabidxmin[0] gt xidxlim[0]) * fabidxmin[0] + $
                        (fabidxmin[0] le xidxlim[0]) * xidxlim[0]
                overlapmax[0] = (fabidxmax[0] le xidxlim[1]) * fabidxmax[0] + $
                        (fabidxmax[0] gt xidxlim[1]) * xidxlim[1]
		if amr.ndim gt 1 then begin
                   overlapmin[1] = (fabidxmin[1] gt yidxlim[0]) * fabidxmin[1] + $
                        (fabidxmin[1] le yidxlim[0]) * yidxlim[0]
                   overlapmax[1] = (fabidxmax[1] le yidxlim[1]) * fabidxmax[1] + $
                        (fabidxmax[1] gt yidxlim[1]) * yidxlim[1]
		endif
		if amr.ndim gt 2 then begin
                   overlapmin[2] = (fabidxmin[2] gt zidxlim[0]) * fabidxmin[2] + $
                        (fabidxmin[2] le zidxlim[0]) * zidxlim[0]
                   overlapmax[2] = (fabidxmax[2] le zidxlim[1]) * fabidxmax[2] + $
                        (fabidxmax[2] gt zidxlim[1]) * zidxlim[1]
		endif

                ; if there is no overlap between this fab and the
                ; image box, then move to the next fab
                if (overlapmin[0] gt overlapmax[0]) then continue
                if amr.ndim gt 1 then $
                   if (overlapmin[1] gt overlapmax[1]) then continue
                if amr.ndim gt 2 then $
                   if (overlapmin[2] gt overlapmax[2]) then continue

                ; there is an overlap, so set the limits on the portion of
                ; the fab we wish to extract
                idxlo = overlapmin - fabidxmin
                idxhi = overlapmax - fabidxmin

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

		; extract the portion of the data in the desired region
		if amr.ndim eq 1 then data = data[idxlo[0]:idxhi[0]] $
		else if amr.ndim eq 2 then data = data[idxlo[0]:idxhi[0], $
						       idxlo[1]:idxhi[1]] $
		else if amr.ndim eq 3 then data = data[idxlo[0]:idxhi[0], $
						       idxlo[1]:idxhi[1], $
						       idxlo[2]:idxhi[2]]

		; add the data * cell volume
		amrsum = amrsum + total(data) * cellvol

	endfor
endfor

; return
return, amrsum

end








