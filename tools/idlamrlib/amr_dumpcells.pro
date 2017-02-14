pro amr_dumpcells, density, xvel, yvel, zvel, maxlevel=maxlevel, $
	fname=fname

; Copyright Mark Krumholz (2001) and Robert Fisher (2002);
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

; Function to determine dumpcells output.

if (keyword_set(fname)) then openw, 1, fname $
else openw, 1, 'dumpcells.out'

; Output dumpcells header information

printf, 1, density.boxmin [0], density.boxmax [0]
printf, 1, density.boxmin [1], density.boxmax [1]
printf, 1, density.boxmin [2], density.boxmax [2]

amridxlo = density.idxlo [*, 0]
amridxhi = density.idxhi [*, 0]

res0 = amridxhi - amridxlo + 1

printf, 1, res0 [0]
printf, 1, res0 [1]
printf, 1, res0 [2]

; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=density.maxlevel

; loop through the levels, starting at the finest
for n=maxlevel, 0, -1 do begin

	; set grid spacing, cell volume, and refinement ratio
	dx = reform(density.gridspacing[*,n])
	cellvol = double(dx[0])
	for i=1, density.ndim-1 do cellvol = cellvol * dx[i]
	if n ne 0 then refratio = density.refratio[n-1] $
	else refratio = 4

	; loop through fabs
	for m=0, density.levels[n].nfab-1 do begin

		; create shortcut for fab limits
		fabidxmin = (*density.levels[n].fabptr)[m].idxlo
		fabidxmax = (*density.levels[n].fabptr)[m].idxhi

		; We want to get the list of all fine fabs that overlay
		; this current coarse fab. We will store the result in
		; overlay_list. Don't do this if we're on maxlevel,
		; though.
		if n lt maxlevel then begin
		   overlay_list = lonarr(density.levels[n+1].nfab+1) - 1
		   overlay_list_ptr = 0
		   for i=0, density.levels[n+1].nfab-1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*density.levels[n+1].fabptr)[i].idxlo $
			   / refratio
			overlaymax = ((*density.levels[n+1].fabptr)[i].idxhi+1) $
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
		densitydata = *(*density.levels[n].fabptr)[m].dataptr
                xveldata    = *(*xvel.levels[n].fabptr)[m].dataptr
                yveldata    = *(*yvel.levels[n].fabptr)[m].dataptr
                zveldata    = *(*zvel.levels[n].fabptr)[m].dataptr

		; If there isn't an overlay, we can skip this next part.
		; If there is, we construct a mask to block out cells that
		; are overlayed by finer data.
		if n lt maxlevel then begin

		   ; initialize the mask
		   mask = densitydata * 0

		   ; loop through the overlaying fabs
		   overlay_list_ptr=0
		   while overlay_list[overlay_list_ptr] ne -1 do begin

			; get limits of the possibly overlaying fab,
			; coarsened to this level
			overlaymin = (*density.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxlo $
			   / refratio
			overlaymax = ((*density.levels[n+1].fabptr)$
			   [overlay_list[overlay_list_ptr]].idxhi+1) $
			   / refratio - 1

			; create an object to record the intersection limits
			intersectmin = lonarr(density.ndim)
			intersectmax = lonarr(density.ndim)

			; loop through dimensions
			for i=0, density.ndim-1 do begin

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
			if density.ndim eq 1 then begin
			   mask[maskmin[0]:maskmax[0]] = $
				mask[maskmin[0]:maskmax[0]] + 1
			endif
			if density.ndim eq 2 then begin
			   mask[maskmin[0]:maskmax[0],*] = $
				mask[maskmin[0]:maskmax[0],*] + 1
			   mask[*,maskmin[1]:maskmax[1]] = $
				mask[*,maskmin[1]:maskmax[1]] + 1
			endif
			if density.ndim eq 3 then begin
			   mask[maskmin[0]:maskmax[0],*,*] = $
				mask[maskmin[0]:maskmax[0],*,*] + 1
			   mask[*,maskmin[1]:maskmax[1],*] = $
				mask[*,maskmin[1]:maskmax[1],*] + 1
			   mask[*,*,maskmin[2]:maskmax[2]] = $
				mask[*,*,maskmin[2]:maskmax[2]] + 1
			endif

			; put a 1 in mask cells that are inside the overlap
			; region in every dimension, a 0 otherwise
			mask = (mask eq density.ndim)

			; apply the mask to the region
			densitydata = (1 - mask) * densitydata

			; increment the pointer
			overlay_list_ptr = overlay_list_ptr + 1

		   endwhile
		endif

                ; iterate over nonmasked data

                nonmasked = where (densitydata, count)

		if (nonmasked[0] eq -1) then continue

                sizenonmasked = size (nonmasked)
                length = sizenonmasked [1]

                ; output nonmasked data

                nx = fabidxmax [0] - fabidxmin [0] + 1
                ny = fabidxmax [1] - fabidxmin [1] + 1
                nz = fabidxmax [2] - fabidxmin [2] + 1

                for i = 0, length - 1 do begin

                  ; get indices of this unmasked cell

                  kdx = nonmasked [i] / (nx * ny)
                  jdx = (nonmasked [i] - kdx * (nx * ny) ) / nx
                  idx = nonmasked [i] - jdx * nx - kdx * (nx * ny)

                  kdx = kdx + fabidxmin [2]
                  jdx = jdx + fabidxmin [1]
                  idx = idx + fabidxmin [0]

                  ; convert to equivalent index on finest level
                  ; max level 7, fixed refinement ratio

                  kdx = kdx * refratio^(7 - n)
                  jdx = jdx * refratio^(7 - n)
                  idx = idx * refratio^(7 - n)

                  printf, 1, idx, jdx, kdx, $
                        densitydata [nonmasked [i]], $
                        xveldata    [nonmasked [i]], $
                        yveldata    [nonmasked [i]], $
                        zveldata    [nonmasked [i]], $
                        dx [0], format = '(3I, 5E)'
                endfor

	endfor
endfor

close, 1

return
end

