function max_amr, amr, min, maxidx=maxidx, minidx=minidx

; Copyright Mark Krumholz (2001)
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


; this routine finds the maximum value in an amr object. It works just
; like the standard idl routine max, except that the maxidx and minidx
; return an (2*ndim + 2) element array as follows:
;    ( level, fabnum, (xyz) index within fab, (xyz) index within level )
; the index within the fab refers to the index within that particular
; fab array, from 0 to number of cells in fab - 1, and the index within
; the level refers to the indices within the level, starting at 0 on
; that level and going to number of cells on that level - 1.

; get machine-specific minimum
mathprec=machar(/double)
globalmax = -mathprec.xmax
maxidx=lonarr(2*amr.ndim+2)
if n_params() eq 2 then begin
	globalmin = mathprec.xmax
	minidx = lonarr(2*amr.ndim+2)
endif

; go through amr structure
for n=0, amr.maxlevel do begin
	for m=0, amr.levels[n].nfab-1 do begin

		; get max
		localmax = max(*(*amr.levels[n].fabptr)[m].dataptr, idx)
		if localmax gt globalmax then begin
			globalmax = localmax
			maxidxsave = idx
			maxidx[0] = n
			maxidx[1] = m
		endif

		; get min if request
		if n_params() eq 2 then begin
			localmin = $
			   min(*(*amr.levels[n].fabptr)[m].dataptr, idx)
			if localmin lt globalmin then begin
				globalmin = localmin
				minidxsave = idx
				minidx[0] = n
				minidx[1] = m
			endif
		endif

	endfor
endfor

; get max index
idxlo = (*amr.levels[maxidx[0]].fabptr)[maxidx[1]].idxlo
ncell = (*amr.levels[maxidx[0]].fabptr)[maxidx[1]].idxhi - $
	(*amr.levels[maxidx[0]].fabptr)[maxidx[1]].idxlo + 1
if amr.ndim eq 1 then begin
	maxidx[2] = maxidxsave
	maxidx[3] = maxidx[2] + idxlo[0]
endif else if amr.ndim eq 2 then begin
	maxidx[2] = maxidxsave mod ncell[0]
	maxidx[3] = maxidxsave / ncell[0]
	maxidx[4] = maxidx[2] + idxlo[0]
	maxidx[5] = maxidx[3] + idxlo[1]
endif else if amr.ndim eq 3 then begin
	maxidx[2] = maxidxsave mod ncell[0]
	maxidx[3] = ((maxidxsave - maxidx[2])/ncell[0]) mod ncell[1]
	maxidx[4] = (maxidxsave - maxidx[2] - ncell[0]*maxidx[3]) / $
			(ncell[0] * ncell[1])
	maxidx[5] = maxidx[2] + idxlo[0]
	maxidx[6] = maxidx[3] + idxlo[1]
	maxidx[7] = maxidx[4] + idxlo[2]
endif else begin
	print, 'Error: ndim must be 1, 2, or 3'
	return, -1
endelse

; get min index
if n_params() eq 2 then begin
	min = globalmin
	idxlo = (*amr.levels[minidx[0]].fabptr)[minidx[1]].idxlo
	ncell = (*amr.levels[minidx[0]].fabptr)[minidx[1]].idxhi - $
		(*amr.levels[minidx[0]].fabptr)[minidx[1]].idxlo + 1
	if amr.ndim eq 1 then begin
		minidx[2] = minidxsave
		minidx[3] = minidx[2] + idxlo[0]
	endif else if amr.ndim eq 2 then begin
		minidx[2] = minidxsave mod ncell[0]
		minidx[3] = minidxsave / ncell[0]
		minidx[4] = minidx[2] + idxlo[0]
		minidx[5] = minidx[3] + idxlo[1]
	endif else if amr.ndim eq 3 then begin
		minidx[2] = minidxsave mod ncell[0]
		minidx[3] = ((minidxsave - minidx[2])/ncell[0]) mod ncell[1]
		minidx[4] = (minidxsave - minidx[2] - ncell[0]*minidx[3]) / $
			(ncell[0] * ncell[1])
		minidx[5] = minidx[2] + idxlo[0]
		minidx[6] = minidx[3] + idxlo[1]
		minidx[7] = minidx[4] + idxlo[2]
	endif else begin
		print, 'Error: ndim must be 1, 2, or 3'
		return, -1
	endelse
endif

return, globalmax

end



