function amr_flatten, amr, level, idxmin=idxmin, idxmax=idxmax, $
	sample=sample, minlevel=minlevel

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


; this routine takes an amr object and interpolates / coarsens it
; onto a single level of refinement, specified by level. The keywords
; idxmin and idxmax specify the minimum and maximum indices (on the
; specified level) of the region to be flattened (default is the
; entire amr volume). The keyword sample specifies that refinement /
; coarsening should be done by sampling. The default is bilinear
; interpolation.

; error check
if (level lt 0) or (level gt amr.maxlevel) then begin
	print, 'Error: level must be between 0 and ', $
		strtrim(string(amr.maxlevel),2)
	return, -1
endif

; read keywords
if not keyword_set(idxmin) then idxmin=amr.idxlo[*,level] $
else begin
	if total(idxmin lt amr.idxlo[*,level]) ne 0 then begin
		print, 'Error: minimum indices on level ', $
			strtrim(string(level),2), ' are: ', $
			reform(amr.idxlo[*,level])
		return, -1
	endif
endelse
if not keyword_set(idxmax) then idxmax=amr.idxhi[*,level] $
else begin
	if total(idxmin gt amr.idxhi[*,level]) ne 0 then begin
		print, 'Error: maximum indices on level ', $
			strtrim(string(level),2), ' are: ', $
			reform(amr.idxhi[*,level])
		return, -1
	endif
endelse
if not keyword_set(minlevel) then minlevel=0

; create object to store flattened data
if amr.ndim eq 1 then data=dblarr(idxmax[0]-idxmin[0]+1) $
else if amr.ndim eq 2 then data=dblarr(idxmax[0]-idxmin[0]+1, $
				       idxmax[1]-idxmin[1]+1) $
else if amr.ndim eq 3 then data=dblarr(idxmax[0]-idxmin[0]+1, $
				       idxmax[1]-idxmin[1]+1, $
				       idxmax[2]-idxmin[2]+1) $
else begin
	print, 'Error: ndim must be 1, 2, or 3'
	return, -1
endelse
datasize=idxmax-idxmin+1

; now proceed through levels and fabs
for n=minlevel, level do begin

	for m=0, amr.levels[n].nfab-1 do begin

		; figure out refinement ratio between this level and target
		; level
		refratio = 1
		for l=n, level-1 do refratio = refratio * amr.refratio[l]

                ; refine fab indices to target level
                fablo = (*amr.levels[n].fabptr)[m].idxlo * refratio
                fabhi = ((*amr.levels[n].fabptr)[m].idxhi+1) * $
                        refratio - 1
                fabsize = fabhi - fablo + 1

                ; determine if this fab is within the domain we're
                ; looking for. If not, skip it.
                if (total(fablo gt idxmax) ne 0) or $
		   (total(fabhi lt idxmin) ne 0) then continue

		; refine fab to target level
		fab = (*amr.levels[n].fabptr)[m].dataptr
		if keyword_set(sample) then begin
			if amr.ndim eq 1 then $
				fabref = rebin(*fab, fabsize[0], /sample) $
			else if amr.ndim eq 2 then $
				fabref = rebin(*fab, fabsize[0], fabsize[1], $
					       /sample) $
			else if amr.ndim eq 3 then $
				fabref = rebin(*fab, fabsize[0], fabsize[1], $
					       fabsize[2], /sample)
		endif else begin
			if amr.ndim eq 1 then $
				fabref = refine(*fab, fabsize[0]) $
;				fabref = rebin(*fab, fabsize[0]) $
			else if amr.ndim eq 2 then $
				fabref = refine(*fab, fabsize[0], fabsize[1]) $
;				fabref = rebin(*fab, fabsize[0], fabsize[1]) $
			else if amr.ndim eq 3 then $
				fabref = refine(*fab, fabsize[0], fabsize[1], $
					fabsize[2])
;				fabref = rebin(*fab, fabsize[0], $
;					fabsize[1], fabsize[2])
		endelse

		; construct max and min within part of data
		; box that overlaps this fab
		datamin = fablo - idxmin
		datamin = (datamin ge 0) * datamin
		datamax = fabhi - idxmin
		datamax = (datamax lt datasize) * datamax + $
			(datamax ge datasize) * (datasize - 1)

		; construct max and min within this fab that
		; overlaps the data box
		fabmin = idxmin - fablo
		fabmin = (fabmin ge 0) * fabmin
		fabmax = idxmax - fablo
		fabmax = (fabmax lt fabsize) * fabmax + $
			(fabmax ge fabsize) * (fabsize - 1)


		; insert fab into data array
		if amr.ndim eq 1 then data[datamin[0] : datamax[0]] = $
			fabref[fabmin[0] : fabmax[0]] $
		else if amr.ndim eq 2 then $
			data[datamin[0] : datamax[0], $
			     datamin[1] : datamax[1]] = $
			fabref[fabmin[0] : fabmax[0], $
			       fabmin[1] : fabmax[1]] $
		else if amr.ndim eq 3 then $
			data[datamin[0] : datamax[0], $
			     datamin[1] : datamax[1], $
			     datamin[2] : datamax[2]] = $
			fabref[fabmin[0] : fabmax[0], $
			       fabmin[1] : fabmax[1], $
			       fabmin[2] : fabmax[2]]

	endfor

endfor

; return data
return, data
end

