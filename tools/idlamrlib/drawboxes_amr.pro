pro drawboxes_amr, amr, plane, slice, anisotropic=anisotropic, $
	xrange=xrange, yrange=yrange, snaptogrid=snaptogrid, $
	overplot=overplot, maxlevel=maxlevel, verbose=verbose, $
	minlevel=minlevel, linestyle=linestyle, psym=psym, $
	symsize=symsize, thick=thick, color=color, _extra=_extra

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


; this routine draws the amr box boundaries seen at a particular
; plane and slice. The keywords are the same as for raster_amr.

; get keywords
if not keyword_set(anisotropic) then isotropic=1 else isotropic=0
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel
if n_elements(minlevel) eq 0 then minlevel=0
if n_elements(linestyle) eq 0 then linestyle=!p.linestyle
if n_elements(psym) eq 0 then psym=!p.psym
if n_elements(symsize) eq 0 then symsize=1.0
if n_elements(thick) eq 0 then thick=!p.thick

; Scalable pixels generally equal ps, in which case we want black (= 0)
; because the background is white. Otherwise use top of the color table.
scalablePixels = !d.flags mod 2
if not scalablePixels then begin
	if n_elements(color) eq 0 then color=!d.table_size-1
endif else begin
	if n_elements(color) eq 0 then color=0
endelse

; figure out which direction we're slicing along
if amr.ndim eq 2 then begin
	dim1=0
	dim2=1
endif else if amr.ndim eq 3 then begin
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
		return
	endelse
endif else begin
	print, 'Error: amr object must have 2 or 3 dimensions.'
	return
endelse

; set range of plot if this is not an overplot
if not keyword_set(overplot) then begin
	if not keyword_set(xrange) then $
		xlim=[amr.boxmin[dim1], amr.boxmax[dim1]] $
	else xlim=xrange
	if not keyword_set(yrange) then $
		ylim=[amr.boxmin[dim2], amr.boxmax[dim2]] $
	else ylim=yrange

	; if asked to do so, snap range to grid
	if keyword_set(xrange) and keyword_set(snaptogrid) then begin
		xlim[0] = floor( (xlim[0] - amr.boxmin[dim1]) / $
				 amr.gridspacing[dim1,0] ) * $
			  amr.gridspacing[dim1, 0] + $
			  amr.boxmin[dim1]
		xlim[1] = ceil( (xlim[1] - amr.boxmin[dim1]) / $
				 amr.gridspacing[dim1,0] ) * $
			  amr.gridspacing[dim1, 0] + $
			  amr.boxmin[dim1]
		if keyword_set(verbose) then print, 'New x range = ', xlim
	endif
	if keyword_set(yrange) and keyword_set(snaptogrid) then begin
		ylim[0] = floor( (ylim[0] - amr.boxmin[dim2]) / $
				 amr.gridspacing[dim2,0] ) * $
			  amr.gridspacing[dim2, 0] + $
			  amr.boxmin[dim2]
		ylim[1] = ceil( (ylim[1] - amr.boxmin[dim2]) / $
				 amr.gridspacing[dim2,0] ) * $
			  amr.gridspacing[dim2, 0] + $
			  amr.boxmin[dim2]
		if keyword_set(verbose) then print, 'New y range: ', ylim
	endif

	; create axes unless this is an overplot
	plot, [0,0], [0,0], /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra
endif

; Turn off device decomposition for scalable pixels
if not scalablePixels then begin
        if !Version.Release ge 5.2 then device, get_decomposed=thisdecomposed
        device, decomposed=0
endif

; now go through the amr structure
for n=minlevel, maxlevel do begin
	for m=0, amr.levels[n].nfab-1 do begin

		; check if this fab intersects the image plane
		if amr.ndim eq 3 then begin
			if (slice lt (*amr.levels[n].fabptr)[m].xlo[dim3]) $
			or (slice gt (*amr.levels[n].fabptr)[m].xhi[dim3]) $
			then continue
		endif

		; we could check if this fab intersects the right region of
		; the image plane, but it's easier to just draw everything
		; and use the clip keyword to cut out what we don't want

		; get limits on fab
		fabmin = [(*amr.levels[n].fabptr)[m].xlo[dim1], $
			  (*amr.levels[n].fabptr)[m].xlo[dim2]]
		fabmax = [(*amr.levels[n].fabptr)[m].xhi[dim1], $
			  (*amr.levels[n].fabptr)[m].xhi[dim2]]

		; draw the box
		plots, [fabmin[0], fabmin[0], fabmax[0], fabmax[0], $
			fabmin[0]], $
		       [fabmin[1], fabmax[1], fabmax[1], fabmin[1], $
			fabmin[1]], noclip=0, thick=thick, $
			linestyle=linestyle, psym=psym, $
			color=color, symsize=symsize

	endfor
endfor

; turn decomposition back to previous setting
if not scalablePixels then begin
        if !Version.Release ge 5.2 then device, decomposed=thisdecomposed
endif

return
end




