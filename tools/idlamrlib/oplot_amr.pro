pro oplot_amr, amr, dir, yval, zval, xrange=xrange, $
	maxlevel=maxlevel, minlevel=minlevel, xvec=xvec, yvec=yvec, $
	_extra=_extra

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


; see plot_amr for a description -- this routine is to plot_amr as
; oplot is to plot. The only distinctive feature is that the xrange
; keyword is accepted. Its effect is to specify the range of the
; data to be plotted, and it has no effect on the axes generated.

; get keywords
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel
if n_elements(minlevel) eq 0 then minlevel=0

; figure out which direction we're taking the line in
if amr.ndim eq 2 then begin
	if dir eq 0 then dim1=1 else $
	if dir eq 1 then dim1=0 else begin
		print, 'Error: direction must be 0 or 1'
		return
	endelse
endif else if amr.ndim eq 3 then begin
	if dir eq 0 then begin
		dim1=1
		dim2=2
	endif else if dir eq 1 then begin
		dim1=0
		dim2=2
	endif else if dir eq 2 then begin
		dim1=0
		dim2=1
	endif else begin
		print, 'Error: dir must be 0, 1, or 2'
		return
	endelse
endif else if amr.ndim ne 1 then begin
	print, 'Error: number of dimensions must be 1, 2, or 3'
	return
endif

; set x range
if not keyword_set(xrange) then xrange=[amr.boxmin[dir], amr.boxmax[dir]]

; create vector to hold data -- it should be as long as the maximum
; number of fine cells in the direction we're looking. Note that we need
; one extra element if the data is edge centered rather than cell
; centered.
xvec = dblarr(amr.idxhi[dir, maxlevel] - amr.idxlo[dir, maxlevel] + 1 + $
		amr.idxtype[dir])
yvec = dblarr(amr.idxhi[dir, maxlevel] - amr.idxlo[dir, maxlevel] + 1 + $
		amr.idxtype[dir])
edgecenter = amr.idxtype[dir]
counter = 0

; create a vector to hold information on where fine cells exist, so we
; should skip coarse cells in that region
fineflag = bytarr(amr.idxhi[dir, maxlevel] - $
		  amr.idxlo[dir, maxlevel] + 1 + edgecenter)
nfineflag = n_elements(fineflag)

; loop through levels, starting at the finest level
for n=maxlevel, minlevel, -1 do begin

	; get grid spacing in on this level
	dx = amr.levels[n].gridspacing

	; loop through fabs
	for m=0, amr.levels[n].nfab-1 do begin

		; check that this fab is within the range we want along
		; the direction of line. If not, skip it.
		if (xrange[0] gt (*amr.levels[n].fabptr)[m].xhi[dir]) or $
		   (xrange[1] lt (*amr.levels[n].fabptr)[m].xlo[dir]) $
			then continue

		; see if this fab intersects the region we want. If not,
		; skip it.
		if amr.ndim gt 1 then $
			if (yval lt (*amr.levels[n].fabptr)[m].xlo[dim1]) or $
			   (yval gt (*amr.levels[n].fabptr)[m].xhi[dim1]) $
				then continue;
		if amr.ndim gt 2 then $
			if (zval lt (*amr.levels[n].fabptr)[m].xlo[dim2]) or $
			   (zval gt (*amr.levels[n].fabptr)[m].xhi[dim2]) $
				then continue;

		; figure out which indices we want in perpendicular directions
		if amr.ndim eq 2 then begin
			if (not edgecenter) then $
				yidx = fix( (yval - $
				       (*amr.levels[n].fabptr)[m].xlo[dim1]) $
			       	       / dx[dim1] - 0.5) $
			else $
				yidx = fix( (yval - $
				       (*amr.levels[n].fabptr)[m].xlo[dim1]) $
			       	       / dx[dim1] )
		endif else begin
			if (not edgecenter) then begin
				yidx = fix( (yval - $
				       (*amr.levels[n].fabptr)[m].xlo[dim1]) $
			       	       / dx[dim1] - 0.5)
				zidx = fix( (zval - $
				       (*amr.levels[n].fabptr)[m].xlo[dim2]) $
			       	       / dx[dim2] - 0.5)
			endif else begin
				yidx = fix( (yval - $
				       (*amr.levels[n].fabptr)[m].xlo[dim1]) $
			       	       / dx[dim1] )
				zidx = fix( (zval - $
				       (*amr.levels[n].fabptr)[m].xlo[dim2]) $
			       	       / dx[dim2] )
			endelse
		endelse

		; create a mask to indicate regions where fine cells
		; have already filled in the values and we should
		; therefore not use coarse cells
		if edgecenter then begin
			coarseflag = rebin(fineflag[0:nfineflag-2], $
					   amr.idxhi[dir,n] - $
					   amr.idxlo[dir,n] + 1, $
					   /sample) 
			coarseflag = [coarseflag, fineflag[nfineflag-1]]
		endif else begin
			coarseflag = rebin(fineflag[0:nfineflag-1], $
					   amr.idxhi[dir,n] - $
					   amr.idxlo[dir,n] + 1, $
					   /sample)
		endelse
		coarseflag = $
		   coarseflag[ (*amr.levels[n].fabptr)[m].idxlo[dir] : $
			       (*amr.levels[n].fabptr)[m].idxhi[dir] ]
		ptsToUse = where(coarseflag eq 0)

		; if no points to use, move to next fab
		if ptsToUse[0] eq -1 then continue

		; store the positions of these points
		npts = n_elements(ptsToUse)
		if edgecenter then $
			xvec[counter : counter+npts-1] = $
				(*amr.levels[n].fabptr)[m].xlo[dir] + $
				dx[dir] * ptsToUse $
		else $
			xvec[counter : counter+npts-1] = $
				(*amr.levels[n].fabptr)[m].xlo[dir] + $
				dx[dir] * (ptsToUse + 0.5)

		; store the values at these points
		if amr.ndim eq 1 then begin
		   yvec[counter : counter+npts-1] = $
		      (*(*amr.levels[n].fabptr)[m].dataptr)[ptsToUse]
		endif else if amr.ndim eq 2 then begin
		   if dir eq 0 then $
		      yvec[counter : counter+npts-1] = $
			(*(*amr.levels[n].fabptr)[m].dataptr)[ptsToUse,yidx] $
		   else $
		      yvec[counter : counter+npts-1] = $
			(*(*amr.levels[n].fabptr)[m].dataptr)[yidx,ptsToUse]
		endif else begin
		   if dir eq 0 then $
		      yvec[counter : counter+npts-1] = $
			(*(*amr.levels[n].fabptr)[m].dataptr) $
			  [ptsToUse,yidx,zidx] $
		   else if dir eq 1 then $
		      yvec[counter : counter+npts-1] = $
			(*(*amr.levels[n].fabptr)[m].dataptr) $
			  [yidx,ptsToUse,zidx] $
		   else $
		      yvec[counter : counter+npts-1] = $
			(*(*amr.levels[n].fabptr)[m].dataptr) $
			  [yidx,zidx,ptsToUse]
		endelse

		; move point counter
		counter = counter + npts

		; record that we have fine cells at these indices
		thisCoarseFlag = bytarr(amr.idxhi[dir, n] - $
					amr.idxlo[dir, n] + 1 + edgecenter)
		thisCoarseFlag[ (*amr.levels[n].fabptr)[m].idxlo[dir] : $
			  (*amr.levels[n].fabptr)[m].idxhi[dir] ] = 1
		thisFineFlag = rebin(thisCoarseFlag[0: $
					amr.idxhi[dir, n] - $
					amr.idxlo[dir, n]], $
				amr.idxhi[dir, maxlevel] - $
				amr.idxlo[dir, maxlevel] + 1, /sample)
		if edgecenter then $
			thisFineFlag = [thisFineFlag, $
					thisCoarseFlag[amr.idxhi[dir, n] - $
					amr.idxlo[dir, n] + 1]]
		fineFlag = fineFlag or thisFineFlag

	endfor
endfor

; extract the part of the data arrays that contains actual data
xvec = xvec[0:counter-1]
yvec = yvec[0:counter-1]

; sort the data
sortlist = sort(xvec)
xvec = xvec[sortlist]
yvec = yvec[sortlist]

; plot the data
oplot, xvec, yvec, _extra=_extra

end






