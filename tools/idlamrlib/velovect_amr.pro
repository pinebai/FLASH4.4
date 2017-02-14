pro velovect_amr, amr_vx, amr_vy, plane, slice, 		$
	anisotropic=anisotropic, overplot=overplot,		$
	verbose=verbose, length=length,	color=color,		$
	maxlevel=maxlevel, interval=interval,			$
	noscale=noscale, xrange=xrange, yrange=yrange,		$
	scalestr=scalestr, snaptogrid=snaptogrid,		$
	xmargin=xmargin, ymargin=ymargin,			$
	charsize=charsize, scalefac=scalefac, _extra = extra

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


; This routine produces a velocity field of a pair of amr objects.
; The variables amr_vx and amr_vy are interpreted as amr objects
; containing the "x" and "y" components of the velocity field.
; The plot is sliced along the plane plane, so, for example, plane = 0
; produces a yz slice, plane = 1 produces xz, and plane = 2 produces
; xy. Slice  specifies the physical location of the slice in data units;
; in other words, setting plane = 2, slice = 1.0e17 shows the xy plane
; sliced at z = 1.0e17. Both plane and slice are ignored for 2d amr
; objects.
; Effects of keywords (other than standard IDL graphics keywords):
;   anisotropic: plots use isotropic axes by default. If anisotropic is
;	set, axes are allowed to be anisotropic and the plot will fill
;	the available area
;   verbose: make the routine print out some data about the amr object
;   maxlevel: causes the routine to only raster up to the specified
;	level of refinement. This can be useful in preventing the
;	production of over large postscript files.
;   overplot: causes the routine to act as an overplot, i.e. to use
;	an existing set of axes rather than drawing new one, and to
;	overwrite rather than erase any existing plots
;   length: the default length of 1.0 causes the highest velocity
;	(in the xy plane) present in the data to correspond to an
;	arrow with a length of one cell on the finest level used.
;	Changing the value from 1.0 shrinks or lengthens it
;	proportionately.
;   interval: this causes the routine to only plot every (interval)
;	cell, rather than every cell
;   noscale: by default, the routine draws an arrow to show the velocity
;	scale. If noarrow is set, it does not.
;   scalestr: this parameter is interpreted as a string and added to the
;	label for the scale arrow. This is useful for adding units.
;   snaptogrid: if set, this keyword adjusts the x and y range so that
;	they fall exactly at coarse cell boundaries.

; get keywords
if not keyword_set(anisotropic) then isotropic=1 else isotropic=0
if n_elements(maxlevel) eq 0 then maxlevel=amr_vx.maxlevel
if not keyword_set(overplot) then overplot = 0
if not keyword_set(length) then length=1.0
if n_elements(color) eq 0 then color=!d.table_size-1
if not keyword_set(interval) then interval=1
if not keyword_set(charsize) then charsize=!p.charsize
if charsize eq 0 then charsize=1.0
if not keyword_set(scalefac) then scalefac=1.0
if not keyword_set(xmargin) then xmargin=!x.margin
if not keyword_set(ymargin) then ymargin=!y.margin
mathprec = machar(/double)
maskval = sqrt(mathprec.xmax)/3.0

; check for consistency of amr_vx and amr_vy
if amr_vx.name ne amr_vy.name then begin
	print, 'Error: amr_vx and amr_vy must be from the same plot file'
	return
endif

; figure out which direction we're slicing along
if amr_vx.ndim eq 2 then begin
	dim1=0
	dim2=1
endif else if amr_vx.ndim eq 3 then begin
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

; get some information about the graphics environment
tempmulti = !p.multi
tempmulti[1:3] = tempmulti[1:3] + (tempmulti[1:3] eq 0)
if overplot then begin
	tempmulti[0] = tempmulti[0] - 1
	if tempmulti[0] lt 0 then tempmulti[0] = $
		tempmulti[1]*tempmulti[2]*tempmulti[3] - 1
endif

; get the row and column for this plot	
if tempmulti[4] eq 0 then begin
	plotrow = tempmulti[2] - tempmulti[0] / tempmulti[1] - 1
	plotcol = tempmulti[0] mod tempmulti[1]
endif else begin
	plotrow = tempmulti[2] - (tempmulti[0] mod tempmulti[2]) - 1
	plotcol = tempmulti[0] / tempmulti[2]
endelse

; Get plot margins in normal coordinates. If we're on a device that uses
; windows, make sure a valid window is open before attempting this.
devname = !d.name
devname = strmid(devname, 0, 3)
devname = strupcase(devname)
case devname of
	'MAC': if !d.window eq -1 then window, 0
	'WIN': if !d.window eq -1 then window, 0
	'X': if !d.window eq -1 then window, 0
	ELSE:
endcase
xmargindev = xmargin * !d.x_ch_size
ymargindev = ymargin * !d.y_ch_size
marginnorm = convert_coord(xmargindev, ymargindev, /device, /to_normal)

; set range of plot range
if overplot then begin
	if n_elements(xrange) eq 0 then xlim=!x.crange else xlim=xrange
	if n_elements(yrange) eq 0 then ylim=!y.crange else ylim=yrange
endif else begin

	; set range of plot range
	if not keyword_set(xrange) then $
		xlim=[amr_vx.boxmin[dim1], amr_vx.boxmax[dim1]] $
	else xlim=xrange
	if not keyword_set(yrange) then $
		ylim=[amr_vx.boxmin[dim2], amr_vx.boxmax[dim2]] $
	else ylim=yrange

	; if asked to do so, snap range to grid
	if keyword_set(xrange) and keyword_set(snaptogrid) then begin
		xlim[0] = floor( (xlim[0] - amr_vx.boxmin[dim1]) / $
				 amr_vx.gridspacing[dim1,0] ) * $
			  amr_vx.gridspacing[dim1, 0] + $
			  amr_vx.boxmin[dim1]
		xlim[1] = ceil( (xlim[1] - amr_vx.boxmin[dim1]) / $
				 amr_vx.gridspacing[dim1,0] ) * $
			  amr_vx.gridspacing[dim1, 0] + $
			  amr_vx.boxmin[dim1]
		if keyword_set(verbose) then print, 'New x range = ', xlim
	endif
	if keyword_set(yrange) and keyword_set(snaptogrid) then begin
		ylim[0] = floor( (ylim[0] - amr_vx.boxmin[dim2]) / $
				 amr_vx.gridspacing[dim2,0] ) * $
			  amr_vx.gridspacing[dim2, 0] + $
			  amr_vx.boxmin[dim2]
		ylim[1] = ceil( (ylim[1] - amr_vx.boxmin[dim2]) / $
				 amr_vx.gridspacing[dim2,0] ) * $
			  amr_vx.gridspacing[dim2, 0] + $
			  amr_vx.boxmin[dim2]
		if keyword_set(verbose) then print, 'New y range: ', ylim
	endif
endelse

; create plot axes if needed
if keyword_set(noarrow) then begin
	if not overplot then $
		plot, [0,0], [0,0], /nodata, isotropic=isotropic, $
			xrange = xlim, yrange = ylim, $
			charsize=charsize, xmargin=xmargin, ymargin=ymargin, $
			_extra = extra
endif else begin
	; in the case with a scale arrow, we must deal with
	; multiple plots (!p.multi settings) manually

	; set the plot box position
	plotpos = fltarr(4)
	plotpos[0] = marginnorm[0,0] + float(plotcol) / tempmulti[1]
	plotpos[1] = marginnorm[1,0] + float(plotrow) / tempmulti[2]
	plotpos[2] = (plotcol+0.7) / tempmulti[1] - marginnorm[0,1]
	plotpos[3] = float(plotrow+1) / tempmulti[2] - marginnorm[1,1]

	if not overplot then $
		plot, [0,0], [0,0], /nodata, _extra = extra, $
			xrange = xlim, yrange = ylim, $
			isotropic = isotropic, $
			charsize=charsize, position = plotpos

endelse

; turn off device decomposition if applicable
if (strupcase(!d.name) eq 'MAC') or (strupcase(!d.name) eq 'WIN') or $
   (strupcase(!d.name) eq 'X') then begin
	device, get_decomposed = saveDecomposed
	device, decomposed = 0
endif

; figure out the length of the highest velocity present in the problem
vxmax = max_amr(amr_vx, vxmin)
vymax = max_amr(amr_vy, vymin)
vmax = max([abs(vxmax), abs(vxmin), abs(vymax), abs(vymin)])

; We need to keep a list of what index ranges are covered by fine grids,
; so we'll know not to put coarse grids on top of them. Since we don't
; know a priori how many fine grids will be in the part of the grid we're
; plotting, we create a linked list for each level to hold these.
if maxlevel ne 0 then begin
	amr_masklist_struct = { xlo:dblarr(2), xhi:dblarr(2),$
				next:ptr_new(), last:ptr_new() }
	maskptr = ptrarr(maxlevel)
endif

; go through amr levels
for n=maxlevel, 0, -1 do begin

	; set the refinement ratio between this level and the finest level
	refratio = 1
	for m=n, maxlevel-1 do refratio = refratio * amr_vx.refratio[n]

	; set grid spacing
	dx = [amr_vx.levels[n].gridspacing[dim1], $
	      amr_vx.levels[n].gridspacing[dim2]]

	; loop through fabs
	for m=0, amr_vx.levels[n].nfab-1 do begin

		; check that this fab is within the range we want
		; if not, skip it

		; check the slice direction first. If the requested slice
		; doesn't cut through this fab, then go to the next one.
		if amr_vx.ndim eq 3 then begin
			if (slice lt (*amr_vx.levels[n].fabptr)[m].xlo[dim3]) $
			or (slice ge (*amr_vx.levels[n].fabptr)[m].xhi[dim3]) $
			then continue
		endif

		; now check in the other direction. First figure this out
		; in data coordinates
		fabmin = [(*amr_vx.levels[n].fabptr)[m].xlo[dim1], $
			  (*amr_vx.levels[n].fabptr)[m].xlo[dim2]]
		fabmax = [(*amr_vx.levels[n].fabptr)[m].xhi[dim1], $
			  (*amr_vx.levels[n].fabptr)[m].xhi[dim2]]
		fabidxmin = [(*amr_vx.levels[n].fabptr)[m].idxlo[dim1], $
			     (*amr_vx.levels[n].fabptr)[m].idxlo[dim2]]
		fabidxmax = [(*amr_vx.levels[n].fabptr)[m].idxhi[dim1], $
			     (*amr_vx.levels[n].fabptr)[m].idxhi[dim2]]
		datalo=fltarr(2)
		datahi=fltarr(2)
		datalo[0] = (fabmin[0] gt xlim[0]) * fabmin[0] + $
			(fabmin[0] le xlim[0]) * xlim[0]
		datalo[1] = (fabmin[1] gt ylim[0]) * fabmin[1] + $
			(fabmin[1] le ylim[0]) * ylim[0]
		datahi[0] = (fabmax[0] le xlim[1]) * fabmax[0] + $
			(fabmax[0] gt xlim[1]) * xlim[1]
		datahi[1] = (fabmax[1] le ylim[1]) * fabmax[1] + $
			(fabmax[1] gt ylim[1]) * ylim[1]

		; if there is no overlap between this fab and the
		; desired region, then move to the next fab
		if (datalo[0] gt fabmax[0]) or (datahi[0] lt fabmin[0]) or $
		   (datalo[1] gt fabmax[1]) or (datahi[1] lt fabmin[1]) $
			then continue

		; If we're here, we're using this fab, so put it
		; into the masking list unless we're on the 0th
		; level. (The 0th level can't mask anything.)
		if n ne 0 then begin
			if not ptr_valid(maskptr[n-1]) then begin
				; initialize list if this is the first element
				maskptr[n-1] = ptr_new(amr_masklist_struct)
				current = maskptr[n-1]
				last = ptr_new()
			endif else begin
				; create the next element and move the
				; last pointer
				next = ptr_new(amr_masklist_struct)
				(*current).next = next
				last = current
				current = next
			endelse

			; store values in this list element
			(*current).xlo = fabmin
			(*current).xhi = fabmax
			(*current).last = last

		endif

		; figure out the index to use for the 3rd coordinate, if
		; applicable
		if amr_vx.ndim eq 3 then begin
			fabmin3 = (*amr_vx.levels[n].fabptr)[m].xlo[dim3]
			fabmax3 = (*amr_vx.levels[n].fabptr)[m].xhi[dim3]
			fabidxmin3 = (*amr_vx.levels[n].fabptr)[m].idxlo[dim3]
			fabidxmax3 = (*amr_vx.levels[n].fabptr)[m].idxhi[dim3]
			dx3 = (fabmax3-fabmin3) / (fabidxmax3-fabidxmin3+1)
			idx3 = round( (slice-fabmin3)/dx3 - 0.5 )
			if idx3 lt 0 then idx3 = 0
			if idx3 gt (fabidxmax3-fabidxmin3) then $
				idx3 = fabidxmax3-fabidxmin3
		endif

		; set up list of indices to grab
		xidxlist = fabidxmin[0] + $
			indgen(fabidxmax[0]-fabidxmin[0]+1)
		xidxsub = where(xidxlist mod interval eq 0)
		if total(xidxsub eq -1) then continue else $
			xidxlist = xidxlist[xidxsub]
		xidxlist = xidxlist - fabidxmin[0]
		yidxlist = fabidxmin[1] + $
			indgen(fabidxmax[1]-fabidxmin[1]+1)
		yidxsub = where(yidxlist mod interval eq 0)
		if total(yidxsub eq -1) then continue else $
			yidxlist = yidxlist[yidxsub]
		yidxlist = yidxlist - fabidxmin[1]

		; construct arrays giving the x and y coordinates of
		; the cell centers
		xcell = fabmin[0] + dx[0] * (xidxlist + 0.5)
		ycell = fabmin[1] + dx[1] * (yidxlist + 0.5)

		; clip list of indices to grab to ensure that it's in
		; the requested region
		xclip = where((xcell ge datalo[0]) and (xcell le datahi[0]))
		yclip = where((ycell ge datalo[1]) and (ycell le datahi[1]))
		if total(xclip) eq -1 then continue else begin
			xidxlist = xidxlist[xclip]
			xcell = xcell[xclip]
		endelse
		if total(yclip) eq -1 then continue else begin
			yidxlist = yidxlist[yclip]
			ycell = ycell[yclip]
		endelse

		; grab the velocity field
		if amr_vx.ndim eq 2 then begin
			;2d case
			dattmp = *(*amr_vx.levels[n].fabptr)[m].dataptr
			datsz = size(dattmp)
			dattmp = reform(dattmp, 1, datsz[1], datsz[2], $
					/overwrite)
			u = dattmp[0, xidxlist, yidxlist]
			dattmp = *(*amr_vy.levels[n].fabptr)[m].dataptr
			datsz = size(dattmp)
			dattmp = reform(dattmp, 1, datsz[1], datsz[2], $
					/overwrite)
			v = dattmp[0, xidxlist, yidxlist]
		endif else begin
			; 3 possible 3d cases
			if plane eq 0 then begin
			   u = (*(*amr_vx.levels[n].fabptr)[m].dataptr) $ 
				[idx3, xidxlist, yidxlist]
			   v = (*(*amr_vy.levels[n].fabptr)[m].dataptr) $ 
				[idx3, xidxlist, yidxlist]
			endif else if plane eq 1 then begin
			   u = (*(*amr_vx.levels[n].fabptr)[m].dataptr) $ 
				[xidxlist, idx3, yidxlist]
			   v = (*(*amr_vy.levels[n].fabptr)[m].dataptr) $ 
				[xidxlist, idx3, yidxlist]
			endif else begin
			   u = (*(*amr_vx.levels[n].fabptr)[m].dataptr) $ 
				[xidxlist, yidxlist, idx3]
			   v = (*(*amr_vy.levels[n].fabptr)[m].dataptr) $ 
				[xidxlist, yidxlist, idx3]
			endelse
		endelse

		u = reform(u, n_elements(xidxlist), $
			   n_elements(yidxlist))
		v = reform(v, n_elements(xidxlist), $
			   n_elements(yidxlist))

		; Mask regions of the velocity field covered by regions
		; in the mask list. Do this by setting velocities at
		; these points to a dummy value.
		for l=n,maxlevel-1 do begin
			checkptr=maskptr[l]

			; go through mask list on this level
			while ptr_valid(checkptr) do begin

				; grab edges of mask region
				masklo = (*checkptr).xlo
				maskhi = (*checkptr).xhi

				; see which regions of the velocity
				; field are within this region
				xmask = (xcell gt masklo[0]) and $
					(xcell lt maskhi[0])
				ymask = (ycell gt masklo[1]) and $
					(ycell lt maskhi[1])

				; construct the mask
				mask = bytarr(n_elements(xmask), $
					      n_elements(ymask))
				for i=0,n_elements(xmask)-1 do $
				   mask[i,*] = xmask[i] and ymask

				; mask out regions of velocity arrays
				v = mask * maskval + (1 - mask) * v
				u = mask * maskval + (1 - mask) * u

				; point to next mask
				checkptr = (*checkptr).next
			endwhile
		endfor

		; if this entire fab is masked, then go to next fab
		if total(where(u ne maskval)) eq -1 then continue

		; get the largest velocity present in this field
		utmp = (u lt maskval) * u
		vtmp = (v lt maskval) * v
		vmaxlocal = max([max(abs(utmp)), max(abs(vtmp))])

		; display vector field
		velovect_mrk, u, v, xcell, ycell, /overplot, $
			length = length * vmaxlocal / $
				 (vmax * refratio * interval), $
			missing=maskval, color=color, $
			cellsize=dx*interval
	endfor
endfor

; free the mask list
for n=1, maxlevel do begin
	current = maskptr[n-1]
	if not ptr_valid(current) then continue
	while ptr_valid((*current).next) do current = (*current).next
	while ptr_valid(current) do begin
		prev = (*current).last
		ptr_free, current
		current = prev
	endwhile
endfor

; draw a scale arrow if requested
if not keyword_set(noscale) then begin

	; figure out the length of the scale arrow
	arrowlength = length * max(amr_vx.gridspacing[[dim1,dim2],maxlevel]) $
		* scalefac
	normlength = convert_coord([0,0], [0,arrowlength], /data, /to_normal)
	normlength = normlength[1,1]-normlength[1,0]

	; figure out the position of the arrow
	arrowpos = fltarr(2)
	arrowpos[0] = (plotcol + 0.77) / tempmulti[1]
	arrowpos[1] = plotpos[1]

	; draw the arrow
	arrow, arrowpos[0], arrowpos[1], arrowpos[0], $
		arrowpos[1]+normlength, /normalized, hsize = -0.2

	; put a numerical label on the scale arrow
	if (abs(vmax) lt 1.0e4) and (abs(vmax) gt 1.0e-4) then $
		arrowstr = string(vmax*scalefac, format = '(f8.2)') $
	else $
		arrowstr = string(vmax*scalefac, format = '(e10.2)')
	arrowstr = strtrim(arrowstr,2)
	arrowstrpos = convert_coord(arrowpos[0], arrowpos[1], $
				    /normal, /to_device)
	arrowstrpos = arrowstrpos[0:1]
	arrowstrpos[0] = arrowstrpos[0] - $
			 strlen(arrowstr) * !d.x_ch_size *charsize / 2.0
	arrowstrpos[1] = arrowstrpos[1] - 1.5 * !d.y_ch_size * charsize
	xyouts, arrowstrpos[0], arrowstrpos[1], arrowstr, charsize=charsize, $
		/device
	if keyword_set(scalestr) then begin
		scalestrpos = convert_coord(arrowpos[0], arrowpos[1],$
			 /normal, /to_device)
		scalestrpos = scalestrpos[0:1]
		scalestrpos[0] = scalestrpos[0] - $
			 strlen(scalestr) * !d.x_ch_size *charsize / 2.0
		scalestrpos[1] = scalestrpos[1] - 2.7 * !d.y_ch_size * charsize
		xyouts, scalestrpos[0], scalestrpos[1], scalestr, $
			charsize=charsize, /device
	endif
endif

; turn on device decomposition if applicable
if (strupcase(!d.name) eq 'MAC') or (strupcase(!d.name) eq 'WIN') or $
   (strupcase(!d.name) eq 'X') then begin
	device, decomposed = saveDecomposed
endif

return
end



