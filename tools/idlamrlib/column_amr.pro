pro column_amr, amr, plane, anisotropic=anisotropic,		$
	xrange=xrange, yrange=yrange, nocolorbar=nocolorbar,	$
	verbose=verbose, scalerange=scalerange,			$
	maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
	pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
	smeardist=smeardist, charsize=charsize,			$
	minlevel=minlevel, fitcolorbar=fitcolorbar,		$
	suppresscolorbar=suppresscolorbar, ymargin=ymargin,	$
	skinnycolorbar=skinnycolorbar, xmargin=xmargin,		$
	imgout=imgout, floor=floor, noplot=noplot,		$
	_extra = extra

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

; get keywords
if not keyword_set(anisotropic) then isotropic=1 else isotropic=0
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel
if n_elements(minlevel) eq 0 then minlevel=0
if not keyword_set(charsize) then charsize=!p.charsize
if not keyword_set(xmargin) then xmargin=!x.margin
if not keyword_set(ymargin) then ymargin=!y.margin

; figure out which direction we're slicing along
if amr.ndim lt 3 then begin
	print, 'Error: amr object must have 3 dimensions.'
	return
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
	return
endelse

; set range of plot range
if not keyword_set(xrange) then xlim=[amr.boxmin[dim1], amr.boxmax[dim1]] $
else xlim=xrange
if not keyword_set(yrange) then ylim=[amr.boxmin[dim2], amr.boxmax[dim2]] $
else ylim=yrange

; if asked to do so, snap range to grid
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

; get some information about the graphics environment
bangmulti = !p.multi
tempmulti = !p.multi
tempmulti[1:3] = tempmulti[1:3] + (tempmulti[1:3] eq 0)
xmargindev = xmargin * !d.x_ch_size
ymargindev = ymargin * !d.y_ch_size

; Convert data to normal coordinates. If we're on a device that uses
; windows, make sure a valid window is open before attempting this.
if not keyword_set(noplot) then begin
    devname = !d.name
    devname = strmid(devname, 0, 3)
    devname = strupcase(devname)
    case devname of
	'MAC': if !d.window eq -1 then window, 0
	'WIN': if !d.window eq -1 then window, 0
	'X': if !d.window eq -1 then window, 0
	ELSE:
    endcase
    marginnorm = convert_coord(xmargindev, ymargindev, /device, /to_normal)
endif

; see how many colors are available
ncolors = !d.table_size

; draw axes
if not keyword_set(noplot) then begin
    if keyword_set(nocolorbar) then $
	plot, [0,0], [0,0], /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, charsize=charsize, _extra = extra $
    else begin
	; in the case with a colorbar, we must deal with multiple plots
	; (!p.multi settings) manually

	; set the row and column for this plot
	if tempmulti[4] eq 0 then begin
		plotrow = tempmulti[2] - tempmulti[0] / tempmulti[1] - 1
		plotcol = tempmulti[0] mod tempmulti[1]
	endif else begin
		plotrow = tempmulti[2] - (tempmulti[0] mod tempmulti[2]) - 1
		plotcol = tempmulti[0] / tempmulti[2]
	endelse

	if not keyword_set(pltsize) then pltsize=0.7

	; set the plot box position
	plotpos = fltarr(4)
	plotpos[0] = marginnorm[0,0] + float(plotcol) / tempmulti[1]
	plotpos[1] = marginnorm[1,0] + float(plotrow) / tempmulti[2]
	plotpos[2] = (plotcol+0.7) / tempmulti[1] - marginnorm[0,1]
	plotpos[3] = float(plotrow+1) / tempmulti[2] - marginnorm[1,1]

	; draw the plot box
	plot, [0,0], [0,0], /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra, $
		position = plotpos, charsize=charsize
    endelse
endif

; set image size equal to number of cells on maxlevel
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

; What we do here depends on whether we have a device with scalable
; pixels or not. If we have scalable pixels, we can just construct
; an image box where each cell of the amr grid on the finest level
; represents a single pixel. We then rely on the scaling of pixels
; to make that image correctly overlay the plot window. If we don't
; have scalable pixels, the number of cells per pixel will in be a
; complicated function of the size of the graphics device and the
; resolution of the amr grid. For reference, X windows does not have
; scalable pixels, postscript does.

scalablePixels = (!d.flags mod 2) or keyword_set(noplot)

if not scalablePixels then begin
	; we do not have scalable pixels

	; get max and min for image box in pixels
	imgdisplo = convert_coord(xlim[0], ylim[0], /data, /to_device)
	imgdisphi = convert_coord(xlim[1], ylim[1], /data, /to_device)
	imgdisplo=round(imgdisplo[0:1])
	imgdisphi=round(imgdisphi[0:1])
	if (imgdisphi[0] - imgdisplo[0]) lt 2 then begin
		print, 'Error: requested x axis range is too small.'
		return
	endif
	if (imgdisphi[1] - imgdisplo[1]) lt 2 then begin
		print, 'Error: requested y axis range is too small.'
		return
	endif

endif

; go through amr structure, filling appropriate values into
; image box
boxmin = [amr.boxmin[dim1], amr.boxmin[dim2]]
boxmax = [amr.boxmax[dim1], amr.boxmax[dim2]]
for n=maxlevel, minlevel, -1 do begin

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

		; apply the threshhold
		if keyword_set(rhocut) then data = data * (data ge rhocut)

		; now sum in the appropriate direction
		column = reform(amr.gridspacing[dim3,n] * $
				total(data, dim3 + 1, /double), $
			        n_elements(xidxlist), n_elements(yidxlist))

		; figure out the index range on the level of the image that
		; the parts of the column image we have left correspond to
		colidxminref = [xidxlist[0], yidxlist[0]]
		colidxmaxref = [xidxlist[n_elements(xidxlist)-1], $
				yidxlist[n_elements(yidxlist)-1]] + 1
		for l=maxlevel-1,n,-1 do begin
		   colidxminref = colidxminref * amr.refratio[l]
		   colidxmaxref = colidxmaxref * amr.refratio[l]
		endfor
		colidxmaxref = colidxmaxref - 1

		; refine the column data to the level of the image
		colref = rebin(column, colidxmaxref[0]-colidxminref[0]+1, $
				colidxmaxref[1]-colidxminref[1]+1, /sample)

		; get index limits within the image box for the refined
		; column data
		imgidxlo = long(colidxminref-imglo)
		imgidxhi = long(colidxmaxref-imglo)

		; If the edge of the image box falls in between the cells
		; of a coarse fab, colidxminref may be smaller than imglo,
		; and colidxmaxref may be larger than imghi. If
		; that occurs, we will only use the column data that falls
		; within the image box.
		if imgidxlo[0] lt 0 then begin
			; clip at low x end
			colref = colref[imglo[0]-colidxminref[0]: $
					colidxmaxref[0]-colidxminref[0], *]
			imgidxlo[0] = 0
			colidxminref[0] = imglo[0]
		endif
		if imgidxlo[1] lt 0 then begin
			; clip at low y end
			colref = colref[*, imglo[1]-colidxminref[1]: $
					colidxmaxref[1]-colidxminref[1]]
			imgidxlo[1] = 0
			colidxminref[1] = imglo[1]
		endif
		if imgidxhi[0] gt imghi[0]-imglo[0] then begin
			; clip at high x end
			colref = colref[0:imghi[0]-colidxminref[0], *]
			imgidxhi[0] = imghi[0]-imglo[0]
		endif
		if imgidxhi[1] gt imghi[1]-imglo[1] then begin
			; clip at high y end
			colref = colref[*, 0:imghi[1]-colidxminref[1]]
			imgidxhi[1] = imghi[1]-imglo[1]
		endif

		; add the refined column data to the image
		img[imgidxlo[0]:imgidxhi[0], $
		    imgidxlo[1]:imgidxhi[1]] = $
		   img[imgidxlo[0]:imgidxhi[0], $
		       imgidxlo[1]:imgidxhi[1]] + colref

	endfor
endfor

; floor if requested
if keyword_set(floor) then begin
	imgnz = where(img ne 0.0)
	imgzer = where(img eq 0.0)
	if ((imgzer[0] ne -1) and (imgnz[0] ne -1)) then $
		img[imgzer] = min(img[imgnz])
endif

; smear the image with a beam if requested
if keyword_set(smearsize) then $
	img = beamsmear(img, amr.gridspacing[dim3,maxlevel], $
			smeardist, smearsize)

; take log of image if requested
if keyword_set(log) then begin
	; avoid taking log of zero or negative
	if (min(img) le 0.0) then begin
		if keyword_set(scalerange) then newmin=scalerange[0] $
		else newmin=min(img[where(img gt 0.0)])
		img[where(img le 0.0)] = newmin
	endif
	img = alog10(img)
endif

; store image
imgout=img

; exit if noplot is set
if keyword_set(noplot) then return

; if using nonscalable pixels, regrid image to correct size in pixels
if not scalablePixels then $
	img = congrid(img, imgdisphi[0]-imgdisplo[0]+1, $
			  imgdisphi[1]-imgdisplo[1]+1)

; set the range to be used to assign colors
if n_elements(scalerange) eq 0 then begin
	maxrange = max(img)
	minrange = min(img)
	if keyword_set(verbose) then begin
		print, 'Local minimum = ', minrange
		print, 'Local maximum = ', maxrange
	endif
endif else begin
	maxrange = scalerange[1]
	minrange = scalerange[0]
	if keyword_set(verbose) then begin
		print, 'User-provided minimum = ', minrange
		print, 'User-provided maximum = ', maxrange
	endif
endelse

; draw the color bar
if not keyword_set(nocolorbar) and not keyword_set(suppresscolorbar) $
	then begin

	; set up color bar string formatting
	if (abs(maxrange) lt 1.0e2) and (abs(maxrange) gt 1.0e-1) then $
		format = '(f9.2)' $
	else $
		format = '(e10.2)'

	; prevent pathological case of maxrange = minrange
	if maxrange eq minrange then begin
		if maxrange ne 0 then begin
			; give ourselves a +- 5% range
			maxrange = maxrange * 1.05
			minrange = minrange * 0.95
		endif else begin
			; just pick a number out of the air
			maxrange = 1.0
			minrange = -1.0
		endelse
	endif

	; set the colorbar position
	colorbarpos = fltarr(4)
	colorbarpos[0] = (plotcol + 0.88) / tempmulti[1]
	colorbarpos[1] = plotpos[1]
	if not keyword_set(skinnycolorbar) then $
		colorbarpos[2] = (plotcol+0.95) / tempmulti[1] $
	else $
		colorbarpos[2] = (plotcol+0.9105) / tempmulti[1]
	if not keyword_set(fitcolorbar) then colorbarpos[3] = plotpos[3] $
	else begin
		loplotpix=convert_coord(plotpos[0], plotpos[1], /normal, $
			/to_device)
		hiplotpix=convert_coord(plotpos[2], plotpos[3], /normal, $
			/to_device)
		axislength=hiplotpix-loplotpix
		axislength=min([axislength[0], axislength[1]])
		yhiplotpix=loplotpix[1]+axislength
		hiplotnorm=convert_coord(0, yhiplotpix, /device, /to_normal)
		colorbarpos[3] = hiplotnorm[1]
	endelse

	; draw the color bar
	if n_elements(ncolorbardiv) eq 0 then ncolorbardiv=6
	colorbar, /vertical, maxrange = maxrange, minrange = minrange, $
		format = format, ncolors=ncolors, position=colorbarpos, $
		charsize = charsize, divisions=ncolorbardiv
endif

; put pixels into specified range
if n_elements(scalerange) ne 0 then begin
	img = (img lt minrange) * minrange + (img ge minrange) * img
	img = (img gt maxrange) * maxrange + (img le maxrange) * img
endif

; scale image to color table
; code to deal with number of image colors lifted from Fanning colorbar.pro
;if scalablePixels then begin
;	oldDevice = !D.NAME

	; What kind of computer are we using? SET_PLOT to appropriate
	; display device.
;	thisOS = !VERSION.OS_FAMILY;
;	thisOS = STRMID(thisOS, 0, 3)
;	thisOS = STRUPCASE(thisOS)
;	CASE thisOS of
;		'MAC': SET_PLOT, thisOS
;		'WIN': SET_PLOT, thisOS
;		ELSE: SET_PLOT, 'X'
;	ENDCASE

	; Here is how many colors we should use.
	ncolors = !D.TABLE_SIZE
;	SET_PLOT, oldDevice

;ENDIF ELSE ncolors = !D.TABLE_SIZE

; scale the image to the color table
img1 = fix((img - minrange) / (maxrange - minrange) * (ncolors-1))

; Turn off device decomposition for scalable pixels
if not scalablePixels then begin
	if !Version.Release ge 5.2 then device, get_decomposed=thisdecomposed
	device, decomposed=0
endif

; Now display the image. How we do this again depends on whether we
; have scalable pixels or not.

if scalablePixels then begin

	; we have scalable pixels

	; get box limits in normal coordinates
	imglonorm = convert_coord(xlim[0], ylim[0], /data, /to_norm)
	imghinorm = convert_coord(xlim[1], ylim[1], /data, /to_norm)

	; show image
	tv, img1, imglonorm[0], imglonorm[1], $
		xsize = imghinorm[0] - imglonorm[0], $
		ysize = imghinorm[1] - imglonorm[1], $
		/norm

endif else begin

	; regrid image to correct size in pixels
	imgdisp = congrid(img1, imgdisphi[0]-imgdisplo[0]+1, $
			  imgdisphi[1]-imgdisplo[1]+1)

	; display
	tv, imgdisp, imgdisplo[0], imgdisplo[1]
endelse

; turn decomposition back to previous setting
if not scalablePixels then begin
	if !Version.Release ge 5.2 then device, decomposed=thisdecomposed
endif

if keyword_set(verbose) then begin
	print, 'imglo0 = ', imglo[0] / !D.X_SIZE
	print, 'imglo1 = ', imglo[1] / !D.Y_SIZE
	print, 'imghi0 = ', imghi[0] / !D.X_SIZE
	print, 'imglo1 = ', imghi[1] / !D.Y_SIZE
endif

; overlay the plot axes again
!p.multi=bangmulti
if keyword_set(nocolorbar) then $
	plot, [0,0], [0,0], /noerase, /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, charsize=charsize, _extra = extra $
else $
	plot, [0,0], [0,0], /noerase, /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra, $
		position = plotpos, charsize=charsize

; increment !p.multi, since the plot call with /noerase set doesn't do
; that automatically
!p.multi[0] = !p.multi[0] + 1
if !p.multi[0] eq tempmulti[1]*tempmulti[2]*tempmulti[3] then !p.multi[0]=0

return
end
