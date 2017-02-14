pro raster_amr, amr, plane, slice, anisotropic=anisotropic, $
	xrange=xrange, yrange=yrange, nocolorbar=nocolorbar, $
	verbose=verbose, zrange=zrange, log=log, $
	maxlevel=maxlevel, localrange=localrange, $
	snaptogrid=snaptogrid, charsize=charsize, $
	fitcolorbar=fitcolorbar, xmargin=xmargin, ymargin=ymargin, $
	imgout=imgout, $
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


; this routine produces a raster plot of an amr object. The image is
; sliced along the plane plane, so, for example, plane = 0 produces a
; yz slice, plane = 1 produces xz, and plane = 2 produces xy. Slice 
; specifies the physical location of the slice in data units; in other
; words, setting plane = 2, slice = 1.0e17 shows the xy plane sliced
; at z = 1.0e17. Both plane and slice are ignored for 2d amr objects.
; Effects of keywords (other than standard IDL graphics keywords):
;   anisotropic: plots use isotropic axes by default. If anisotropic is
;	set, axes are allowed to be anisotropic and the plot will fill
;	the available area
;   nocolorbar: supresses the production of a colorbar
;   verbose: make the routine print out some data about the amr object
;   maxlevel: causes the routine to only raster up to the specified
;	level of refinement. This can be useful in preventing the
;	production of overly large postscript files.
;   zrange: specifies the minimum and maximum for the range of the
;	color bar
;   localrange: if set, this specifies that the min and max of the color
;	bar are to be taken from the local data, not the global data
;   snaptogrid: if set to n, this keyword adjusts the x and y range so
;	that they fall exactly at cell boundaries on level n

; get keywords
if not keyword_set(anisotropic) then isotropic=1 else isotropic=0
if n_elements(maxlevel) eq 0 then maxlevel=amr.maxlevel
if not keyword_set(charsize) then charsize=!p.charsize
if not keyword_set(xmargin) then xmargin=!x.margin
if not keyword_set(ymargin) then ymargin=!y.margin

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

; set range of plot range
if not keyword_set(xrange) then xlim=[amr.boxmin[dim1], amr.boxmax[dim1]] $
else xlim=xrange
if not keyword_set(yrange) then ylim=[amr.boxmin[dim2], amr.boxmax[dim2]] $
else ylim=yrange

; if asked to do so, snap range to grid
if keyword_set(xrange) and (n_elements(snaptogrid) ne 0) then begin
	xlim[0] = floor( (xlim[0] - amr.boxmin[dim1]) / $
			 amr.gridspacing[dim1,n] ) * $
		  amr.gridspacing[dim1, n] + $
		  amr.boxmin[dim1]
	xlim[1] = ceil( (xlim[1] - amr.boxmin[dim1]) / $
			 amr.gridspacing[dim1,n] ) * $
		  amr.gridspacing[dim1, n] + $
		  amr.boxmin[dim1]
	if keyword_set(verbose) then print, 'New x range = ', xlim
endif
if keyword_set(yrange) and (n_elements(snaptogrid) ne 0) then begin
	ylim[0] = floor( (ylim[0] - amr.boxmin[dim2]) / $
			 amr.gridspacing[dim2,n] ) * $
		  amr.gridspacing[dim2, n] + $
		  amr.boxmin[dim2]
	ylim[1] = ceil( (ylim[1] - amr.boxmin[dim2]) / $
			 amr.gridspacing[dim2,n] ) * $
		  amr.gridspacing[dim2, n] + $
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

; see how many colors are available
ncolors = !d.table_size

; draw axes
if keyword_set(nocolorbar) then $
	plot, [0,0], [0,0], /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, charsize=charsize, $
		xmargin=xmargin, ymargin=ymargin, _extra = extra $
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

; What we do here depends on whether we have a device with scalable
; pixels or not. If we have scalable pixels, we can just construct
; an image box where each cell of the amr grid on the finest level
; represents a single pixel. We then rely on the scaling of pixels
; to make that image correctly overlay the plot window. If we don't
; have scalable pixels, the number of cells per pixel will in be a
; complicated function of the size of the graphics device and the
; resolution of the amr grid. For reference, X windows does not have
; scalable pixels, postscript does.

scalablePixels = !d.flags mod 2

if scalablePixels then begin

	; we have scalable pixels, so set image index limits equal to
	; index limits on fines

	; now set image size equal to number of cells on maxlevel
	; covered by the requested range
	imglo = [floor((xlim[0]-amr.boxmin[dim1]) / $
			amr.gridspacing[dim1,maxlevel]), $
		 floor((ylim[0]-amr.boxmin[dim2]) / $
			amr.gridspacing[dim2,maxlevel])]
	imghi = [ceil((xlim[1]-amr.boxmin[dim1]) / $
			amr.gridspacing[dim1,maxlevel]), $
		 ceil((ylim[1]-amr.boxmin[dim2]) / $
			amr.gridspacing[dim2,maxlevel])] - 1

endif else begin

	; we do not have scalable pixels

	; get max and min for image box in pixels
	imglo = convert_coord(xlim[0], ylim[0], /data, /to_device)
	imghi = convert_coord(xlim[1], ylim[1], /data, /to_device)
	imglo=round(imglo[0:1])
	imghi=round(imghi[0:1])
	if (imghi[0] - imglo[0]) lt 2 then begin
		print, 'Error: requested x axis range is too small.'
		return
	endif
	if (imghi[1] - imglo[1]) lt 2 then begin
		print, 'Error: requested y axis range is too small.'
		return
	endif

endelse

; construct the image box
img=fltarr(imghi[0]-imglo[0]+1, imghi[1]-imglo[1]+1)

; go through amr structure, filling appropriate values into
; image box
boxmin = [amr.boxmin[dim1], amr.boxmin[dim2]]
for n=0, maxlevel do begin

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

	; if there is a 3rd coordinate, figure out the index corresponding
	; to it on this level
	if amr.ndim eq 3 then zidx = floor((slice - amr.boxmin[dim3]) / $
		amr.gridspacing[dim3,n])

	for m=0, amr.levels[n].nfab-1 do begin

		; figure out which region of this fab overlaps with
		; the image box

		; check the slice direction first. If the requested slice
		; doesn't cut through this fab, then go to the next one.
		if amr.ndim eq 3 then begin
			if (zidx lt (*amr.levels[n].fabptr)[m].idxlo[dim3]) $
			or (zidx gt (*amr.levels[n].fabptr)[m].idxhi[dim3]) $
			then continue else $
			idx3 = zidx - (*amr.levels[n].fabptr)[m].idxlo[dim3]
		endif

		; set up some shorthands
		fabidxmin = [(*amr.levels[n].fabptr)[m].idxlo[dim1], $
			     (*amr.levels[n].fabptr)[m].idxlo[dim2]]
		fabidxmax = [(*amr.levels[n].fabptr)[m].idxhi[dim1], $
			     (*amr.levels[n].fabptr)[m].idxhi[dim2]]
		
		; find the overlap between our target index range and the
		; index range stored in this fab
		overlapmin = lonarr(2)
		overlapmax = lonarr(2)
		overlapmin[0] = (fabidxmin[0] gt xidxlim[0]) * fabidxmin[0] + $
			(fabidxmin[0] le xidxlim[0]) * xidxlim[0]
		overlapmin[1] = (fabidxmin[1] gt yidxlim[0]) * fabidxmin[1] + $
			(fabidxmin[1] le yidxlim[0]) * yidxlim[0]
		overlapmax[0] = (fabidxmax[0] le xidxlim[1]) * fabidxmax[0] + $
			(fabidxmax[0] gt xidxlim[1]) * xidxlim[1]
		overlapmax[1] = (fabidxmax[1] le yidxlim[1]) * fabidxmax[1] + $
			(fabidxmax[1] gt yidxlim[1]) * yidxlim[1]

		; if there is no overlap between this fab and the
		; image box, then move to the next fab
		if (overlapmin[0] gt overlapmax[0]) or $
		   (overlapmin[1] gt overlapmax[1]) then continue

		; there is an overlap, so set the limits on the portion of
		; the fab we wish to extract
		idxlo = overlapmin - fabidxmin
		idxhi = overlapmax - fabidxmin

		; We must convert limits in data coordinates to indices in
		; the image array. How we do this depends on whether we
		; have scalable pixels or not.

		if scalablePixels then begin

			; we have scalable pixels

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

		endif else begin

			; we do not have scalable pixels

			; get the physical position of the min and max
			; of the overlap box
			datamin = boxmin + (overlapmin) * dx
			datamax = boxmin + (overlapmax + 1) * dx
			imgdatalo = convert_coord(datamin[0], datamin[1], $
				/data, /to_device)
			imgdatahi = convert_coord(datamax[0], datamax[1], $
				/data, /to_device)
			imgidxlo = round(imgdatalo-imglo)
			imgidxhi = round(imgdatahi-imglo)

			; crop the image box index range to fit what we
			; have available
			imgidxlo = (imgidxlo gt 0) * imgidxlo
			imgidxhi = (imgidxhi le (imghi-imglo)) * imgidxhi + $
				   (imgidxhi gt (imghi-imglo)) * $
					(imghi-imglo)


		endelse

		; put data into image array
		if amr.ndim eq 2 then $
		; 2d case
		img[imgidxlo[0]:imgidxhi[0], imgidxlo[1]:imgidxhi[1]] = $
			congrid( (*(*amr.levels[n].fabptr)[m].dataptr)$
				     [idxlo[0]:idxhi[0], idxlo[1]:idxhi[1]], $
				  imgidxhi[0]-imgidxlo[0]+1, $
				  imgidxhi[1]-imgidxlo[1]+1) $
		else begin

		; 3 possible 3d cases
		if plane eq 0 then $
		img[imgidxlo[0]:imgidxhi[0], $
		    imgidxlo[1]:imgidxhi[1]] = $
			congrid( reform( $
				    (*(*amr.levels[n].fabptr)[m].dataptr)$
				     [idx3, idxlo[0]:idxhi[0], $
				      idxlo[1]:idxhi[1]], $
				    idxhi[0]-idxlo[0]+1, $
				    idxhi[1]-idxlo[1]+1), $
				  imgidxhi[0]-imgidxlo[0]+1, $
				  imgidxhi[1]-imgidxlo[1]+1) $
		else if plane eq 1 then $
		img[imgidxlo[0]:imgidxhi[0], $
		    imgidxlo[1]:imgidxhi[1]] = $
			congrid( reform( $
				    (*(*amr.levels[n].fabptr)[m].dataptr)$
				     [idxlo[0]:idxhi[0], idx3, $
				      idxlo[1]:idxhi[1]], $
				    idxhi[0]-idxlo[0]+1, $
				    idxhi[1]-idxlo[1]+1), $
				  imgidxhi[0]-imgidxlo[0]+1, $
				  imgidxhi[1]-imgidxlo[1]+1) $
		else $
		img[imgidxlo[0]:imgidxhi[0], $
		    imgidxlo[1]:imgidxhi[1]] = $
			congrid( reform( $
				    (*(*amr.levels[n].fabptr)[m].dataptr)$
				     [idxlo[0]:idxhi[0], $
				      idxlo[1]:idxhi[1], idx3], $
				    idxhi[0]-idxlo[0]+1, $
				    idxhi[1]-idxlo[1]+1), $
				  imgidxhi[0]-imgidxlo[0]+1, $
				  imgidxhi[1]-imgidxlo[1]+1)

		endelse

	endfor
endfor

; if requested, take the log of the image
if keyword_set(log) then img=alog10(img)

; set output img
imgout=img

; set the range to be used to assign colors
if n_elements(zrange) eq 0 then begin
	if not keyword_set(localrange) then begin
		maxrange = max_amr (amr, minrange)
		if keyword_set(log) then begin
			maxrange=alog10(maxrange)
			minrange=alog10(minrange)
		endif
		if keyword_set(verbose) then begin
			print, 'Global minimum = ', minrange
			print, 'Global maximum = ', maxrange
		endif
	endif else begin
		maxrange = max(img)
		minrange = min(img)
		if keyword_set(verbose) then begin
			print, 'Local minimum = ', minrange
			print, 'Local maximum = ', maxrange
		endif
	endelse
endif else begin
	maxrange = zrange[1]
	minrange = zrange[0]
	if keyword_set(verbose) then begin
		print, 'User-provided minimum = ', minrange
		print, 'User-provided maximum = ', maxrange
	endif
endelse

; draw the color bar
if not keyword_set(nocolorbar) then begin

	; set up color bar string formatting
	if (abs(maxrange) lt 1.0e4) and (abs(maxrange) gt 1.0e-4) then $
		format = '(f8.2)' $
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
	colorbarpos[2] = (plotcol+0.95) / tempmulti[1]
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
	colorbar, /vertical, maxrange = maxrange, minrange = minrange, $
		format = format, ncolors=ncolors, position=colorbarpos, $
		charsize = charsize
endif

; put pixels into specified range
if n_elements(zrange) ne 0 then begin
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

	; we do not have scalable pixels
	tv, img1, imglo[0], imglo[1]
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
		isotropic=isotropic, _extra = extra, charsize=charsize, $
		xmargin=xmargin, ymargin=ymargin $
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
