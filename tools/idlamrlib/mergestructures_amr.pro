pro mergestructures_amr, amr1, amr2, outamr1, outamr2, maxlevel=maxlevel

; Copyright Mark Krumholz and Robert Fisher (2002)
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


; This routine expects two inputs, amr1 and amr2, which are amr objects
; covering the same physical domain. It returns two amr objects,
; outamr1 and outamr2. These variables should be uninitialized when
; passed into the mergestructures_amr.
;
; This routine takes as input two amr objects that cover the same
; physical domain but may have different grid structures. It creates
; a new structure such that, at every physical point in space, the
; finest level covering that point is equal to the min of the finest
; levels of for that point in each of the two input objects, amr1 and
; amr2. The routine copies the appropriate data from each of the input
; amr objects into the new structure.


; Check for keyword
if (n_elements(maxlevel) eq 0) then maxlevel=amr1.maxlevel

; First check that the two objects cover the same domain and have
; the same gridspacing and periodicity on every level.
ndim=amr1.ndim
if (ndim ne amr2.ndim) then begin
  print, 'ERROR : Dimensions not equal.'
  return
endif

lowcheck = amr1.boxmin eq amr2.boxmin
highcheck = amr1.boxmax eq amr2.boxmax

if (total(lowcheck+highcheck) ne 2*ndim) then begin
  print, 'ERROR : Domain sizes not equal.'
  return
endif

outmaxlevel = min([amr1.maxlevel, amr2.maxlevel, maxlevel])
for n=0,outmaxlevel do begin
	if (total(amr1.gridspacing[*,n] eq amr2.gridspacing[*,n]) ne ndim) $
		then begin
                  print, "ERROR : Grid spacings not equal on every level."
                  return
        endif

;	if (total(amr1.periodicity[*,n] eq amr2.periodicity[*,n]) ne ndim) $
;		then begin
;                   print, 'ERROR : Periodicity not same on every level.'
;                   print, n
;                   return
;        endif
endfor

; Create new amr objects to hold the results, and set the number
; of levels in each object to be equal to min of the number of levels
; in either of the two source objects.
outamr1 = amr1
outamr2 = amr2
if (amr1.maxlevel gt outmaxlevel) then begin
	outamr1.maxlevel = outmaxlevel
	outamr1.idxlo = outamr1.idxlo[*,0:outmaxlevel]
	outamr1.idxhi = outamr1.idxhi[*,0:outmaxlevel]
;	outamr1.periodicity = outamr1.periodicity[*,0:outmaxlevel]
	outamr1.gridspacing = outamr1.gridspacing[*,0:outmaxlevel]
	outamr1.levels = outamr1.levels[0:outmaxlevel]
endif
if (amr2.maxlevel gt outmaxlevel) then begin
	outamr2.maxlevel = outmaxlevel
	outamr2.idxlo = outamr2.idxlo[*,0:outmaxlevel]
	outamr2.idxhi = outamr2.idxhi[*,0:outmaxlevel]
;	outamr2.periodicity = outamr2.periodicity[*,0:outmaxlevel]
	outamr2.gridspacing = outamr2.gridspacing[*,0:outmaxlevel]
	outamr2.levels = outamr2.levels[0:outmaxlevel]
endif

; Clean out some of the level information in each of the new amr objects.
for n=0, outmaxlevel do begin
	outamr1.levels[n].nfab = 0
	outamr1.levels[n].fabptr = ptr_new()
endfor

; Now start looping through levels and fabs of amr1. The strategy for
; each fab is to find all its intersection with all other fabs on that
; level in amr2. Each intersection will generate a new fab in outamr.
; We make a first pass to find all the intersections, since we need to
; know how many there are so we can allocate memory for them. Then
; we go back through and actually copy the data.

for n=0, outmaxlevel do begin

   print, 'Working on level = ', n

   ; Create an intersection list
   nintersect = 0
   intersectlist = lonarr(long(amr1.levels[n].nfab) * $
	long(amr2.levels[n].nfab),2)

   ; Loop through fabs to find intersections
   for m1=0, amr1.levels[n].nfab-1 do begin
	for m2=0, amr2.levels[n].nfab-1 do begin

	   ; Check for intersection
	   check1 = (*amr1.levels[n].fabptr)[m1].idxlo gt $
		    (*amr2.levels[n].fabptr)[m2].idxhi
	   check2 = (*amr1.levels[n].fabptr)[m1].idxhi lt $
		    (*amr2.levels[n].fabptr)[m2].idxlo

	   ; If intersection is found, record it
	   if (total(check1 + check2) eq 0) then begin
		intersectlist[nintersect,0]=m1
		intersectlist[nintersect,1]=m2
		nintersect=nintersect+1
	   endif
	endfor
   endfor

   ; Allocate memory for correct number of fab pointers in output amr objects
   if (nintersect gt 0) then begin
	outamr1.levels[n].fabptr = $
	    ptr_new(replicate((*amr1.levels[0].fabptr)[0],nintersect))
	outamr2.levels[n].fabptr = $
	    ptr_new(replicate((*amr1.levels[0].fabptr)[0],nintersect))
   endif

   ; Set the correct number of fabs
   outamr1.levels[n].nfab = nintersect
   outamr2.levels[n].nfab = nintersect
   outamr1.levels[n].ngrids = nintersect
   outamr2.levels[n].ngrids = nintersect

   ; Now go through intersections and copy the data
   for m=0, nintersect-1 do begin

	; Shortcuts for fab numbers
	fab1=intersectlist[m,0]
	fab2=intersectlist[m,1]

	; Grab the index limits in the intersecting fabs
	idxlo1 = (*amr1.levels[n].fabptr)[fab1].idxlo
	idxlo2 = (*amr2.levels[n].fabptr)[fab2].idxlo
	idxhi1 = (*amr1.levels[n].fabptr)[fab1].idxhi
	idxhi2 = (*amr2.levels[n].fabptr)[fab2].idxhi

	; Figure out the index range of the intersection region
	intersectidxlo = (idxlo1 gt idxlo2)*idxlo1 + (idxlo1 le idxlo2)*idxlo2
	intersectidxhi = (idxhi1 lt idxhi2)*idxhi1 + (idxhi1 ge idxhi2)*idxhi2

	; Figure out the physical limits associated with these index limits
	intersectxlo = amr1.boxmin + amr1.gridspacing[*,n] * $
		(intersectidxlo - 0.5)
	intersectxhi = amr1.boxmin + amr1.gridspacing[*,n] * $
		(intersectidxhi + 0.5)

	; Figure out the corresponding index ranges in each source fab
	fabidxlo1 = intersectidxlo - idxlo1
	fabidxhi1 = intersectidxhi - idxlo1
	fabidxlo2 = intersectidxlo - idxlo2
	fabidxhi2 = intersectidxhi - idxlo2

	; Grab the intersection data from each input amr object
	if (ndim eq 1) then begin
	   data1 = (*(*amr1.levels[n].fabptr)[fab1].dataptr) $
			[fabidxlo1:fabidxhi1]
	   data2 = (*(*amr2.levels[n].fabptr)[fab2].dataptr) $
			[fabidxlo2:fabidxhi2]
	endif else if (ndim eq 2) then begin
	   data1 = (*(*amr1.levels[n].fabptr)[fab1].dataptr) $
			[fabidxlo1[0]:fabidxhi1[0], $
			 fabidxlo1[1]:fabidxhi1[1]]
	   data2 = (*(*amr2.levels[n].fabptr)[fab2].dataptr) $
			[fabidxlo2[0]:fabidxhi2[0], $
			 fabidxlo2[1]:fabidxhi2[1]]
	endif else begin
	   data1 = (*(*amr1.levels[n].fabptr)[fab1].dataptr) $
			[fabidxlo1[0]:fabidxhi1[0], $
			 fabidxlo1[1]:fabidxhi1[1], $
			 fabidxlo1[2]:fabidxhi1[2]]
	   data2 = (*(*amr2.levels[n].fabptr)[fab2].dataptr) $
			[fabidxlo2[0]:fabidxhi2[0], $
			 fabidxlo2[1]:fabidxhi2[1], $
			 fabidxlo2[2]:fabidxhi2[2]]
	endelse

	; Copy the data into the output amr objects
	(*outamr1.levels[n].fabptr)[m].idxlo = intersectidxlo
	(*outamr2.levels[n].fabptr)[m].idxlo = intersectidxlo
	(*outamr1.levels[n].fabptr)[m].idxhi = intersectidxhi
	(*outamr2.levels[n].fabptr)[m].idxhi = intersectidxhi
	(*outamr1.levels[n].fabptr)[m].idxtype = $
		(*amr1.levels[0].fabptr)[0].idxtype
	(*outamr2.levels[n].fabptr)[m].idxtype = $
		(*amr2.levels[0].fabptr)[0].idxtype
	(*outamr1.levels[n].fabptr)[m].xlo = intersectxlo
	(*outamr2.levels[n].fabptr)[m].xlo = intersectxlo
	(*outamr1.levels[n].fabptr)[m].xhi = intersectxhi
	(*outamr2.levels[n].fabptr)[m].xhi = intersectxhi
	(*outamr1.levels[n].fabptr)[m].dataptr = ptr_new(data1)
	(*outamr2.levels[n].fabptr)[m].dataptr = ptr_new(data2)

   endfor

endfor

; Find highest level that has any cells on it
for n=outmaxlevel,0,-1 do begin
    if (outamr1.levels[n].nfab ne 0) then begin
	newmaxlevel = n
	break
    endif
endfor

; Delete unneeded levels
if (newmaxlevel ne outmaxlevel) then begin
    if (outamr1.nsink ne 0) then begin
	newAmrDescriptor1 = { name:outamr1.name, $
	    componentName:outamr1.componentName, $
	    version:outamr1.version, $
	    ndim:outamr1.ndim, $
	    time:outamr1.time, $
	    maxlevel:newmaxlevel, $
	    boxmin:outamr1.boxmin, boxmax:outamr1.boxmax, $
	    refratio:outamr1.refratio[0:newmaxlevel], $
	    idxlo:outamr1.idxlo[*,0:newmaxlevel], $
	    idxhi:outamr1.idxhi[*,0:newmaxlevel], $
	    idxtype:outamr1.idxtype, $
;	    periodicity:outamr1.periodicity[*,0:newmaxlevel], $
	    levelsteps:outamr1.levelsteps[0:newmaxlevel], $
	    gridspacing:outamr1.gridspacing[*,0:newmaxlevel], $
	    coordtype:outamr1.coordtype, $
	    levels:outamr1.levels[0:newmaxlevel], $
	    nsink:outamr1.nsink, sinkparticles:outamr1.sinkparticles }
    endif else begin
	newAmrDescriptor1 = { name:outamr1.name, $
	    componentName:outamr1.componentName, $
	    version:outamr1.version, $
	    ndim:outamr1.ndim, $
	    time:outamr1.time, $
	    maxlevel:newmaxlevel, $
	    boxmin:outamr1.boxmin, boxmax:outamr1.boxmax, $
	    refratio:outamr1.refratio[0:newmaxlevel], $
	    idxlo:outamr1.idxlo[*,0:newmaxlevel], $
	    idxhi:outamr1.idxhi[*,0:newmaxlevel], $
	    idxtype:outamr1.idxtype, $
;	    periodicity:outamr1.periodicity[*,0:newmaxlevel], $
	    levelsteps:outamr1.levelsteps[0:newmaxlevel], $
	    gridspacing:outamr1.gridspacing[*,0:newmaxlevel], $
	    coordtype:outamr1.coordtype, $
	    levels:outamr1.levels[0:newmaxlevel], $
	    nsink:outamr1.nsink }
    endelse
    if (outamr2.nsink ne 0) then begin
	newAmrDescriptor2 = { name:outamr2.name, $
	    componentName:outamr2.componentName, $
	    version:outamr2.version, $
	    ndim:outamr2.ndim, $
	    time:outamr2.time, $
	    maxlevel:newmaxlevel, $
	    boxmin:outamr2.boxmin, boxmax:outamr2.boxmax, $
	    refratio:outamr2.refratio[0:newmaxlevel], $
	    idxlo:outamr2.idxlo[*,0:newmaxlevel], $
	    idxhi:outamr2.idxhi[*,0:newmaxlevel], $
	    idxtype:outamr2.idxtype, $
;	    periodicity:outamr2.periodicity[*,0:newmaxlevel], $
	    levelsteps:outamr2.levelsteps[0:newmaxlevel], $
	    gridspacing:outamr2.gridspacing[*,0:newmaxlevel], $
	    coordtype:outamr2.coordtype, $
	    levels:outamr2.levels[0:newmaxlevel], $
	    nsink:outamr2.nsink, sinkparticles:outamr2.sinkparticles }
    endif else begin
	newAmrDescriptor2 = { name:outamr2.name, $
	    componentName:outamr2.componentName, $
	    version:outamr2.version, $
	    ndim:outamr2.ndim, $
	    time:outamr2.time, $
	    maxlevel:newmaxlevel, $
	    boxmin:outamr2.boxmin, boxmax:outamr2.boxmax, $
	    refratio:outamr2.refratio[0:newmaxlevel], $
	    idxlo:outamr2.idxlo[*,0:newmaxlevel], $
	    idxhi:outamr2.idxhi[*,0:newmaxlevel], $
	    idxtype:outamr2.idxtype, $
;	    periodicity:outamr2.periodicity[*,0:newmaxlevel], $
	    levelsteps:outamr2.levelsteps[0:newmaxlevel], $
	    gridspacing:outamr2.gridspacing[*,0:newmaxlevel], $
	    coordtype:outamr2.coordtype, $
	    levels:outamr2.levels[0:newmaxlevel], $
	    nsink:outamr2.nsink }
    endelse
    outamr1=newAmrDescriptor1
    outamr2=newAmrDescriptor2
endif

end

