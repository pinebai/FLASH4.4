function refine, data, d1, d2, d3

; This function does linear refinement as it should be done in amr, i.e.
; assuming that data points are located at cell centers and using
; linear extrapolation past the last coarse data point if necessary. The
; parameters d1, d2, and d3 are the dimensions of the result. They must
; be equal to or greater than the dimensions of the data in every
; dimension, and must be integer multiples of the dimensions of data.

; get dimensions of data we've been passed
sz=size(data)
ndim=sz[0]
if (ndim lt 1) or (ndim gt 3) then return, -1
l1=sz[1]
if ndim gt 1 then l2=sz[2]
if ndim gt 2 then l3=sz[3]

; check for errors
if d1 lt l1 then return, -1
if ndim gt 1 then if d2 lt l2 then return, -1
if ndim gt 2 then if d3 lt l3 then return, -1
if (d1/l1) ne (float(d1)/float(l1)) then return, -1
if ndim gt 1 then if (d2/l2) ne (float(d2)/float(l2)) then return, -1
if ndim gt 2 then if (d3/l3) ne (float(d3)/float(l3)) then return, -1

; get refinement ratios
f1=d1/l1
if ndim gt 1 then f2=d2/l2
if ndim gt 2 then f3=d3/l3

; compute x gradients
if ndim eq 1 then begin
	grad=dblarr(l1+1)
	grad[1:l1-1]=(data[1:l1-1]-data[0:l1-2])/f1
	grad[0]=grad[1]
	grad[l1]=grad[l1-1]
endif else if ndim eq 2 then begin
	grad=dblarr(l1+1,l2)
	grad[1:l1-1,*]=(data[1:l1-1,*]-data[0:l1-2,*])/f1
	grad[0,*]=grad[1,*]
	grad[l1,*]=grad[l1-1,*]
endif else begin
	grad=dblarr(l1+1,l2,l3)
	grad[1:l1-1,*,*]=(data[1:l1-1,*,*]-data[0:l1-2,*,*])/f1
	grad[0,*,*]=grad[1,*,*]
	grad[l1,*,*]=grad[l1-1,*,*]
endelse

; create refined data in the x direction
if ndim eq 1 then xrefdata=dblarr(d1)
if ndim eq 2 then xrefdata=dblarr(d1,l2)
if ndim eq 3 then xrefdata=dblarr(d1,l2,l3)

; do nothing if f1 eq 1
if f1 ne 1 then begin

   ; first deal with the points where we can interpolate
   for i=0,f1-1 do begin
	idxlist=ceil(f1/2.) + lindgen(l1-1)*f1 + i
	if ndim eq 1 then begin
	   if (f1/2) eq (f1/2.) then $
	      xrefdata[idxlist]= $
		data[0:l1-2]+grad[1:l1-1]*(i+0.5) $
	   else $
	      xrefdata[idxlist]= $
		data[0:l1-2]+grad[1:l1-1]*(i+1)
	endif else if ndim eq 2 then begin
	   if (f1/2) eq (f1/2.) then $
		xrefdata[idxlist,*]= $
	   	  data[0:l1-2,*]+grad[1:l1-1,*]*(i+0.5) $
	   else $
		xrefdata[idxlist,*]= $
	   	  data[0:l1-2,*]+grad[1:l1-1,*]*(i+1)
	endif else if ndim eq 3 then begin
	   if (f1/2) eq (f1/2.) then $
		xrefdata[idxlist,*,*]= $
		  data[0:l1-2,*,*]+grad[1:l1-1,*,*]*(i+0.5) $
	   else $
		xrefdata[idxlist,*,*]= $
		  data[0:l1-2,*,*]+grad[1:l1-1,*,*]*(i+1)
	endif
   endfor
   if (f1/2.) ne (f1/2) then begin
	idxlist=floor(f1/2.) + lindgen(l1)*f1
	if ndim eq 1 then xrefdata[idxlist] = data $
	else if ndim eq 2 then xrefdata[idxlist,*] = data $
	else xrefdata[idxlist,*,*] = data
   endif

   ; now deal with the endpoints, where we must extrapolate
   if (f1/2.) eq (f1/2) then begin
      if ndim eq 1 then begin
	xrefdata[0:ceil(f1/2.)-1]= $
		data[0]+grad[0]*(dindgen(ceil(f1/2.))+0.5-ceil(f1/2.))
	xrefdata[d1-ceil(f1/2.):d1-1]= $
		data[l1-1]+grad[l1]*(dindgen(ceil(f1/2.))+0.5)
      endif else if ndim eq 2 then begin
	for j=0,l2-1 do begin
		xrefdata[0:ceil(f1/2.)-1,j]= $
			data[0,j]+grad[0,j]*(dindgen(ceil(f1/2.))+0.5-f1/2)
		xrefdata[d1-ceil(f1/2.):d1-1,j]= $
			data[l1-1,j]+grad[l1,j]*(dindgen(ceil(f1/2.))+0.5)
	endfor
      endif else begin
	for j=0,l2-1 do begin
	   for k=0,l3-1 do begin
		xrefdata[0:ceil(f1/2.)-1,j,k]= $
			data[0,j,k]+grad[0,j,k]*(dindgen(ceil(f1/2.))+0.5-f1/2)
		xrefdata[d1-ceil(f1/2.):d1-1,j,k]= $
			data[l1-1,j,k]+grad[l1,j,k]*(dindgen(ceil(f1/2.))+0.5)
	   endfor
	endfor
      endelse

   endif else begin

      if ndim eq 1 then begin
	xrefdata[0:floor(f1/2.)-1]= $
		data[0]+grad[0]*(dindgen(floor(f1/2.))-floor(f1/2.))
	xrefdata[d1-floor(f1/2.):d1-1]= $
		data[l1-1]+grad[l1]*(dindgen(floor(f1/2.))+1)
      endif else if ndim eq 2 then begin
	for j=0,l2-1 do begin
		xrefdata[0:floor(f1/2.)-1,j]= $
			data[0,j]+grad[0,j]*(dindgen(floor(f1/2.))-f1/2)
		xrefdata[d1-floor(f1/2.):d1-1,j]= $
			data[l1-1,j]+grad[l1,j]*(dindgen(floor(f1/2.))+1)
	endfor
      endif else begin
	for j=0,l2-1 do begin
	   for k=0,l3-1 do begin
		xrefdata[0:floor(f1/2.)-1,j,k]= $
			data[0,j,k]+grad[0,j,k]*(dindgen(floor(f1/2.))-f1/2)
		xrefdata[d1-floor(f1/2.):d1-1,*,*]= $
			data[l1-1,j,k]+grad[l1,j,k]*(dindgen(floor(f1/2.))+1)
	   endfor
	endfor
      endelse
   endelse

endif else xrefdata = data

; now repeat the process in the y direction using the data that has
; already been refined in x, if we are 2 or 3d
if ndim eq 1 then return, xrefdata

; compute y gradients
if ndim eq 2 then begin
	grad=dblarr(d1,l2+1)
	grad[*,1:l2-1]=(xrefdata[*,1:l2-1]-xrefdata[*,0:l2-2])/f2
	grad[*,0]=grad[*,1]
	grad[*,l2]=grad[*,l2-1]
endif else begin
	grad=dblarr(d1,l2+1,l3)
	grad[*,1:l2-1,*]=(xrefdata[*,1:l2-1,*]-xrefdata[*,0:l2-2,*])/f2
	grad[*,0,*]=grad[*,1,*]
	grad[*,l2,*]=grad[*,l2-1,*]
endelse

; create refined data in the y direction
if ndim eq 2 then yrefdata=dblarr(d1,d2)
if ndim eq 3 then yrefdata=dblarr(d1,d2,l3)

; do nothing if f2 eq 1
if f2 ne 1 then begin

   ; first deal with the points where we can interpolate
   for i=0,f2-1 do begin
	idxlist=ceil(f2/2.) + lindgen(l2-1)*f2 + i
	if ndim eq 2 then begin
	   if (f2/2) eq (f2/2.) then $
		yrefdata[*,idxlist]= $
	   	  xrefdata[*,0:l2-2]+grad[*,1:l2-1]*(i+0.5) $
	   else $
		yrefdata[*,idxlist]= $
	   	  xrefdata[*,0:l2-2]+grad[*,1:l2-1]*(i+1)
	endif else if ndim eq 3 then begin
	   if (f2/2) eq (f2/2.) then $
		yrefdata[*,idxlist,*]= $
		  xrefdata[*,0:l2-2,*]+grad[*,1:l2-1,*]*(i+0.5) $
	   else $
		yrefdata[*,idxlist,*]= $
		  xrefdata[*,0:l2-2,*]+grad[*,1:l2-1,*]*(i+1)
	endif
   endfor
   if (f2/2.) ne (f2/2) then begin
	idxlist=floor(f2/2.) + lindgen(l2)*f2
	if ndim eq 2 then yrefdata[*,idxlist] = xrefdata $
	else yrefdata[*,idxlist,*] = xrefdata
   endif

   ; now deal with the endpoints, where we must extrapolate
   if (f2/2.) eq (f2/2) then begin
      if ndim eq 2 then begin
	for i=0,d1-1 do begin
		yrefdata[i,0:ceil(f2/2.)-1]= $
			xrefdata[i,0]+grad[i,0]*(dindgen(ceil(f2/2.))+0.5-f2/2)
		yrefdata[i,d2-ceil(f2/2.):d2-1]= $
			xrefdata[i,l2-1]+grad[i,l2]*(dindgen(ceil(f2/2.))+0.5)
	endfor
      endif else begin
	for i=0,d1-1 do begin
	   for k=0,l3-1 do begin
		yrefdata[i,0:ceil(f2/2.)-1,k]= $
		   xrefdata[i,0,k]+grad[i,0,k]*(dindgen(ceil(f2/2.))+0.5-f2/2)
		yrefdata[i,d2-ceil(f2/2.):d2-1,k]= $
		   xrefdata[i,l2-1,k]+grad[i,l2,k]*(dindgen(ceil(f2/2.))+0.5)
	   endfor
	endfor
      endelse

   endif else begin

      if ndim eq 2 then begin
	for i=0,d1-1 do begin
		yrefdata[i,0:floor(f2/2.)-1]= $
			xrefdata[i,0]+grad[i,0]*(dindgen(floor(f2/2.))-f2/2)
		yrefdata[i,d2-floor(f2/2.):d2-1]= $
			xrefdata[i,l2-1]+grad[i,l2]*(dindgen(floor(f2/2.))+1)
	endfor
      endif else begin
	for i=0,d1-1 do begin
	   for k=0,l3-1 do begin
		yrefdata[i,0:floor(f2/2.)-1,k]= $
		   xrefdata[i,0,k]+grad[i,0,k]*(dindgen(floor(f2/2.))-f2/2)
		yrefdata[i,d2-floor(f2/2.):d2-1,k]= $
		   xrefdata[i,l2-1,k]+grad[i,l2,k]*(dindgen(floor(f2/2.))+1)
	   endfor
	endfor
      endelse
   endelse

endif else yrefdata = xrefdata

if ndim eq 2 then return, yrefdata

; repeat for z direction

; compute z gradients
grad=dblarr(d1,d2,l3+1)
grad[*,*,1:l3-1]=(yrefdata[*,*,1:l3-1]-yrefdata[*,*,0:l3-2])/f3
grad[*,*,0]=grad[*,*,1]
grad[*,*,l3]=grad[*,*,l3-1]

; create refined data in the z direction
zrefdata=dblarr(d1,d2,d3)

; do nothing if f3 eq 1
if f3 ne 1 then begin

   ; first deal with the points where we can interpolate
   for i=0,f3-1 do begin
	idxlist=ceil(f3/2.) + lindgen(l3-1)*f3 + i
	if (f3/2) eq (f3/2.) then $
		zrefdata[*,*,idxlist]= $
		  yrefdata[*,*,0:l3-2]+grad[*,*,1:l3-1]*(i+0.5) $
	   else $
		zrefdata[*,*,idxlist]= $
		  yrefdata[*,*,0:l3-2]+grad[*,*,1:l3-1]*(i+1)
   endfor
   if (f3/2.) ne (f3/2) then begin
	idxlist=floor(f3/2.) + lindgen(l3)*f3
	zrefdata[*,*,idxlist] = yrefdata
   endif

   ; now deal with the endpoints, where we must extrapolate
   if (f3/2.) eq (f3/2) then begin
	for i=0,d1-1 do begin
	   for j=0,d2-1 do begin
		zrefdata[i,j,0:ceil(f3/2.)-1]= $
		   yrefdata[i,j,0]+grad[i,j,0]*(dindgen(ceil(f3/2.))+0.5-f3/2)
		zrefdata[i,j,d3-ceil(f3/2.):d3-1]= $
		   yrefdata[i,j,l3-1]+grad[i,j,l3]*(dindgen(ceil(f3/2.))+0.5)
	   endfor
	endfor

   endif else begin

	for i=0,d1-1 do begin
	   for j=0,d2-1 do begin
		zrefdata[i,j,0:floor(f3/2.)-1]= $
		   yrefdata[i,j,0]+grad[i,j,0]*(dindgen(floor(f3/2.))-f3/2)
		zrefdata[i,j,d3-floor(f3/2.):d3y-1]= $
		   yrefdata[i,j,l3-1]+grad[i,j,l3]*(dindgen(floor(f3/2.))+1)
	   endfor
	endfor
   endelse

endif else zrefdata = yrefdata

return, zrefdata

end




 