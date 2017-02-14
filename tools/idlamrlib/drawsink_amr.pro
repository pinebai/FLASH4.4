pro drawsink_amr, amr, dir, box=box, circle=circle, psym=psym, color=color, $
	vx=vx, vy=vy, length=length

; This routine draws the sink particles associated with an AMR object
; on the current plot, which is assumed to be projected in the direction
; dir. If the keyword box is set, it also draws an accretion
; box with a radius of box cells around the object. The symbol for the
; sink particle defaults to an asterisk (psym=2), but can be overriden by
; setting the psym keyword.

if not keyword_set(psym) then psym=2
if n_elements(color) eq 0 then color=!p.color
if not keyword_set(length) then length=1.0
if dir eq 0 then begin
	dir1 = 1
	dir2 = 2
endif else if dir eq 1 then begin
	dir1 = 0
	dir2 = 2
endif else if dir eq 2 then begin
	dir1 = 0
	dir2 = 1
endif else begin
	print, 'Error: dir must be 0, 1, or 2'
	return
endelse

if keyword_set(vx) then begin
	vxmax = max_amr(vx, vxmin)
	vymax = max_amr(vy, vymin)
	vmax = max([abs(vxmax), abs(vxmin), abs(vymax), abs(vymin)])
endif

maxlevel=amr.maxlevel
dx=amr.gridspacing[*,maxlevel]

for i=0, amr.nsink+amr.nstar-1 do begin

	if (i le amr.nsink-1) then begin
		pos = [amr.sinkparticles[i].x[dir1], $
		       amr.sinkparticles[i].x[dir2]]
		pos3d = amr.sinkparticles[i].x
		v = [amr.sinkparticles[i].p[dir1], $
		       amr.sinkparticles[i].p[dir2]]/amr.sinkparticles[i].m
	endif else begin
		pos = [amr.starparticles[i-amr.nsink].x[dir1], $
		       amr.starparticles[i-amr.nsink].x[dir2]]
		pos3d = amr.starparticles[i].x
		v = [amr.starparticles[i].p[dir1], $
		       amr.starparticles[i].p[dir2]]/amr.starparticles[i].m
	endelse

	plots, pos[0], pos[1], color=color, psym=psym

	if keyword_set(box) then begin
		loc=coord_to_fab(pos3d, amr)
		loc1=loc[5+dir1,maxlevel]
		loc2=loc[5+dir2,maxlevel]
		plots, amr.boxmin[dir1] + $
			dx[dir1]*[loc1-box,loc1-box,loc1+box+1,loc1+box+1,loc1-box], $
		       amr.boxmin[dir2] + $
			dx[dir2]*[loc2-box,loc2+box+1,loc2+box+1,loc2-box,loc2-box], $
		       color=color
	endif

	if keyword_set(circle) then begin
		loc=coord_to_fab(pos3d, amr)
		loc1=loc[5+dir1,maxlevel]
		loc2=loc[5+dir2,maxlevel]
		loclow=[amr.boxmin[dir1],amr.boxmin[dir2]] + $
			[dx[dir1],dx[dir2]]*[loc1,loc2]
		lochi=[amr.boxmin[dir1],amr.boxmin[dir2]] + $
			[dx[dir1],dx[dir2]]*[loc1+1,loc2+1]
		tmp=intarr(2*circle+1,2*circle+1)
		for n=-circle,circle do for m=-circle,circle do $
			if (n^2+m^2 le circle^2) then tmp[n+circle,m+circle]=1

		for n=-circle,circle-1 do begin

			; draw the horizontal segments
			row1=tmp[*,n+circle]
			row2=tmp[*,n+1+circle]
			segment=row1 xor row2
			for m=-circle,circle do begin
				if segment[m+circle] eq 1 then begin
				   plots, [loclow[0],lochi[0]]+dx[dir1]*m, $
				      [lochi[1],lochi[1]]+dx[dir2]*n, color=color
				endif
			endfor

			; draw the vertical segments
			col1=tmp[n+circle,*]
			col2=tmp[n+1+circle,*]
			segment=col1 xor col2
			for m=-circle,circle do begin
				if segment[m+circle] eq 1 then begin
				   plots, [lochi[0],lochi[0]]+dx[dir1]*n, $
				      [loclow[1],lochi[1]]+dx[dir2]*m, color=color
				endif
			endfor
		endfor

		; outer boundaries
		segment=tmp[*,0]
		for m=-circle,circle do begin
			if segment[m+circle] eq 1 then begin
			   plots, [loclow[0],lochi[0]]+dx[dir1]*m, $
			      [loclow[1],loclow[1]]+dx[dir2]*(-circle), color=color
			endif
		endfor
		segment=tmp[*,2*circle]
		for m=-circle,circle do begin
			if segment[m+circle] eq 1 then begin
			   plots, [loclow[0],lochi[0]]+dx[dir1]*m, $
			      [lochi[1],lochi[1]]+dx[dir2]*circle, color=color
			endif
		endfor
		segment=tmp[0,*]
		for m=-circle,circle do begin
			if segment[m+circle] eq 1 then begin
			   plots, [loclow[0],loclow[0]]+dx[dir1]*(-circle), $
			      [loclow[1],lochi[1]]+dx[dir2]*m, color=color
			endif
		endfor
		segment=tmp[2*circle,*]
		for m=-circle,circle do begin
			if segment[m+circle] eq 1 then begin
			   plots, [lochi[0],lochi[0]]+dx[dir1]*circle, $
			      [loclow[1],lochi[1]]+dx[dir2]*m, color=color
			endif
		endfor

	endif

	; draw velocity arrows
	if keyword_set(vx) then begin

		; figure out the length of the arrow
		arrowlength = length * sqrt(total(v^2))/vmax * $
			max(vx.gridspacing[[dir1,dir2],maxlevel])

		; draw arrow
		arrow, pos[0], pos[1], $
			pos[0] + arrowlength*v[0]/sqrt(total(v^2)), $
			pos[1] + arrowlength*v[1]/sqrt(total(v^2)), $
			/data, color=color, hsize=-0.3

	endif

endfor

end

