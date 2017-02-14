function courant
; return the block number of where the courant condition is determined

; common in the variables from the read-in routine
common size, tot_blocks, nvar, nxb, nyb, nzb, ntopx, ntopy, ntopz
common time, time, dt
common tree, lrefine, nodetype, gid, coord, size, bnd_box
common vars, unk, unk_names


; if gamc is not defined, make one up -- this usually means it is
; constant, so the zone that sets the courant condition won't be
; affected
if n_elements(gamc) EQ 0 then gamc = 1.4

; compute the sound speed
c_s = reform(sqrt(unk[var_index('gamc'),*,*,*,*] * $
                  unk[var_index('pres'),*,*,*,*] / $
                  unk[var_index('dens'),*,*,*,*]))

dx = fltarr(nxb,nyb)
dy = fltarr(nxb,nyb)

cfl = fltarr(tot_blocks)

; make the loop index long so we can deal with more than 32768 blocks
i = (lonarr(1))[0]

; there should be an easier IDLish way to do this w/o loops -- but I'm
; too lazy right now
for i = 0l, tot_blocks-1l do begin
    dx[*,*] = size[0,i] * reform(replicate(1.,nxb*nyb),nxb,nyb)/float(nxb)
    dy[*,*] = size[1,i] * reform(replicate(1.,nxb*nyb),nxb,nyb)/float(nyb)

    cfl[i] = min(1./(abs(unk[var_index('velx'),i,*,*,*])/dx + $
                     abs(unk[var_index('vely'),i,*,*,*])/dy + $
                     c_s[i,*,*]*sqrt((1./dx)^2 + (1./dy)^2)))

endfor

min_cfl = (where(cfl EQ min(cfl) AND nodetype EQ 1))[0]

undefine, cfl
undefine, c_s

return, min_cfl
end









