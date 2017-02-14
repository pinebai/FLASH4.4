function radial_average_polar, filename, variable, RCOORD = rc, TIME = time

; radial average 2-d spherical data


read_amr, filename, VAR_NAME=variable, TREE=tree, DATA=var, PARAMETERS=params

time = params.time

rmax = max(tree[*].bndBox[1,0])
rmin = min(tree[*].bndBox[0,0])

sample = 1.

var = reform(temporary(var))

r = fltarr(params.totBlocks,params.nxb)

; compute the radial distance of each point
for i = 0l, params.totBlocks-1 do begin

    rcoords = (tree[i].bndBox[1,0] - tree[i].bndBox[0,0])* $
      ((findgen(params.nxb) + .5)/float(params.nxb)) + tree[i].bndBox[0,0]

    r[i,*] = rcoords

endfor

index_good = where(tree[*].nodeType EQ 1)

svar = var[index_good,*,*]
sr   = r[index_good,*]

; averge the theta out of each block, since they are at the same
; radius
vravg = total(var,3)/params.nyb

; create our radial bin array
npts = params.ntopx*params.nxb*2^(max(tree[*].lrefine)-1)

rl = findgen(npts)*(rmax - rmin)/npts + rmin
rr = (findgen(npts) + 1.0)*(rmax - rmin)/npts + rmin

rc = 0.5*(rl + rr)

count = lonarr(npts)

var_avg = fltarr(npts)

count[*] = 0
var_avg[*] = 0

print, 'made it here'

; this is a really crappy way to do this
for n = 0, params.totBlocks-1 do begin

    for i = 0, params.nxb-1 do begin

        index = where(r[n,i] GE rl AND r[n,i] LT rr)

        count[index] = count[index] + 1
        var_avg[index] = var_avg[index] + vravg[n,i]

    endfor

endfor

print, 'made it here too'

var_avg = var_avg/(count > 1)

return, var_avg

end
        







