function radial_average, filename, variable, XCENTER = xctr, YCENTER = yctr, $
                         RCOORD = rBinCenters



read_amr, filename, VAR_NAME=variable, TREE=tree, DATA=var, PARAMETERS=params

xmax = max(tree[*].bndBox[1,0])
xmin = min(tree[*].bndBox[0,0])

ymax = max(tree[*].bndBox[1,1])
ymin = min(tree[*].bndBox[0,1])

print, xmin, xmax, ymin, ymax

sample = 1.

if n_elements(xctr) EQ 0 then xctr = .5*(xmax - xmin)
if n_elements(yctr) EQ 0 then yctr = .5*(ymax - ymin)

print, xctr, yctr

var = reform(temporary(var))

s = size(var)
print, s

r = fltarr(s[1],s[2],s[3])

; compute the radial distance of each point
for i = 0l, s[1]-1 do begin

    xcoords = (tree[i].bndBox[1,0] - tree[i].bndBox[0,0])* $
      ((findgen(params.nxb) + .5)/float(params.nxb)) + tree[i].bndBox[0,0]
    x = xcoords # replicate(1.,params.nyb)

    ycoords = (tree[i].bndBox[1,1] - tree[i].bndBox[0,1])* $
      ((findgen(params.nyb) + .5)/float(params.nyb)) + tree[i].bndBox[0,1]
    y = replicate(1.,params.nxb) # ycoords

    r[i,*,*] = sqrt((x - xctr)^2 + (y - yctr)^2)

endfor

index_good = where(tree[*].nodeType EQ 1)

svar = var[index_good,*,*]
sr   = r[index_good,*,*]

s = size(svar)

r1d = reform(sr,s[1]*s[2]*s[3])
var1d = reform(svar,s[1]*s[2]*s[3])

print, 'position of peak = ', r1d[where(var1d EQ max(var1d))]

samp_index = indgen(s[1]*s[2]*s[3])

r1d = temporary(r1d[samp_index])
var1d = temporary(var1d[samp_index])

; bin it
lrefine_max = max(tree[*].lrefine)
nXzones = params.ntopx*params.nxb*2^(lrefine_max - 1)
nYzones = params.ntopy*params.nyb*2^(lrefine_max - 1)

nbins = floor(sample*sqrt(float(nXzones)^2 + float(nYzones)^2)) + 1


rmin = 0.e0
rmax = max(r1d)

dr = (rmax - rmin)/float(nbins)
rBinEdge = (findgen(nbins+1))*dr + rmin

rBinCenters = 0.5*(rBinEdge - shift(rBinEdge,1)) + shift(rBinEdge,1)
rBinCenters = rBinCenters[1:*]

varBin = fltarr(nbins)
varBin[*] = 0.0

binCount = lonarr(nbins)
binCount[*] = 0


; do a spiffy binning.  Start by sorting the data in r order
indexSorted = sort(r1d)

r1d = temporary(r1d[indexSorted])
var1d = temporary(var1d[indexSorted])

oldIndex = (where(r1d GT rBinEdge[0]))[0]

print, 'using ', nbins, ' bins'

for i = 0l, nbins-1 do begin

    index = (where(r1d[oldIndex:*] GT rBinEdge[i+1]))[0]

    if (i mod 5000 EQ 0) then print, i

; if index = -1, then there are no radial values larger than the
; current bin limit.
    if (index EQ -1) then begin
        varBin[i] = total(var1d[oldIndex:(size(r1d))[1]-1])
        binCount[i] = (size(r1d))[1] - oldIndex
        if (i LT nbins-1) then varBin[i+1:*] = 0.
        goto, break_point
    endif

    varBin[i] = total(var1d[oldIndex:oldIndex+index])
    binCount[i] = index
    
    oldIndex = oldIndex + index
endfor

break_point: print, 'done binning'

for i = 0l, nbins-1 do begin
    
    if (binCount[i] NE 0) then varBin[i] = varBin[i]/binCount[i]

endfor

print, 'number of zones binned = ', total(binCount), (size(samp_index))[1]
return, varBin

end












