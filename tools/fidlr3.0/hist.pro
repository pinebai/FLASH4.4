pro hist, filename, $
          var_name, $
          BIN_MIN=bin_min, $
          BIN_MAX=bin_max, $
          NBINS=nbins, $
          HIST_SCALE=max_scale, $
          LOG=log, $
          POSTSCRIPT=ipost, $
          OUTFILE=outfile, $
          AUTO=auto

if (n_elements(ipost) EQ 0) then ipost = 0
if (n_elements(outfile) EQ 0) then outfile = 'hist.ps'

if ipost EQ 1 then begin
    set_plot, 'PS'

    output_color = 'black'

    xsize = 7.5
    ysize = 9.5
    
    device, file = outfile, xsize=xsize, ysize=ysize, $
      xoff=0.5, yoff=.75, /inch, /color, bits_per_pixel = 8
    
endif else begin
    set_plot, 'X'

    output_color = 'white'
endelse


;------------------------------------------------------------------------------
; start off by determining the file type (HDF or HDF5) and then read
; in just the list of variables
;------------------------------------------------------------------------------
itype = determine_file_type(filename)

varnames = get_var_list(filename)

; store the number of native variable names -- some features will not
; use the derived variables
native_num = (size(varnames))[1]


;-----------------------------------------------------------------------------
; check if the variable name requested exists
;-----------------------------------------------------------------------------

if n_elements(var_name) EQ 0 then begin
    print, '******************************************************************'
    print, 'Error: please specify the variable name to bin when calling hist'
    print, 'valid variable names are: ', varnames
    print, '******************************************************************'
    return
endif


;-----------------------------------------------------------------------------
; read in the file -- just the variable desired
;-----------------------------------------------------------------------------

print, 'about to read file = ', filename, ' variable = ', var_name

read_amr, filename, VAR_NAME=var_name, $
  DATA=unk, TREE=tree, PARAMETERS=params

help, unk

;-----------------------------------------------------------------------------
; check to see whether minima and maxima were given, create the bins
;-----------------------------------------------------------------------------
if (n_elements(bin_min) EQ 0) then bin_min = min(unk[0,*,*,*,*])
if (n_elements(bin_max) EQ 0) then bin_max = max(unk[0,*,*,*,*])
if (n_elements(nbins) EQ 0) then nbins = 25
if (n_elements(log) EQ 0) then log = 0

if (n_elements(auto) EQ 0) then auto = 0

if (auto NE 0) then begin
    bin_min = min(unk[0,*,*,*,*])
    bin_max = max(unk[0,*,*,*,*])
endif


print, 'nbins = ', nbins, bin_max, bin_min

if (log EQ 0) then begin
    dbin = (bin_max - bin_min)/nbins
    bins = (findgen(nbins+1))*dbin + bin_min
endif else begin
    dbin = (alog10(bin_max) - alog10(bin_min))/nbins
    bins = 10.0^((findgen(nbins+1))*dbin + alog10(bin_min))
endelse

print, dbin, bins
count = lonarr(nbins)
count[*] = 0l


print, 'made it here'
print, count
print, bins

;------------------------------------------------------------------------------
; put all of the variable data into a vector and sort it
;------------------------------------------------------------------------------
nzones=0l
nzones = params.nxb*params.nyb*params.nzb*params.totBlocks
tvar = reform(unk[0,*,*,*,*],nzones)

indexSorted = sort(tvar)

tvar = temporary(tvar(indexSorted))

;------------------------------------------------------------------------------
; create the weights array -- store the refinement level of each zone
; in an array identical to tvar -- these will be used to ensure that
; a coarse block contributes more to the PDF than a fine block
;------------------------------------------------------------------------------
weights = replicate(1, params.totBlocks, params.nxb, params.nyb, params.nzb)

if (params.nyb EQ 1 AND params.nzb EQ 1) then begin
 ndim = 1
endif else if (params.nzb EQ 1) then begin
 ndim = 2
endif else begin
 ndim = 3
endelse

lrefine_max = max(tree[*].lrefine)

for blk = 0l, params.totBlocks-1 do begin
    weights[blk,*,*,*] = 2^((lrefine_max-tree[blk].lrefine)*ndim)
endfor

; sort the weights too so they match their counterpart in tvar
weights = reform(temporary(weights), nzones)
weights = temporary(weights(indexSorted))

;------------------------------------------------------------------------------
; do the binning
;------------------------------------------------------------------------------

; set the initial index into tvar at the location were bin_min falls
oldIndex = (where(tvar GE bins[0]))[0]

print, 'oldIndex = ', oldIndex

for i = 0, nbins-1 do begin

; make where's job easier as we go along by reducing the number of 
; elements it needs to look at, since we are sorted
    index = (where(tvar[oldIndex:*] GT bins[i+1]))[0]

; if index = -1, then there are no values in tvar that are greater
; than the current bin limit.  All the remaining elements are in the 
; current bin, and the counts for all other bins are 0.
    if (index EQ -1) then begin
        count[i] = total(weights[oldIndex:(size(tvar))[1]-1])
        if (i LT nbins-1) then count[i+1:*] = 0
        goto, break_point
    endif

    print, 'i = ', i, ' oldIndex = ', oldIndex, ' index = ', index

    count[i] = total(weights[oldIndex:index+oldIndex])
    
    oldIndex = oldIndex + index
endfor

; for some reason, IDL 5.3 does not have the break statement, so we
; need to resort to a goto
break_point: print, 'done binning'

; now normalize
normal_factor = total(count)
count = temporary(count)/normal_factor

if (max_scale EQ -1) then begin
    max_count = 1.1*max(count)
endif else begin
    max_count = max_scale
endelse

;------------------------------------------------------------------------------
; dump out the results
;------------------------------------------------------------------------------
if (log EQ 0) then begin
    center = (bins - shift(bins,1))/2. + shift(bins,1)
    center = center[1:*]
endif else begin
    center = 10.0^((alog10(bins) - alog10(shift(bins,1)))/2. + alog10(shift(bins,1)))
    center = center[1:*]
endelse

help, bins
help, center

for i = 0, nbins-1 do begin
    print, 'range: ', bins[i], ' to ', bins[i+1], ':', count[i], center[i]
endfor


; get the colortable
iclrmap = color_index('Grayscale', MIN_VALUE=colorMin, MAX_VALUE=colorMax)
loadct, iclrmap, FILE = xflash_dir + 'flash_colors.tbl', /SILENT

colors = replicate(color(output_color), nbins)
print, colors
;bar_plot, count, COLORS=colors, BARNAMES=string(center)

;plot, center, count

if (log EQ 0) then begin

; create the axes -- the data is normalized to fall between 0 and 1
    plot, [bins[0],bins[nbins-1]], [0.0,max_count], $
      PSYM = 1, XSTYLE=1, YSTYLE=1, $
      XRANGE=[bins[0],bins[nbins]], YRANGE=[0,max_count]

; loop over all the bins and plot polygons
    dsmall = dbin*0.05

    for i = 0, nbins-1 do begin
        polyfill, [bins[i]+dsmall, bins[i]+dsmall, bins[i+1]-dsmall, $
                   bins[i+1]-dsmall, bins[i]+dsmall], $
              [0, count[i], count[i], 0, 0]
    endfor

endif else begin

; create the axes -- the data is normalized to fall between 0 and 1
    plot, [alog10(bins[0]),alog10(bins[nbins-1])], [0.0,max_count], $
      PSYM = 1, XSTYLE=1, YSTYLE=1, $
      XRANGE=[alog10(bins[0]),alog10(bins[nbins])], YRANGE=[0,max_count]

; loop over all the bins and plot polygons
    for i = 0, nbins-1 do begin
        polyfill, [alog10(1.005*bins[i]), alog10(1.005*bins[i]), alog10(0.995*bins[i+1]), alog10(0.995*bins[i+1]), alog10(1.005*bins[i])], $
              [0, count[i], count[i], 0, 0]
    endfor

endelse 
  
if ipost EQ 1 then device, /close

end





