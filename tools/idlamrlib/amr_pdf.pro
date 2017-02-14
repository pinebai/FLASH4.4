function amr_pdf, amr, nbin, binctr=binctr, logctr=logctr, $
	max=max, min=min

; set keywords to defaults
if n_elements(max) eq 0 then max=(1.0+1.0d-10)*max_amr(amr)
if n_elements(min) eq 0 then min=min_amr(amr)

; set up bins
binlim=min*(max/min)^(findgen(nbin+1)/nbin)
pdf=dblarr(nbin)
binctr=min*(max/min)^((findgen(nbin)+0.5)/nbin)
meanval=amr_mean(amr)
logctr=alog(binctr/meanval)

; get the volume of the domain
vol=1.0d0
for i=0, amr.ndim-1 do vol=vol*(amr.boxmax[i]-amr.boxmin[i])

; go through the bins
for i=0,nbin-1 do begin
    masklo=amr_ge(amr, binlim[i])
    maskhi=amr_lt(amr, binlim[i+1])
    mask=amr_and(masklo, maskhi)
    pdf[i]=amr_sum(mask)/vol
    amr_free, masklo
    amr_free, maskhi
    amr_free, mask
endfor

; return
return, pdf

end