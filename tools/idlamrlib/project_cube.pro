
; project a cube containing 'value' centered on 'center' with size 'size' onto a regular
; grid 'img' framed by imglim under transformation QSC, subsampling
; 'subsample'^3 times
pro project_cube,img,imglim,QSC, density, center, size, ssinput
; make sure cube can project at all
;poscam = (QSC ## [center,double(1.0)])[0:2]
;if (total(poscam lt imglim[0,*]) or total(poscam  gt imglim[1,*])) then begin
;    print,'In PROJECT_CUBE cube ',center,size,' falls outsize camera',imglim
;    return
;endif

ximgsize = (size(img))[1]
yimgsize = (size(img))[2]
xcellsize = double((imglim[1,0]-imglim[0,0]))/ double(ximgsize)
ycellsize = double((imglim[1,1]-imglim[0,1]))/ double(yimgsize)

; increase subsampling if amr cell is larger than projection cell
sscorr = 1
if size[0] gt xcellsize then begin
    sscorr = max([sscorr,fix(size[0] / xcellsize + 0.5)])
endif
if size[0] gt ycellsize then begin
    sscorr = max([sscorr,fix(size[0] / ycellsize + 0.5)])
endif
if size[1] gt ycellsize then begin
    sscorr = max([sscorr,fix(size[1] / ycellsize + 0.5)])
endif
if size[1] gt xcellsize then begin
    sscorr = max([sscorr,fix(size[1] / xcellsize + 0.5)])
endif
subsample = ssinput
if sscorr gt 1 then subsample = ssinput * sscorr
subsample = 1.0
if subsample gt 100 then begin
    print,'Subsample wanted to be ',subsample,' but was clipped to 100'
    subsample = 100
endif

ss3 = double(subsample*subsample*subsample)
fractvolume = double(size[0])*double(size[1])*double(size[2])/ss3
totalvolume = double(size[0])*double(size[1])*double(size[2])
if (subsample / 2) eq 1 then begin
    ; odd subsampling
    ssindlo = -(subsample -1)/2.0
    sindhi = (subsample-1)/2.0
endif else begin
    ; even subsampling
    ssindlo = -subsample/2.0 + 0.5
    ssindhi = subsample/2.0 - 0.5
endelse
ssinc = double(1.0)
nloops = ss3
if nloops gt double(1000000.0) then begin
    ssinc = double(nloops / 1000000.0)
endif

;poscvect = [size * (double(ssinc) / double(subsample)), 1.0]
;poscvect = (QSC ## poscvect)[0:2]

poscvectx = size*double(ssinc)/double(subsample)*(QSC##[1,0,0,1])[0:2]
poscvecty = size*double(ssinc)/double(subsample)*(QSC##[0,1,0,1])[0:2]
poscvectz = size*double(ssinc)/double(subsample)*(QSC##[0,0,1,1])[0:2]
poscstart = (QSC##([center,1.0]))[0:2]

;poscam = poscstart + ssi*poscvectx + ssj*poscvecty + ssk*poscvecty
;BAD  poscam = poscstart + [ssi,ssj,ssk]*(QSC##(size*ssinc/subsample))

imcrementmass = double(density)*double(fractvolume)/double(xcellsize*ycellsize)
for ssi = double(ssindlo), double(ssindhi),double(ssinc) do begin
    dx = ssi*poscvectx
    for ssj = ssindlo, ssindhi,ssinc do begin
        dy = ssj*poscvecty
        for ssk = ssindlo, ssindhi,ssinc do begin
            dz = ssk*poscvectz
            ;npi = float(ssi) / float(subsample)
            ;npj = float(ssj) / float(subsample)
            ;npk = float(ssk) / float(subsample)
            ;possim = [npi*size[0], npj*size[1], npk*size[2]]+center
            ;poscam = (QSC ## [possim,double(1.0)])[0:2]
            poscam = poscstart + dx+dy+dz
            ;if total(poscam lt imglim[0,*]) then begin
                ; falls outside low end
            ;endif else if total(poscam gt imglim[1,*]) then begin
                ; falls outside high end
            ;endif else begin
                ; falls inside camera box
            xind = (poscam[0]-imglim[0,0]) / xcellsize
            yind = (poscam[1]-imglim[0,1]) / ycellsize
            if ( (xind ge 0 and yind ge 0) and $
                 (xind lt ximgsize and yind lt yimgsize) ) then begin
                ;img[xind,yind] = img[xind,yind] + value/ss3/double(xcellsize*ycellsize)
                img[xind,yind] = img[xind,yind] + $
                  imcrementmass
            endif
                ; cell# = position-boxstart / cellsize
            ;endelse
        endfor
    endfor
endfor
      
end ; PROJECT_CUBE
