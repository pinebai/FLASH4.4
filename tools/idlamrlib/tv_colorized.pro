pro tv_colorized,image,xin,yin,cin,overridect=overridect,_extra=extra
;oldcolor = intarr(256,3)
;tvlct,oldcolor,/get
;loadct,colormap
x=0
if keyword_set(xin) then x = xin
y=0
if keyword_set(yin) then y = yin
channel=0
if keyword_set(cin) then channel = cin

c = intarr(256,3)
backup=c
tvlct,backup,/get
ct=13
if keyword_set(overridect) then begin
    ct = overridect
endif
loadct,ct,/silent
tvlct,c,/get
tvlct,backup
imc = [[[image]],[[image]],[[image]]]
for i=0,2 do imc[*,*,i] = c[image[*,*],i]
tv,imc,x,y,channel,true=3,_extra=extra
;tvlct,oldcolor
end
