pro rho1d
;;;;;;;;use float ncreat!!!!!!!!! otherwise short upto only 2^15
   nfile = 100
   nstart = 0

  for k1=0,nfile do begin
     if (k1 LT 10) then begin
         infile  = strcompress(string(format='("dumps/sb_hdf5_chk_000"+(I))',k1),/remove_all)
     endif else begin
     if (k1 LT 100)  then begin
         infile  = strcompress(string(format='("dumps/sb_hdf5_chk_00"+(I))',k1),/remove_all)
     endif else begin
         infile  = strcompress(string(format='("dumps/sb_hdf5_chk_0"+(I))',k1),/remove_all)
     endelse
     endelse

    datfile = strcompress(string(format='("dumps/rho_"+(I))',k1),/remove_all)

    openw, lun, datfile,/get_lun

    a=loaddata(infile,'dens',XCOORDS=x, TIME=t)

    tmp = size(a)
    nx = tmp[1]
    dx = x[1]-x[0]

    print, 'working on ', infile
    print,datfile
    print,t

    printf,lun, FORMAT='("#"(1I,1E20.10))',nx,t
   
    for i=0,nx-1 do begin 
       printf,lun,FORMAT='(2E20.10)',x[i],a[i]
    endfor

    close, lun
    free_lun, lun
endfor 

return
end
