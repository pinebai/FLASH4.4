pro psclose

; MAKE SURE WE HAVE POSTSCRIPT DEVICE OPEN...
if (!d.name ne 'PS') then begin
message, 'DEVICE is not set to PS!', /INFO
return
endif

; CLOSE THE POSTSCRIPT DEVICE...
device, /close
set_plot, 'X'

end; psclose 
