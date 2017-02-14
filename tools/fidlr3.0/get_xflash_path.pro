function get_xflash_path
; OK:  checked by LBR 2006/04/05
; get the xflash directory from the xflash_dir environmental variable
xflash_dir = getenv('XFLASH3_DIR')

if (strlen(xflash_dir) eq 0) then begin
  continue = dialog_message("You must set the environment variable 'XFLASH3_DIR' to the directory containing fidlr3.pro")
  widget_control, info.mainBase, /destroy
endif else begin
  ; make sure the path ends with a `/'
  if strmid(xflash_dir,strlen(xflash_dir)-1,1) NE '/' then $
    xflash_dir = xflash_dir + '/'
  return, xflash_dir
endelse

end
