;--------------------------------------------------------------------------
; function sci_notat:
;
;    return a string with the argument in scientifix notation ready
;    for displaying on the graphics screen with xyouts
;--------------------------------------------------------------------------
function sci_notat, number

number = double(number)

; get the exponent for the scientific notation representation
; use floor to convert to an integer because of round-off errors 
;    ex, print, alog10(1.d-4)  + 4.d0 is NE 0, whereas 
;        print, alog10(1.d-12) + 12.d0 = 0, so the result is
;        inconsistent -- This only affects sphere, not orion.
;
exponent = floor(alog10(abs(number)))

if abs(exponent) GT 2 then begin

; compute the prefix
    multiplier = number / 10.d0^exponent

    number_str = string(multiplier, format = '(f4.1)') + 'x' + '10!U' + $
      strcompress(string(exponent, format = '(i8)'), /remove_all) + '!N'

endif else begin
    
    number_str = strcompress(string(number, format = '(f8.2)'), /remove_all)

endelse

return, number_str
end
