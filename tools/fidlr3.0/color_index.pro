; store the names of the colortables defined for xflash in
; flash_colors.tbl.  Given the string name of the color, return the
; index into the table

function color_index, name, GET_NAMES=colors, $
                      MIN_VALUE=colorMin, MAX_VALUE=colorMax

colors = ['Red Temp', $
          'Blue -> Red', $
          'Rainbow', $
          'Fire', $
          'Menacing', $
          'RT', $
          'Grayscale']

index = (where(name EQ colors))[0]

; the first few colors in the table are reserved for xflash standard
; colors.   Set the range of usable colors in the table.

case index of 
    0: begin
        colorMin = 255
        colorMax = 120
    end
    1: begin
        colorMin = 255
        colorMax = 20
    end
    else: begin
        colorMin = 255
        colorMax = 12
    end
endcase

return, index

end

