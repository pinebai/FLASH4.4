function get_dat_vars, filename

openr, lun, filename, /get_lun

; read the header line
line = ' '
readf, lun, line

; the first character is a '#', eliminate it.  Also eliminate all the
; spaces at the end

line = strtrim(strmid(line, 1))

varname = strarr(100)

name = 0

; some variable names have a space in them (how unfortunate), so we
; search for '  ' as the delimiter
i = strpos(line, '  ')

; extract the variable names
while (i GE 0) do begin

    print, i
    varname[name] = strmid(line, 0, i)
    print, varname[name]

    line = strtrim(strmid(line, i+1), 1)

    name = name + 1

    i = strpos(line, '  ')
endwhile

; there is one more variable left
varname[name] = strtrim(line, 1)

; the first variable is time -- skip this, as it is always the
; dependent variable
varname = varname[1:name]

return, varname

end
