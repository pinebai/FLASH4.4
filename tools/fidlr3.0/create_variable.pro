; given the solution array (unk) and a variable name (atomic or
; derived), extract the desired variable (varname) via a simple copy
; or deriving it if necessary


function create_variable, unk, varname, debug

IF keyword_set(debug) THEN  print, 'looking to create variable', varname
case varname of

    'tot_vel': begin
        temp_arr = sqrt(reform(unk[var_index('velx'),*,*,*,*]^2 + $
                               unk[var_index('vely'),*,*,*,*]^2))
    end

    'snd_spd': begin
        temp_arr = sqrt(reform(unk[var_index('gamc'),*,*,*,*] * $
                               unk[var_index('pres'),*,*,*,*] / $
                               unk[var_index('dens'),*,*,*,*]))
    end
    
    'mach': begin
        temp_arr = sqrt(reform( (unk[var_index('velx'),*,*,*,*]^2 + $
                                 unk[var_index('vely'),*,*,*,*]^2) / $
                                (unk[var_index('gamc'),*,*,*,*] * $
                                 unk[var_index('pres'),*,*,*,*] / $
                                 unk[var_index('dens'),*,*,*,*])))
    end
    
    'int_ener': begin
        temp_arr = reform(unk[var_index('ener'),*,*,*,*] - .5* $
                          (unk[var_index('velx'),*,*,*,*]^2 + $
                           unk[var_index('vely'),*,*,*,*]^2))
    end
    
    'ekin/eint': begin
        temp_arr = reform(0.5*(unk[var_index('velx'),*,*,*,*]^2 + $
                               unk[var_index('vely'),*,*,*,*]^2) / $
                          (unk[var_index('ener'),*,*,*,*] - $
                           0.5*(unk[var_index('velx'),*,*,*,*]^2 + $
                                unk[var_index('vely'),*,*,*,*]^2)))
    end

    'egrav/eint': begin
        temp_arr = reform(unk[var_index('gpot'),*,*,*,*]/ $
                          (unk[var_index('ener'),*,*,*,*] - .5* $
                          (unk[var_index('velx'),*,*,*,*]^2 + $
                           unk[var_index('vely'),*,*,*,*]^2)))
    end

    else: begin
        IF keyword_set(debug) THEN print, 'all else failed, var_index = ', var_index(varname)
        temp_arr = reform(unk[var_index(varname),*,*,*,*])
    end
    
endcase
; hack: if there's only 1 block, e.g., if only 1 level of refinement,
; then make sure the temp_arr still has the block index...
; the rest of the code assumes the block index exists
unk_size = size(unk)
if (unk_size[2] EQ 1) then begin ; only one block
    temp_arr_size = size(temp_arr)
    dimensions = [1, temp_arr_size[1: temp_arr_size[0]]]
    temp_arr = reform(temp_arr, dimensions)
end


return, temp_arr
end
