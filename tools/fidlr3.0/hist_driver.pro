pro hist_driver, FILE_INFO = fileInfo, $
                 VARIABLE_INFO = variable, $
                 OPTIONS = options, $
                 NBINS = nbins, $
                 HIST_SCALE = max_scale, $
                 OUTPUT = output, $
                 KNOWN_VARIABLES = varnames

; loop over all the files specified by fileInfo and produce the
; histograms

; FLASH files end with a 4-digit file number that is incremented as
; files are produced.  To get the base, strip the end off of the
; filename.  For now, ignore the possiblity of a .gz ending.
fileBase = strmid(fileInfo.prototype,0,strlen(fileInfo.prototype)-4)

pathEnd = strpos(fileBase, '/', /REVERSE_SEARCH)
fileBaseClean = strcompress(strmid(fileBase, pathEnd+1, $
                                   strlen(fileBase)-pathEnd))

;==============================================================================
; loop over the files
;==============================================================================
for ifile = fileInfo.startSuffix, fileInfo.endSuffix, fileInfo.step do begin

    filename = fileBase + string(ifile, format = '(i4.4)')

; ---- create the output filename ---------------------------------------------

    outfile = fileBaseClean

    outfile = outfile + strcompress(variable.name, /REMOVE_ALL) + $
      string(ifile, format = '(i4.4)')
   
    case output.type of
        1: outfile = outfile + '_hist.ps'
        2: begin
                
; newer version of IDL no longer support GIFs
            if (!VERSION.RELEASE LE 5.3) then begin
                outfile = outfile + '.gif'
            endif else begin
                outfile = outfile + '.png'
            endelse
        end
        else:
    endcase
 

    hist, filename, $
          variable.name, $
          BIN_MIN=float(variable.min), $
          BIN_MAX=float(variable.max), $
          NBINS=nbins, $
          HIST_SCALE=max_scale, $
          LOG=options.log, $
          POSTSCRIPT=output.type, $
          OUTFILE=outfile, $
          AUTO=variable.auto

endfor

end



