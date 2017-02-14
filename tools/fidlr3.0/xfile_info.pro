; contour widget event handler
pro xfile_info_event, ev

widget_control, ev.top, get_uvalue=info
widget_control, ev.id, get_uvalue=uval

; take action pased on which widget was changed
case uval of

    'done': begin
        widget_control, ev.top, /destroy
    end

    else:
endcase

end




pro xfile_info, filename

file_information, filename, /SILENT, $
  HDF_TYPE=hdfType, $
  CORNERS=corners, $
  DATE=date, $
  FILE_TYPE=fileType, $
  RUN_COMMENT=runComment, $
  FLASH_VERSION=flashVersion, $
  DIMENSIONALITY=ndim, $
  NXB=nxb, NYB=nyb, NZB=nzb, $
  PRECISION=precision, $
  NUMBER_OF_VARIABLES=nvar, $
  NUMBER_OF_PARTICLES=numParticles, $
  BUILD_DATE=buildDate, $
  BUILD_DIR=buildDir, $
  BUILD_MACHINE=buildMachine, $
  SETUP_CALL=setupCall




; create the main widget base
main_base = widget_base(/column, title = 'file information')

; divide this base into vertical sections
filename_base = widget_base(main_base, /column, /frame)
info_base = widget_base(main_base, /column, /frame)
build_base = widget_base(main_base, /column, /frame)
format_base = widget_base(main_base, /column, /frame)
var_base = widget_base(main_base, /column, /frame)
dimens_base = widget_base(main_base, /column, /frame)
end_base = widget_base(main_base, /column)

ifilename = widget_label(filename_base, /align_left, $
                         value='filename:        ' + filename)

iversion = widget_label(info_base, /align_left, $
                        value='FLASH version:   ' + flashVersion[0])

idate = widget_label(info_base, /align_left, $
                     value='execution date:  ' + strcompress(date[0]))

iblank4 = widget_label(info_base, /align_left, value=' ')

icomment = widget_label(info_base, /align_left, $
                        value='run comment:     ' + $
                        strcompress(runComment[0]))



ibdate = widget_label(build_base, /align_left, $
                      value='build date:      ' + strcompress(buildDate[0]))

ibdir = widget_label(build_base, /align_left, $
                     value='build directory: ' + strcompress(buildDir[0]))

ibmachine = widget_label(build_base, /align_left, $
                     value='build machine:   ' + strcompress(buildMachine[0]))

isetup = widget_label(build_base, /align_left, $
                     value='setup call:      ' + strcompress(setupCall[0]))



ihdf = widget_label(format_base, /align_left, $
                    value='format:          ' + hdfType[0])

invar = widget_label(var_base, /align_left, $
                     value='# of variables: ' + strcompress(string(nvar)))

inpart = widget_label(var_base, /align_left, $
                      value='# of particles: ' + $
                      strcompress(string(numParticles)))

ipres = widget_label(var_base, /align_left, $
                     value='precision:       ' + precision)

icorners = widget_label(var_base, /align_left, $
                        value='corner data:     ' + corners)


idim = widget_label(dimens_base, /align_left, $
                    value='# of dimens:    ' + strcompress(string(ndim)))

inxb = widget_label(dimens_base, /align_left, $
                    value='nxb:            ' + strcompress(string(nxb)))

inyb = widget_label(dimens_base, /align_left, $
                    value='nyb:            ' + strcompress(string(nyb)))

inzb = widget_label(dimens_base, /align_left, $
                    value='nzb:            ' + strcompress(string(nzb)))



done = widget_button(end_base, value = 'Done', uvalue = 'done')


;------------------------------------------------------------------------
; draw the widget
;------------------------------------------------------------------------
widget_control, main_base, /realize


; setup a structure to hold the widget information
info = {main_base:main_base, $  ; main base
        info_base:info_base, $  ; information base
        end_base:end_base, $    ; end base widget
        done:done}              ; button for done


; register the info structure with the widget base
widget_control, main_base, set_uvalue=info, /no_copy

xmanager, 'xfile_info', main_base

end
