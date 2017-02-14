pro psopen, filename, $
LANDSCAPE=landscape, $
PORTRAIT=portrait, $
XSIZE=xsize, $
YSIZE=ysize, $
FONTSIZE=fontsize, $
INCHES=inches, $
COLOR=color, $
ENCAPSULATED=encapsulated, $
_EXTRA=_extra

; LOOK FOR A POSTSCRIPT EXTENSION...
filename = strtrim(filename,2)
dotpos = strpos(filename,'.',/reverse_search)
if (dotpos ge 0) then begin
filelen = strlen(filename)
extension = strupcase(strmid(filename,dotpos,filelen-dotpos))
if (extension eq '.PS') OR (extension eq '.EPS') $
then filename = strmid(filename,0,dotpos)
endif

; APPEND THE PROPER EXTENSION...
if keyword_set(ENCAPSULATED) $
then filename = filename+'.eps' $
else filename = filename+'.ps'

; SET THE OUTPUT DEVICE FOR IDL GRAPHICS TO POSTSCRIPT...
set_plot,'PS', COPY=keyword_set(COLOR), INTERPOLATE=keyword_set(COLOR)

; LANDSCAPE, PORTRAIT, BITS_PER_PIXEL, ENCAPSULATED , INCHES KEYWORDS
; STICK UNTIL YOU CHANGE THEM... MUST SET EXPLICITLY EVERY TIME...

; IF NEITHER LANDSCAPE NOR PORTRAIT IS SET, ASSUME PORTRAIT...
if not(keyword_set(LANDSCAPE) OR keyword_set(PORTRAIT)) $
then portrait=1
if not keyword_set(BITS_PER_PIXEL) $
then bits_per_pixel=8
if not keyword_set(FONT) $ ; DEFAULT IS COURIER
then font='Courier'
if not keyword_set(FONTSIZE) $ ; DEFAULT IS 12 POINT...
then fontsize=12
if not keyword_set(CHARSIZE) $ ; DEFAULT IS 140 PIXELS WIDE, 180 TALL...
then charsize=[140,180]
if not keyword_set(INCHES) AND (N_elements(INCHES) eq 0) $
then inches=1



; A TWO-ELEMENT VECTOR TO SPECIFY THE FONT SIZE AND LINE SPACING
; (LEADING) OF VECTOR AND TRUETYPE FONTS, AND THE LINE SPACING OF
; DEVICE FONTS. THE WAY THAT THE VALUE OF THIS VECTOR DETERMINES
; CHARACTER SIZE IS NOT COMPLETELY INTUITIVE.

;!!!!!!!!!!!!!
; CAN WE USE BOTH SET_CHAR_SIZE *AND* FONT_SIZE ????
;!!!!!!!!!!!!!

; IF INCHES EXPLICITLY SET TO ZERO, SIZES ARE IN CENTIMETERS...
inch2cm = keyword_set(INCHES) + 2.54*(inches eq 0)

; SET UP LANDSCAPE ORIENTATION...
if keyword_set(LANDSCAPE) then begin

if not keyword_set(XSIZE) OR not keyword_set(YSIZE) then begin
xsize = 10.5*inch2cm
ysize = 8.0*inch2cm
endif

; DEFAULT IS TO KEEP THE IMAGE CENTERED ON PAGE...
; IF XOFFSET AND YOFFSET ARE SENT IN VIA _EXTRA, THOSE
; VALUES WILL OVERRIDE THE ONES BELOW... THAT'S PREFERABLE..
xoffset = (8.5*inch2cm-ysize)/2.
yoffset = (11.*inch2cm-xsize)/2. + xsize

endif

; SET UP PORTRAIT ORIENTATION...
if keyword_set(PORTRAIT) then begin

if not keyword_set(XSIZE) OR not keyword_set(YSIZE) then begin
xsize = 8.0*inch2cm
ysize = 10.5*inch2cm
endif
xoffset = (8.5*inch2cm-xsize)/2.
yoffset = (11.*inch2cm-ysize)/2.

endif

; SET THE DEVICE PARAMETERS...
device, FILENAME=filename, $
LANDSCAPE=keyword_set(LANDSCAPE), $
PORTRAIT=keyword_set(PORTRAIT), $
XSIZE=xsize, $
YSIZE=ysize, $
XOFFSET=xoffset, $
YOFFSET=yoffset, $
COLOR=keyword_set(COLOR), $
INCHES=keyword_set(INCHES), $
BITS_PER_PIXEL=bits_per_pixel, $
ENCAPSULATED=keyword_set(ENCAPSULATED), $
SET_FONT=font, /TT_FONT, $
FONT_SIZE=fontsize, SET_CHARACTER_SIZE=charsize, $
_EXTRA=_extra

end; psopen 
