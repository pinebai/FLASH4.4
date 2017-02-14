;============================================================
; read a FLASH3 particles dump file into IDL
;  only supports HDF5, single precision
;============================================================
PRO read_particles, filename, $
                    PROPERTY_REQUESTED=propertyRequested, $
                    VERBOSE=verbose, $
                    NUM_PARTICLES=numParticles, $
                    PARTICLES=particles, $
                    STORED_PROPS = propertyFound

; check arguments
if n_elements(particles) gt 0 then undefine, particles
if (n_elements(verbose) EQ 0) then verbose = 0

; check to see which property to read; if propertyRequested is not defined, read them all
if  n_elements(propertyRequested) EQ 0 then begin
  propertyRequested = 'none'
endif else begin 
  ; check that variable name has four chars -- append with blanks
  ; if necessary.  Format string pads with blanks on the right if
  ; too short
  propertyRequested = string(propertyRequested, FORMAT='(A-4)')
endelse 

if (verbose) then print, ' Property requested = ', propertyRequested

;-----------------------------------------------------
;  begin reading
;---------------------------------------------------

file_identifier = H5F_OPEN(filename)
file_contents = h5_parse(filename)
tag_contents = tag_names(file_contents)
tag_where = where(tag_contents EQ "PARTICLE_NAMES")
; NOTE that tag_names removes white space and capitalizes everything!

if ( tag_where EQ -1) then begin
 if (verbose) then print, 'NO particles in the file'
 H5F_CLOSE, file_identifier  ; abort condition
 return
endif 

; ---------------------------------
;  at this point we are sure there are particles here, so find out
;  the particle property names
;------------------------------------

; get the names of the particle variables
dataset = H5D_OPEN(file_identifier, "particle names")
propertyFound =  H5D_READ(dataset)
H5D_CLOSE, dataset
propertyFound = reform(temporary(propertyFound))
; remove blank spaces
propertyFound = strcompress(propertyFound,/REMOVE_ALL)

; check for unique names
testUnique = propertyFound;
testUnique = testUnique[sort(testUnique)]
testUnique = testUnique[uniq(testUnique)]
sizeUnique = n_elements(testUnique)
sizeAll    = n_elements(propertyFound)
if (sizeUnique NE sizeAll) then begin
  print,' ERROR! some properties do not have unique names! Please check output of STORED_PROPS'
  print,'     Cannot continue....'
  numParticles = 0
  return
endif 

;-------------------------------------------------------------------
; we can also now set the number of particle properties
;------------------------------------------------------------------
numProperties = (size(propertyFound))[1]
if (propertyRequested NE 'none') then  begin 
  prop = (where(propertyFound EQ propertyRequested))[0]
  if (prop EQ -1) then begin
    print, 'ERROR: requested property ',propertyRequested,' not found in dataset'
    print, '       therefore reading in all variables'
  endif else begin 
    print, ' Reading in only the property ', propertyRequested
  endelse 
  nprop_readin = prop
endif else begin 
  prop = -1
  nprop_readin = -1
endelse 

;--------------------------------------------------------------------
; set up a structure -- this was done automatically in flash2 format
; thanks to dan for figuring out this execute string hack
;--------------------------------------------------------------------
; clever way to programatically define a structure whose attributes
; are named after the particle attributes found in the FLASH file
;--------------------------------------------------------------------
create_struct_string = 'particle = {'
for p = 0, numProperties-1 do begin
  comma_string = ', '
  if (p EQ 0) then comma_string = ' '
  create_struct_string = create_struct_string + comma_string + propertyFound[p] + ':double(0) '        
endfor 

create_struct_string = create_struct_string + '}'        
result = execute(create_struct_string)
if (result EQ  0) then print,' Error in creating particle structure!'

;-------------------------------------------------------------------
; read in the tracer particles
;-------------------------------------------------------------------
dataset = H5D_OPEN(file_identifier, "tracer particles")
particles_read = H5D_READ(dataset)
H5D_CLOSE, dataset

; 'particles_read' will be a 1d array if there is only one particle
; If this is the case, reform it to a 2D-array for our purposes below.
numDims = (size(particles_read))[0]
if (numDims EQ 1) then begin
  particles_read = reform(particles_read, (size(particles_read, /DIMENSIONS))[0], 1)
endif
print, "size(particles_read) is:", size(particles_read)

; 'particles_read' is a 2d array where the 1st dimension holds the
; different data "columns" for each particle and the 2nd dimension
; represents each particle's index. Therefore the length of the 2nd
; dimension is the number of particles in the simulation.
numParticles = (size(particles_read, /DIMENSIONS))[1]

; now assign the data to the structure.  gee this is tedious!
; make a big 1D-array out of the single "row" of the data structure
particles = replicate(particle,numParticles)
for i = 0L, numParticles-1L do begin 
  for j = 0, numProperties-1 do begin 
    assign_string = 'particles[i].' + propertyFound[j] + ' = particles_read[j,i]'
    result = execute(assign_string)
  endfor 
endfor 

H5F_CLOSE, file_identifier

end
