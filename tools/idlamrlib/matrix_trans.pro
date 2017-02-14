function PAAA_M, $
                 convert, $
                 point, angles, $
                 matrix
if convert eq 'toM' then begin
    azimuth = angles[0]
    altitude = angles[1]
    pitch = angles[2]
    QT = [[1.0,0.0,0.0,-point[0]], $
          [0.0,1.0,0.0,-point[1]], $
          [0.0,0.0,1.0,-point[2]], $
          [0.0,0.0,0.0,1.0]]
    QA = [[ cos(azimuth),sin(azimuth),0,0],$
          [-sin(azimuth),cos(azimuth),0,0],$
          [            0,           0,1,0],$
          [            0,           0,0,1 ]
    QP = [[1,            0,             0,0],$
          [0, cos(altitude),sin(altitude),0],$
          [0,-sin(altitude),cos(altitude),0],$
          [0,            0,             0,1]]
    QR = [[ cos(pitch),sin(pitch),0,0],$
          [-sin(pitch),cos(pitch),0,0],$
          [          0,         0,1,0],$
          [          0,         0,0,1]]
    matrix = QR##(QP##(QA##QT))
    
endif
if convert eq 'toPAAA' then begin
    ctheta = acos(matrix[2,2])
    if ctheta gt 1.0 then ctheta = 1.0
    if ctheta lt -1.0 then ceheta = -1.0
    theta = acos(ctheta)
    stheta = sqrt(1.0 - ctheta^2)
    if stheta ne 0.0 then begin
        sphi = matrix[2,1] / stheta
        cphi = matrix[0,2] / stheta
        if abs(sphi) gt abs(cphi) then begin
            phi = acos(cphi)
        endif else begin
            phi = asin(sphi)
        endelse
    endif else begin
        ; theta ~= 0.0
    endelse
endif
end

function AA_M, $
               convert, $
               axis, angle, $
               matrix

if convert eq 'toM' then begin
    print,'Converting to matrix'
    matrix = dblarr(3,3)
    x = axis[0]
    y = axis[1]
    z = axis[2]
    corr = sqrt(x^2+y^2+z^2)
    if corr ne 1.0 then begin
        x = x/corr
        y = y/corr
        z = z/corr
    endif
    matrix[0,0] = 1+(1-cos(angle))*(x*x-1)
    matrix[1,0] = -z*sin(angle)+(1-cos(angle))*x*y
    matrix[2,0] = y*sin(angle)+(1-cos(angle))*x*z
    matrix[0,1] = z*sin(angle)+(1-cos(angle))*x*y
    matrix[1,1] = 1 + (1-cos(angle))*(y*y-1)
    matrix[2,1] = -x*sin(angle)+(1-cos(angle))*y*z
    matrix[0,2] = -y*sin(angle)+(1-cos(angle))*x*z
    matrix[1,2] = x*sin(angle)+(1-cos(angle))*y*z
    matrix[2,2] = 1 + (1-cos(angle))*(z*z-1)
endif else if convert eq 'toAA' then begin
    print,'Converting to axis and angle'
    ; (matrix[2,0]+matrix[0,2])/2.0 = (1-cos(angle))*x*z
    ; (matrix[2,0]-matrix[0,2])/2.0 = y*sin(angle)
    ; (matrix[1,0]+matrix[0,1])/2.0 = (1-cos(angle))*x*y
    ; (matrix[1,0]-matrix[0,1])/2.0 = -z*sin(angle)
    ; (matrix[2,1]+matrix[1,2])/2.0 = (1-cos(angle))*y*z
    ; (matrix[2,1]-matrix[1,2])/2.0 = -x*sin(angle)
    if total(matrix eq [[1,0,0],[0,1,0],[0,0,1]]) eq 9 then begin
        x=1
        y=0
        z=0
        angle=0
        axis = [x,y,z]
        return,1
    endif
    if total(matrix eq 0.0) eq 9 then begin
        x=0
        y=0
        z=0
        angle=0
        axis = [x,y,z]
        return,1
    endif
    xovery = (matrix[1,2]-matrix[2,1]) / (matrix[2,0]-matrix[0,2])
    yoverz = (matrix[0,2]-matrix[2,0]) / (matrix[1,0]-matrix[0,1])
    zoverx = (matrix[1,0]-matrix[0,1]) / (matrix[2,1]-matrix[1,2])
    a = xovery
    b = yoverz
    c = zoverx
    erra = abs(alog(a))
    errb = abs(alog(b))
    errc = abs(alog(c))
    if (erra le errb) and (erra le errc) then begin
        ; x-y relationship if best-conditioned
        x = 1.0
        y = 1.0 / a
        z = c
        corr = sqrt(x^2 + y^2 + z^2)
        x = x/corr
        y = y/corr
        z = z/corr
        sa = (matrix[1,0]-matrix[0,1])/2.0/(-z)
        sa = min([sa,1.0])
        sa = max([sa,-1.0])
        angle = asin( sa )
    endif else if (errb le erra) and (errb le errc) then begin
        ; y-z relationship is best conditioned
        y = 1.0
        z = 1.0 / b
        x = a
        corr = sqrt(x^2 + y^2 + z^2)
        x = x/corr
        y = y/corr
        z = z/corr
        sa = (matrix[2,1]-matrix[1,2])/2.0/(y)
        sa = min([sa,1.0])
        sa = max([sa,-1.0])
        angle = asin( sa )
    endif else if (errc le erra) and (errc le errb) then begin
        ; z-x relationship is best conditioned
        z = 1.0
        x = 1.0/c
        y = b
        corr = sqrt(x^2 + y^2 + z^2)
        x = x/corr
        y = y/corr
        z = z/corr
        sa = (matrix[2,0]-matrix[0,2])/2.0/(-x)
        sa = min([sa,1.0])
        sa = max([sa,-1.0])
        angle = asin( sa )
    endif
    
    if (not a) and (not b) then begin
        ; x=y=0
        x=0
        y=0
        z=1
        sa = (matrix[1,0]-matrix[0,1])/2.0/(-z)
        sa = min([sa,1.0])
        sa = max([sa,-1.0])
        angle = asin( sa )
    endif
    if (not b) and (not c) then begin
        x=1
        y=0
        z=0
        sa = (matrix[2,1]-matrix[1,2])/2.0/(y)
        sa = min([sa,1.0])
        sa = max([sa,-1.0])
        angle = asin( sa )
    endif
    if (not c) and (not a) then begin
        x=0
        y=1
        z=0
        sa = (matrix[2,0]-matrix[0,2])/2.0/(-x)
        sa = min([sa,1.0])
        sa = max([sa,-1.0])
        angle = asin( sa )
    endif       

    axis=[x,y,z]
    ;print,bork
endif

end

function PAA_M,convert,point,axis,angle,matrix
; converts either a rotation point,axis and angle into a matrix, or a
; matrix into a point,axis and angle
if convert eq 'toPAA' then begin
    rm = matrix[0:2,0:2]
    rm4 = dblarr(4,4)
    rm4[0:2,0:2] = rm
    rm4[3,3] = 1.0
    irm4 = invert(rm4)
    tm = irm4##matrix
    point=-1.0*transpose(tm[3,0:2])
    axis = dblarr(3)
    angle = 0.0
    junk = aa_m('toAA',axis,angle,rm)
    print,bork
endif

if convert eq 'toM' then begin
    rm = dblarr(3,3)
    junk = aa_m('toM',axis,angle,rm)
    tm = [ [1.0,0.0,0.0,-point[0]], $
           [0.0,1.0,0.0,-point[1]], $
           [0.0,0.0,1.0,-point[2]], $
           [0.0,0.0,0.0,      1.0]]
    rm4 = dblarr(4,4)
    rm4[0:2,0:2] = rm
    rm4[3,3] = 1.0
    matrix = rm4##tm
    print,bork
endif

end

function PPA_M,convert,point1,point2,angle,matrix
; converts a rotation about the line between 2 points into a matrix,
; or back again
if convert eq 'toM' then begin
    point = (point1+point2)/2.0
    centerline = (point1-point2)/2.0
    ncenter = centerline / sqrt(centerline[0]^2 + centerline[1]^2 + centerline[2]^2)
    dest = [cos(angle),sin(angle),0.0]
    rotangle = acos( dest[0]*ncenter[0] + dest[1]*ncenter[1] + dest[2]*ncenter[2])
    rotaxis = crossp(ncenter,dest)
    ; need to rotate about point so that centerline paralell to dest
    junk = paa_m('toM',point,rotaxis,rotangle,matrix)
endif

if convert eq 'toPPA' then begin
    rotaxis=dblarr(3)
    rotpoint=dblarr(3)
    rotangle = 0.0
    junk = paa_m('toPAA',rotpoint,rotaxis,rotangle,matrix)
    im = invert(matrix)
    ncenter = im##([rotpoint[0]+0.5,rotpoint[1],rotpoint[2],1.0]) - im##([rotpoint[0]-0.5,rotpoint[1],rotpoint[2],1.0]) 
    lp = sqrt(point[0]^2 + point[1]^2 + point[2]^2)/2.0
    ; so that you don't get something like X+1=X-1 for X large
    point1 = point-lp*axis
    point2 = point+lp*axis
endif
end

function PPP_M,convert,point1,point2,point3,matrix
if convert eq 'toM' then begin
    point = (point1+point2+point3)/3.0
    v1 = point-point1
    v2 = point-point2
    axis = crossp(v1,v2)
    if total(axis eq 0.0) eq 3 then begin
        v3 = point-point3
        axis = crossp(v1,v3)
    endif
endif

if convert eq 'toPPP' then begin
endif
end
