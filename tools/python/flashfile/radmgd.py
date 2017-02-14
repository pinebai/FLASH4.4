import numpy as np

def numGroups(ff):
    count = 0
    while True:
        count += 1
        if not ff.exists("r%03d" % count):
            break
        
    return count - 1

def computeFlux(ff, opac, alpha=1):
    ng = numGroups(ff)
    ni = len(ff.leaves)*ff.nb[0]
    flux = np.zeros(ni)

    for n in xrange(ng): 
        erad = ff.getVar("r%03d" % (n+1))
        dens = ff.getVar("dens")
        urad = erad*dens
        coord = ff.coord()

        for i in xrange(ni):
            # Compute the gradient of u for group g, cell i:
            grad = 0.0
            if i == 0:
                dxr = coord[i+1]-coord[i]
                grad = (urad[i+1]-urad[i])/dxr

            elif i == ni-1:
                dxl = coord[i]-coord[i-1]
                grad = (urad[i]-urad[i-1])/dxl

            else:            
                dxr = coord[i+1]-coord[i]
                dxl = coord[i]-coord[i-1]

                gradr = (urad[i+1]-urad[i])/dxr
                gradl = (urad[i]-urad[i-1])/dxl
                grad = 0.5*(gradr + gradl)

            # Compute the flux limited diffusion coefficient:
            dcoef = 1.0/(3.0*opac(n,i) + abs(grad)/(3.0e+10*alpha*urad[i]+1.0e-100))

            # Compute the flux and add it to the total for this cell:
            flux[i] += -dcoef*grad

    return flux
