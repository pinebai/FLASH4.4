from FlashFile import FlashFile
import numpy as np

class FlashFile2d(FlashFile): 
    
    def __init__(self, fn):
        FlashFile.__init__(self,fn)

    def uniformGrid(self, var):
        
        xmin = self.realScalar("xmin", rtp=True)
        xmax = self.realScalar("xmax", rtp=True)
        ymin = self.realScalar("ymin", rtp=True)
        ymax = self.realScalar("ymax", rtp=True)

        # Get lrefine_max/min:
        # TODO: this will break for NOFBS
        lrmin = 1
        lrmax = 1
        try:
            lrmin   = self.integerScalar("lrefine_min", rtp=True)
            lrmax   = self.integerScalar("lrefine_max", rtp=True)
        except:
            pass

        nblockx = self.integerScalar("nblockx", rtp=True)
        nblocky = self.integerScalar("nblocky", rtp=True)

        ncx = (nblockx * self.nb[0] * 2**(lrmax-1))
        ncy = (nblocky * self.nb[1] * 2**(lrmax-1))

        # Create array which will store:
        arr = np.zeros((ncx, ncy))

        # Load all of the data:
        data = self.h5file.getNode("/" + var)[:,0,:,:]
        
        # Load an array containing the refinement level of each block:
        ref = self.h5file.getNode("/refine level")[:]

        # Loop over leaf blocks:
        for i,l in enumerate(self.leaves):
            xs = self.bboxes[l,0,0]
            ys = self.bboxes[l,1,0]
            xe = self.bboxes[l,0,1]
            ye = self.bboxes[l,1,1]

            istart = int(round(xs/(xmax-xmin)*ncx))
            jstart = int(round(ys/(ymax-ymin)*ncy))
            iend   = int(round(xe/(xmax-xmin)*ncx))
            jend   = int(round(ye/(ymax-ymin)*ncy))

            # print i,l,istart, jstart, iend, jend, ref[l]
            
            for m in xrange(istart,iend):
                for n in xrange(jstart,jend):
                    fac = 2**(lrmax-ref[l])
                    
                    lm = (m-istart)/fac
                    ln = (n-jstart)/fac

                    # print m,n,lm,ln, lrmax, ref[l], data.shape

                    arr[m,n] = data[l, ln, lm]

        return arr


    def getVar(self, name):
        var = self.h5file.getNode("/" + name)[:,0,0,:]
        ret = []
        for i, l in enumerate(self.leaves):
            ret += list(var[l,:])

        return np.array(ret)

    def blockBounds(self):
        ret = [self.bboxes[0,0,0]]

        # Loop over leaf blocks:
        for i,l in enumerate(self.leaves):
            ret.append(self.bboxes[l,0,1])
        return np.array(ret)
