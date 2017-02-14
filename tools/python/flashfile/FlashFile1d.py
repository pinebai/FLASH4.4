from FlashFile import FlashFile
import numpy as np

class FlashFile1d(FlashFile): 
    
    def __init__(self, fn):
        FlashFile.__init__(self,fn)

    def coord(self):
        """Get 1D array of cell center coordinates"""
        ret = []

        # Loop over leaf blocks:
        for i,l in enumerate(self.leaves):
            start = self.bboxes[l,0,0]
            size  = self.bboxes[l,0,1] - self.bboxes[l,0,0]
            dx = size/self.nb[0]
            ret += list(np.linspace(start + dx/2, start+size - dx/2, self.nb[0]))

        return np.array(ret)

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
