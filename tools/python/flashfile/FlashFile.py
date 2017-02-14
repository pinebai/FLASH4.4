import tables
import numpy as np

class FlashFile:
    def __init__(self, fn):
        self.h5file = tables.openFile(fn, "r")

        # Load block bounding boxes:
        self.bboxes = self.h5file.getNode("/bounding box")[:,:]

        # Load block coordinates:
        self.coords = self.h5file.root.coordinates[:,:]

        # Generate a list of leaf blocks:
        block_types = self.h5file.getNode("/node type")[:]

        # The leaf nodes have type 1:
        self.leaves = np.where(block_types == 1)[0]

        # Load number of cells in each block:
        self.nb = ( self.integerScalar("nxb"),
                    self.integerScalar("nyb"),
                    self.integerScalar("nzb"), )

        self.time = self.realScalar("time")
        self.cycle = self.integerScalar("nstep")

    def integerScalar(self, name, rtp=False):
        n = "%-80s" % name
        t = self.h5file.getNode("/integer scalars")
        if rtp == True:
            t = self.h5file.getNode("/integer runtime parameters")

        wl = t.getWhereList("name == n")
        if len(wl) != 1:
            raise ValueError("Could not find integer \"%s\"" % name)
        return t[wl]["value"]


    def stringScalar(self, name, rtp=False):
        n = "%-80s" % name
        t = self.h5file.getNode("/string scalars")
        if rtp == True:
            t = self.h5file.getNode("/string runtime parameters")

        wl = t.getWhereList("name == n")
        if len(wl) != 1:
            raise ValueError("Could not find string \"%s\"" % name)
        ret = t[wl]["value"]
        if hasattr(ret, '__len__'): ret = ret[0]
        return ret


    def realScalar(self, name, rtp=False):
        n = "%-80s" % name
        t = self.h5file.getNode("/real scalars")
        if rtp == True:
            t = self.h5file.getNode("/real runtime parameters")

        wl = t.getWhereList("name == n")
        if len(wl) != 1:
            raise ValueError("Could not find real number \"%s\"" % name)
        return t[wl]["value"]

    def varnames(self):
        node = self.h5file.getNode("/unknown names")
        return node[:,0]

    def exists(self, name):
        return name in self.h5file.getNode("/")

    def __del__(self):
        self.h5file.close()

    def var(self, name):
        node = self.h5file.getNode("/" + name)
        return node
