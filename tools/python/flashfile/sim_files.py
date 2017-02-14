import os

def simFiles(basename, ftype="chk", directory='.'):
    files = os.listdir(directory)
    files.sort()
    
    ret = []
    for fn in files:
        if fn.startswith(basename + "_hdf5_" + ftype + "_"):
            ret.append(directory + '/' + fn)

    return ret
