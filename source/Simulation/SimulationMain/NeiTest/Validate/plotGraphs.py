#! /usr/bin/env python

"""plotGraphs.py -- Creates several graphs for the Ionize unit test,
which we use to validate the test."""

import Gnuplot, Gnuplot.funcutils, glob, re, sys


def CreateSingleGraph(listDataFiles):

    """Reads data from ASCII space delimited files, and plots
    graphs using the gnuplot package."""

    #We add an extra step.
    #The graphs in the FLASH2 user-guide show a population fraction for
    #each element normalised to 1.  We will grab the population
    #fraction data for each ion at x=0 (way before the temperature step),
    #sum the total, and generate a factor to normalise the population fraction
    #within each element to 1.
    sumTotal = 0.0
    for fileName in listDataFiles:
        fileObj = open(fileName, 'r')
        firstLine = fileObj.readline()
        (x,y) = firstLine.split()
        sumTotal += float(y)
        fileObj.close()
        
    normaliseFactor = 1.0 / sumTotal


    #Now we have the normalisation factor, we can plot the required graphs.
    g = Gnuplot.Gnuplot(debug=1)

    elementStr = "%s" % (listDataFiles[0][0:2])
    elementStr = elementStr.strip()
    
    titleStr = "Population of %s ions" % (elementStr)
    psFileStr = "%s.ps" % (elementStr)

    g.title(titleStr)
    g('set key outside')
    g('set xlabel "x(cm)"')
    g('set ylabel "Population Fraction"')
    g('set log xy')
    g('set xrange [8E6:1E8]')
    g('set yrange [7E-6:1]')
    g('set mxtics 10')
    g('set mytics 10')
    g('set format x "10^{%L}"')
    g('set format y "10^{%L}"')

    
    for fileName in listDataFiles:
        fileObj = open(fileName, 'r')
        listFileLines = fileObj.readlines()

        listTuples=[]
        for line in listFileLines:
            (x,y) = line.split()
            listTuples.append((float(x),normaliseFactor*float(y)))
        fileObj.close()

        #This is a hack that seems to work.  Just keep replotting until
        #we have included all Gnuplot.Data objects.
        d = Gnuplot.Data(listTuples, title=fileName, with_='lines lw 5')
        g.replot(d)

    g.hardcopy(psFileStr, fontsize=20, enhanced=1, color=1)
    print ('\n******** Saved plot to postscript file "%s" ********\n' % psFileStr)



if __name__ == '__main__':

    listDatFiles = glob.glob('*.dat')
    listElements = ['he','c ','n ','o ','ne','mg','si','s ','ar','ca','fe','ni']

    #We wish to plot a graph for each element.
    for element in listElements:

        #Pattern match the .dat files for this element.
        elementREstring = element
        elementREstring += "\d+.dat"
        elementRE = re.compile(elementREstring, re.IGNORECASE)

        listValidNames = []
        for fileName in listDatFiles:            
            m = elementRE.match(fileName)
            if m:
                listValidNames.append(m.group())

        listValidNames.sort()

        #Plot every single item in list on a single graph.
        CreateSingleGraph(listValidNames)
        #sys.exit(0)
