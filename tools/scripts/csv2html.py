#!/usr/bin/env python

import csv, sys

htmlHead = """
<HTML>
<HEAD>
<TITLE>Title of Page</TITLE>
<style type="text/css">
tr.evenrow {
   background-color: gray; 
   color: black;
}
tr.oddrow {
   background-color: white; 
   color: black;
}
</style>
</HEAD>
<BODY>
"""
htmlFoot = """
</BODY>
</HTML>
"""

fieldOrder = ["Name","Email","Webpage","SentInfo"]

def csv2html(csvfile=None,htmlfile=None):
    ifd = open(csvfile,"rb")
    reader = csv.DictReader(ifd,fieldnames=fieldOrder)
    anslist = []
    # generate the rows
    for row in reader:
        ans = []
        for f in fieldOrder:
            ans.append("    <TD>")
            if row.get(f,None):
               ans.append(row[f])
            else: 
               ans.append("NOEXIST")
            ans.append("</TD>\n")
        anslist.append("".join(ans))
    # Now write out the html
    ofd = open(htmlfile,"w")
    ofd.write(htmlHead)
    ofd.write("\n<TABLE RULES=GROUPS FRAME=BOX>\n")
    ofd.write("<THEAD>\n%s</THEAD>\n<TBODY>\n" % anslist[0])
    del anslist[0]
    currclass = "evenrow"
    for row in anslist:
        if currclass == "evenrow":
           currclass = "oddrow"
        else: currclass = "evenrow"
        ofd.write('<TR class="%s">\n%s</TR>\n' % (currclass,row))
    ofd.write("</TBODY>\n</TABLE>\n")
    ofd.write(htmlFoot)

if __name__=="__main__":
   if len(sys.argv) != 3:
      print "Usage: csv2html.py <csvfile> <htmlfile>"
      sys.exit(1)
   else:
      csv2html(csvfile=sys.argv[1],htmlfile=sys.argv[2])

