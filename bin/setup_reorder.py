#!/usr/bin/env python

import sys, re, getopt, os

################ Global variables

# Dictionary mapping DIMENSION to list of regexps matching name of array
# add to this array to handle more variables

ArrayNames = { "FIVE": [],
               "FOUR": [],
               "THREE": [],
               "TWO": []
             }

regexp = None # Regexp object which identifies the variables
comment = re.compile("^\s*!") # identify comment lines
#quotes = re.compile("""((?:["][^"]*["])|(?:['][^']*['])|(?:[!].*)|(?:[&]\s*))""") 
# use to split a line to quoted strings, comments and others
# I commented this out b/c Murali doensn't seem to use it anywhere - nttaylor

reordRE = re.compile("^\s*!!\s*REORDER[(](?P<num>[2|3|4|5])[)]:\s*(?P<list>.*)$", re.I)

######### Compute variables to be reordered

# Given a filename gets the list of variables to be reordered
# returns info in dictionary with keys "FOUR","FIVE","FLAGS"
def getREORDdict(filename):
    global reordRE
    """Handle one file"""
    reorddata = {"TWO":[], "THREE":[], "FOUR" : [],"FIVE":[],"FLAGS":{} }
    if not os.path.isfile(filename): return reorddata
    for x in file(filename).readlines():
        m = reordRE.match(x)
        if m:
            if m.group("num") == "2":
                key = "TWO"
            elif m.group("num") == "3":
                key = "THREE"
            elif m.group("num") == "4":
                key = "FOUR"
            else: key="FIVE"
            vnames = [x.strip() for x in m.group("list").split(",")]
            for v in vnames:
                if v.startswith("."): 
                    reorddata["FLAGS"][v[1:]] = 1 # handle flags
                    continue
                reorddata[key].append(v)
    
    return reorddata

############ Rewrite Array Access

def makeregexp():
    global ArrayNames,regexp

    rlist = []
    for suf,lst in ArrayNames.items():
        if lst:
           rlist.append(r"(?P<%s>%s)" % (suf,"|".join([r"(?:%s)"%x for x in lst])))
    nameregexp = r"(?P<name>%s)\s*[(]" % "|".join(rlist) 

    # make the regular expression
    regexp = re.compile(nameregexp)

def replfunc(mobj):
    d = mobj.groupdict()
    s,e = mobj.span("name")
    ws,we = mobj.span() # whole match
    pre,post ="",""
    # find char before and after name
    if s>0: pre = mobj.string[s-1]
    if e<len(mobj.string): post = mobj.string[e]
    # if surrounded on either end by alphanumeric or _ return unchanged
    if pre.isalnum() or post.isalnum() or pre == "_" or post == "_":
       return mobj.string[ws:we]
    ppname = ""
    for x in ["TWO","THREE","FOUR","FIVE"]:
        if d.get(x,None): ppname = "ARRAY%s" % x
    return "%s(%s," % (ppname,d["name"])

# Given one line replace array access
# also changes in fortran comments as well as inside strings
# also works in generator mode
def Array(inlines):
    global regexp,comment

    makeregexp()

    for line in inlines:
        if comment.match(line):
           yield line
        else: 
           yield regexp.sub(replfunc,line)
    return

############ house keeping code follows

def usage():
   print >> sys.stderr, "Usage: %s [options] " % sys.argv[0]
   print >> sys.stderr
   print >> sys.stderr, "Processes input freeform source code and reorder access to certain arrays"
   print >> sys.stderr
   print >> sys.stderr, "-i <filename>, --input=<filename>"
   print >> sys.stderr, "           Which file has input source."
   print >> sys.stderr, "           Use '-' for stdin"
   print >> sys.stderr 
   print >> sys.stderr, "-o <filename>, --output=<filename>"
   print >> sys.stderr, "           Which file should have output."
   print >> sys.stderr, "           Use '-' for stdout"
   print >> sys.stderr, "           Use '.EXT' for same name as input except extension is EXT"
   print >> sys.stderr 
   print >> sys.stderr, "--five=<arrname>, --four=<arrname>"
   print >> sys.stderr, "           Add specified name to list of array names to process"
   print >> sys.stderr, "           five implies this is a 5d array, four means 4d array"
   print >> sys.stderr, "           Usually 4d=space+varIndex, 5d=4d+BlockIndex"
   print >> sys.stderr 
   print >> sys.stderr, "--auto"
   print >> sys.stderr, "           Find out the names from the !!REORDER instructions in the file"
   sys.exit(1)

def main():
    global ArrayNames
    opts,rest = getopt.getopt(sys.argv[1:],"ha:i:o:",["help","input=","output=","five=","four=","auto"])
    input = "-"
    output = "-"
    auto = False
    for k,v in opts:
        if k in ["-i","--input"]:
           input = v
        elif k in ["-o","--output"]:
           output = v
        elif k in ["--five","--four", "--three", "--two"]:
           ArrayNames[k[2:].upper()].append(v)
        elif k in ["--auto","-a"]:
           auto = True
        elif k in ["-h","--help"]:
           usage()
        else: 
           print >> sys.stderr, "Unrecognized option pair (%s,%s)" % (k,v)
           usage()

    # now for the magic stuff
    if input == "-":
       if auto: print >>sys.stderr, "auto does not work with stdin. Ignoring"
       print >>sys.stderr, "Awaiting input from stdin"
       ifd = sys.stdin
    else: ifd = open(input,"r")
    if output == "-":
       ofd = sys.stdout
    elif output.startswith("."):
       output = input[:input.index(".")]+output
       ofd = open(output,"w")
    else: ofd = open(output,"w")


    # Compute additional Array names if requested
    if auto:
       reordnames = getREORDdict(input)
       if "ERROR" in reordnames["FLAGS"]:
          print >>sys.stderr, "File has an ERROR flag in the !!REORDER line. Quitting"
          sys.exit(1)
       ArrayNames["TWO"].extend(reordnames["TWO"])
       ArrayNames["THREE"].extend(reordnames["THREE"])
       ArrayNames["FOUR"].extend(reordnames["FOUR"])
       ArrayNames["FIVE"].extend(reordnames["FIVE"])

       

    # Get the list of output lines
    outputlines = list(Array(ifd)) 
    ifd.close()

    #PASS 1:
    # Do the macro substitution ourselves
    reg = re.compile("ARRAY(?:TWO|THREE|FOUR|FIVE)[(]")
    hugeline = "".join(outputlines) # file = 1 long string

    start = 0
    changes = [] # triples (x,y,slots) replace [x:y] part of string with info from slots
    # search for macro occurence
    m = reg.search(hugeline,start)
    while m:
      s = m.start()
      e= m.end()
      count = 1
      stack = [-1]
      slots = []
      slotstart = e
      ptr = slotstart
      
      # parse each char keeping track of how deep in ( nesting we are
      # and if we find a , while at outer level of ( nesting 
      # then we have one arg of macro and save it into slots
      while count > 0:
            c = hugeline[ptr]
            if c == "(": count += 1
            if c == ")": count -= 1
            if (c == "," and  count == 1) or (c==')' and count==0):
               slots.append(hugeline[slotstart:ptr])
               slotstart = ptr+1
            ptr += 1
      
      changes.append((s,ptr,slots))
      # start searching from where we finished
      start = ptr
      m = reg.search(hugeline,start)
    #END PASS 1

    #BEGIN PASS 2
    #do another pass to get dimension and pointer declarations.
    dimensionReg = re.compile("dimension[(]",re.I)
    start = 0
    m = dimensionReg.search(hugeline,start)
    ArrayNamesKeys = ArrayNames.keys()
    while m:  
      s = m.start()
      e = m.end()
      count = 1
      stack = [-1]
      slots = []
      slotstart = e
      ptr = slotstart
      slots.append("dimension")
      #gather up the dimension definition
      while count > 0:
            c = hugeline[ptr]
            if c == "(": count += 1
            if c == ")": count -= 1
            if (c == "," and  count == 1) or (c==')' and count==0):
               slots.append(hugeline[slotstart:ptr])
               slotstart = ptr+1
            ptr += 1
      
      
      #locate that :: so we can figure out if we have an array or not
      count = 2
      colonPtr = ptr
      while count > 0:
         c = hugeline[colonPtr]
         if c == ":": count -= 1
         colonPtr += 1
      

        
      #from here find the first alpha character.  This will be the actual
      #name we're concerned with.
      count = 0
      reordering = False
      while count == 0:
        c = hugeline[colonPtr]
        if c.isalpha() : 
           for key in ArrayNamesKeys:
             for var in ArrayNames[key]:
               targetReg = re.compile("%s\\b" % var)
               if targetReg.match(hugeline,colonPtr):
                 changes.append((s,ptr,slots))
                 reordering = True
                 break
             if reordering : break
             #end iteration over arraynames and keys
           count += 1                        
        colonPtr += 1 #end while count == 0
      
        
      #Check to make sure that all the arrays in this declaration are
      #supposed to be reordered.  If any of them aren't return an error.
      endFound = False
      while (not endFound) and reordering:
        ampersandFound = False
        newlineFound = False
        commaFound = False
        while not (newlineFound or commaFound):
          c = hugeline[colonPtr]  
          if c == "," : commaFound = True
          if c == "&" : ampersandFound = True
          if c == "\n" : 
              newlineFound = True
          colonPtr += 1
          
        if commaFound:
          count = 0
          foundMatch = False
          while count == 0:
            c = hugeline[colonPtr]
            if c.isalpha() :
              for key in ArrayNamesKeys:
                for var in ArrayNames[key]:
                  targetReg = re.compile("%s\\b" % var)
                  if targetReg.match(hugeline,colonPtr):
                    foundMatch = True
                    break
                if foundMatch : break
               #end ArrayNames and keys iterations
              count += 1                        
            colonPtr += 1 #end while count == 0
          if not foundMatch:
              print "ERROR: A variable specified in a dimension declaration that is to be reordered, is not declared in the REORDER declaration of the file %s, and is not specified on the command line!  Aborting!" % input
              sys.exit(-1)

        if newlineFound and (not ampersandFound): 
          endFound = True

      #Go onto the next Fortran line.
      start = ptr
      m = dimensionReg.search(hugeline,start)
    #END PASS 2

    changes = sorted(changes)
    # now write out the whole thing
    curr = 0
    for (s,e,what) in changes:
        ofd.write(hugeline[curr:s])
        aname = what[0]
        changed = what[2:5]
        changed.append(what[1])
        changed.extend(what[5:])
        #what = (arrayname, var,x,y,z...)
        #changed = [x,y,z,varname...]
        # Writeout the re-ordered code
        ofd.write("%s(%s)" % (aname,",".join(changed)))
        curr = e
    # write out last segment of file
    ofd.write(hugeline[curr:])
    ofd.close()

if __name__=="__main__": main()


