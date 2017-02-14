#!/usr/bin/env python

import sys, re, getopt

################ Global variables

comment = re.compile("^\s*!") # identify comment lines
quotes = re.compile("""((?:["][^"]*["])|(?:['][^']*['])|(?:[!].*)|(?:[&]\s*))""") # use to split a line to quoted strings, comments and others

############ Remove continuation lines in Fortran code to make free form
# generator mode: takes a generator and returns another
def rmContLinesF90(inlines):
  buffer = ""
  comment = "" # contains comment lines in between cont lines
  prevline = "" # The previous line

  lineno = 0
  for line in inlines:
    lineno += 1
    # remove spaces in current line
    scnt = len(line)-len(line.lstrip()) # how many initial spaces
    sline = line.strip() 
    if not buffer: # start of a new statement
       # flush comment 
       if comment: 
          yield comment
          comment = ""
       if sline.endswith("&"): # statement not complete
          buffer += scnt*" "+ sline[:-1].strip() # add current line to buffer
       else: # statement complete
          yield line
    else: # have an incomplete statement in buffer
       if sline.startswith("!"): # comment inside continuation
          comment += "\n%s" % line
       elif sline.startswith("&") or sline.startswith("."): # matching continuation line
          buffer += sline[1:].strip() # remove initial &
       else: # attach it with space
          buffer += " " + sline

       # if something left in buffer and needs to flushed out do it
       if buffer and not buffer.endswith("&"): 
          yield buffer
          if buffer[-1] != "\n": yield "\n"
          buffer = ""
       elif buffer: # ends with &
          buffer = buffer[:-1].rstrip()
  # end of for loop
  return # nothing more to yield

#### Replace . in first char with & in last of prev line

# Given a line, add a & at the end
def addCont(line):
    global quotes
    # preprocessor stuff
    if line[-1] == "\n": 
       rline = line[:-1] # remove trailing \n
    else: rline = line
    parts = [x for x in quotes.split(rline) if x] # all non-trivial parts of line
    if not parts: return line # we got an empty line
    # if last part is a comment put & before it else put it after
    if parts[-1][0] == "!": # The last part is a comment
       suffix = parts[-1]
       del parts[-1]
    else:
       suffix = ""
    # now comments are in suffix
    # check if we already have an &
    if parts and parts[-1][0] == "&": return line
    # account for removed \n
    ans = "%s & %s\n" % ("".join(parts),suffix)
    return ans

# convert F77 continuation to F90 continuation
def ContLines(inlines):
  prev = [] # cotains one statement (possibly split into multiple lines)
  comments = []
  for line in inlines:
      sline = line.lstrip()
      scnt = len(line)-len(sline)
      if sline.startswith("!"): # comment line
         comments.append(line)
         continue
      if line.startswith("c"): # comment char for F77
         comments.append("!"+line[1:])
         continue
      if sline.startswith(".") or sline.startswith("&") or sline.startswith("$"): # continuation char
         prev.append(" "*scnt +"&" +sline[1:])
         continue
      # now we have a regular line
      if prev:
         # go back in buffer till last non-pre processor line
         for x in prev[:-1]: # add & to all but last line
             yield addCont(x)
         if prev[-1].endswith("\n"):
            yield prev[-1] # now the last line
         else: yield prev[-1]+"\n"
         prev = []
      # put out the comments
      if comments:
         for x in comments: yield x
         comments = []
      # Handle current line
      prev.append(line)
  # we have used up input
  if prev:
     for x in prev[:-1]: yield addCont(x)
     if prev[-1].endswith("\n"):
        yield prev[-1] # now the last line
     else: yield prev[-1]+"\n"
     prev = []
  if comments:
     for x in comments: yield x
     comments = []
  return

# Add a reorder instruction
def addReorderLine(inlines):
   done = 0
   for line in inlines:
       if line.startswith("#include") and not done:
          done = 1
          yield '!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]\n'
          yield '!!REORDER(4): recvar[xyz]f\n'
       yield line
   return

############ house keeping code follows

def usage():
   print >> sys.stderr, "Usage: %s [options] " % sys.argv[0]
   print >> sys.stderr
   print >> sys.stderr, "Convert fixed format source code to free format"
   print >> sys.stderr
   print >> sys.stderr, "-i <filename>, --input=<filename>"
   print >> sys.stderr, "           Which file has input source."
   print >> sys.stderr, "           Use '-' for stdin"
   print >> sys.stderr 
   print >> sys.stderr, "-o <filename>, --output=<filename>"
   print >> sys.stderr, "           Which file should have output."
   print >> sys.stderr, "           Use '-' for stdout"
   print >> sys.stderr, "           Use '.EXT' for same name as input except extension is EXT"
   sys.exit(1)

def main():
    opts,rest = getopt.getopt(sys.argv[1:],"hi:o:",["help","input=","output="])
    input = "-"
    output = "-"
    for k,v in opts:
        if k in ["-i","--input"]:
           input = v
        elif k in ["-o","--output"]:
           output = v
        elif k in ["-h","--help"]:
           usage()
        else: 
           print "Unrecognized option pair (%s,%s)" % (k,v)
           usage()

    # Setup input and output file descriptors
    if input == "-":
       print >>sys.stderr, "Awaiting input from stdin"
       ifd = sys.stdin
    else: ifd = open(input,"r")
    if output == "-":
       ofd = sys.stdout
    elif output.startswith("."):
       output = input[:input.index(".")]+output
       ofd = open(output,"w")
    else: ofd = open(output,"w")

    # All the action happens here
    ofd.writelines(addReorderLine(ContLines(ifd))) # Here is the magic line
    ifd.close()
    ofd.close()


if __name__=="__main__": main()

