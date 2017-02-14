__doc__ = "This module is to be used to pre-process text files"
__all__ = ["preProcess"]

import re
import globals
from globals import * # for GVars
   
class preProcess:

   rif = re.compile(r"^\s*IF\s*(?P<cond>.*)\s*$")
   relif = re.compile(r"^\s*(?:ELSEIF|ELIF)\s*(?P<cond>.*)\s*$")
   ranyif = re.compile(r"^\s*(?:IF|ELSEIF|ELIF)\s*(?P<cond>.*)\s*$")
   rerr = re.compile(r"^\s*SETUPERROR\s*(?P<msg>.*)\s*$")
   relse  = re.compile(r"^\s*ELSE\s*(?:#.*)?$")
   rendif = re.compile(r"^\s*END\s*IF\s*(?:#.*)?$")
   ruse   = re.compile(r"^\s*USESETUPVARS\s+(?P<vars>.*)(?:#.*)?$")

   rD     = re.compile(r"^\s*D\s+")
   rParam = re.compile(r"^\s*PARAMETER\s+")

   normalkeywords = ["SETUPERROR"] # keywords to be ignored inside a block
   
   def __init__(self,values={},ignorePP=False):
      """values=Dictionary mapping variable to value used to evaluate
      conditionals"""
      self.initvalues = {}
      self.initvalues.update(values)
      self.values = {}
      self.filename = ""
      self.lineno = 0
      self.ignorePP = ignorePP
   
   def processFile(self,filename):
      fd = open(filename,"r")
      self.filename = filename
      # if the first line of the file starts with ##python, then we compile it as python and run its genLines method to get the lines
      if fd.readline().startswith("##python:genLines"):
          # read the rest of the file (doesnt include ##python:genLines line) and compile it into code
          code = compile(fd.read(), filename, "exec")
          fd.close()
          # initialize a namespace with the prelude code
          ns = __import__("genLinesPrelude")
          # execute the code object within the prelude'd namespace
          exec code in ns.__dict__
          # helper to split multi-line strings and flatten nested iterables into one long iterator of strings
          def flatten(it):
              if type(it) == type(""):
                  for x in it.split('\n'):
                      yield x
              else:
                  for x in it:
                      for y in flatten(x):
                          yield y
          # def guarddict(d):
          #     def baduser(msg):
          #         e = Exception()
          #         e.dumbuser = True
          #         e.message = msg
          #         raise e
          #     class A(dict):
          #         def __init__(self,d):
          #             self.d = d
          #             print d
          #             dict.__init__(self)
          #         def __getattribute__(self,x):
          #             return getattr(self.d,x)
          #         def __getitem__(self,x):
          #             if self.d.has_key(x):
          #                 return self.d[x]
          #             else:
          #                 baduser('The setup variable "%s" is required by the file "%s".' % (x, filename))
          #         def __setitem__(self,x,y):
          #             raise Exception('Config scripts may not write to the setup variable dictionary.')
          #     return A(d)
          # # execute the genLines function to get a back an iterable of lines
          # lines = flatten(ns.genLines(guarddict(self.initvalues.copy())))
          lines = flatten(ns.genLines(self.initvalues.copy()))
      else: # just a regular config file
          fd.seek(0)
          lines = fd.readlines()
          fd.close()
      # exception goes to caller
      if self.ignorePP:
          # In this case, this method has been called expressly for
          # the purpose of generating documentation, not because we
          # are really setting up a problem, so we use "filterLines"
          # to strip everything from the Config file but parameter
          # names and values (prefixed by the keyword "PARAMETER")
          # and parameter descriptions (prefixed by a capital 'D')
          ans = self.filterLines(lines)
      else:
          ans = self.processLines(lines)
      return ans

   def printVars(self):
       for k,v in self.values.items():
           if k.startswith("with"): continue
           if type(v)==type(""):
              print "%s='%s'" % (k,v)
           else: print '%s=%s' % (k,v)

   def PPKeyword(self,line):
      """Returns "IF,ELSE,ENDIF,USE,None" depending on syntax of line"""
      if self.rif.match(line):
         return "IF"
      if self.relif.match(line):
         return "ELIF"
      elif self.relse.match(line):
         return "ELSE"
      elif self.rendif.match(line):
         return "ENDIF"
      elif self.ruse.match(line):
         return "USESETUPVARS"
      elif self.rerr.match(line):
         return "SETUPERROR"
      else: return None

   def docKeyword(self, line):
      if self.rD.match(line):
         return True
      if self.rParam.match(line):
         return True
      return False
   
   def raiseError(self,line):
      """Print user defined error message"""
      m = self.rerr.match(line)
      if not m: return
      GVars.out.put("Setup halted by rule defined in %s (line %d)" % (self.filename, self.lineno), globals.ERROR)
      raise SetupError(m.group("msg"))

   def checkUse(self,line):
      """Confirm that variables named in line are present. Else raise Warning"""
      m = self.ruse.match(line)
      if not m: return
      vars = [x.strip() for x in m.group("vars").split(",") ]
      badlist = [x for x in vars if not self.values.has_key(x) ]
      # all declared variables (except for the with* ones)
      initvars = [x for x in self.values.keys() if not x.startswith("with")]
      if badlist: # found missing variables
         msg = ["\nPre-Processor Warning: \nRequired uninitialized variables will be initialized to ''"]
         msg.append("File: %s\nUnnitialized Variables: %s" % (self.filename,", ".join(badlist)))
         msg.append("Initialized Variables: %s" % ", ".join(initvars))
         msg.append("NOTE: * Variable Names are case sensitive. ")
         msg.append("      * The withUNIT variables are not shown in the initialized list.")
         GVars.out.push()
         GVars.out.put("\n".join(msg),globals.PPWARN)
         GVars.out.pop()
         for x in badlist: self.values[x] = ""
          
   def handlePPCmd(self,m,line):
       if m=="SETUPERROR": # need to raise an error
          self.raiseError(line)
       else:
          raise SetupError("Programming error. Got unhandled pre-processing command" % m)
        
   def evalCondn(self,line):
      """Evaluates condition and returns True/False/None. 
      returns None if line has bad syntax"""
      m = self.ranyif.match(line) # either IFDEF or ELSEIF
      if not m: return None
      cond = m.group("cond")
      valcopy = {}
      valcopy.update(self.values)
      try:
         value = eval(cond,{},valcopy)
      except Exception,e: 
         msg = 'File: %s\nLine: %d\nCondition: %s\nVariables: %s\nDetails: %s'
         raise SetupError(msg % (self.filename,self.lineno,cond,self.values,str(e)))
      else:
         if valcopy != self.values: # Warn user about side effect
            msg = ["PreProcessor SIDEEFFECT Warning: Condition has side effect"]
            msg.append("File: %s\nLine : %d\nCode: %s\nBefore: %s\nAfter: %s"
                       % (self.filename,self.lineno,cond,self.values,valcopy))
            msg.append("Trashing Side effects")
            GVars.out.push()
            GVars.out.put("\n".join(msg),globals.WARN)
            GVars.out.pop()
      return value
   

   def filterLines(self,lines):
      # Remove everything except for lines starting with 'D '
      # (param descriptions) and the keyword 'PARAMETER'
      ans = []
      self.lineno = 0
      for line in lines:
         self.lineno += 1
         if self.docKeyword(line):
             ans.append((self.lineno, line))
      return ans

   # The magic code here is explained here. We evaluate arbitrarily nested 
   # conditionals without using any recursion!!! (only a stack of pairs of truth values)
   # * start by pushing True to the stack
   # * normal lines go through iff top of Stack is True
   # * An IF line is evaluated and truth value pushed on stack
   # * An ENDIF pop's the truth value off the stack
   # * An ELSE negates top of stack
   # all other variables are just for good error messages
   #
   # Seems to be a cool way to do it. 
   #
   # To handle ELIF statements we need to push a pair everytime
   #    (Found a true block so far?, truthvalue of current block)
   #   --- GMK
   def PreProcessLines(self,lines):
      stack = [(True,True)]
      iflines = []
      ans = []
      self.lineno = 0
      nontrivial = 0
      self.values.update(self.initvalues)
      if not self.filename: self.filename = "<String List>"
      for line in lines:
         self.lineno += 1
         state = stack[-1][1] # current truth value part of top of stack
         m = self.PPKeyword(line)
         # common case first (regular line)
         if (not m):
            if state: # Replace setupvariables with values
               ans.append((self.lineno, line))
            continue
         if m in self.normalkeywords: # Got a PPCmd to handle (not conditional)
            if state: # not to be ignored
               self.handlePPCmd(m,line)
            continue
         # Now m is one of the conditional commands
         nontrivial = 1
         if m == "IF": # need to evaluate an if
            rv = self.evalCondn(line) # exception passed through to caller
            stack.append((rv,rv)) # push truth value pair to stack
            iflines.append(self.lineno) # note location of IF
         elif m == "ELSE": # switch the most recent state
            a,b = stack[-1]
            b = not a
            stack[-1] = (a,b) 
         elif m == 'ELIF': # need to evaluate condition
            a,b = stack[-1]
            if not a: # if not already found a true block 
               rv = self.evalCondn(line) # exception passed through to caller
               a,b = rv,rv
            else: # ignore this block
               a,b= True,False
            stack[-1] = a,b
         elif m == "ENDIF": # pop the most recent state
            stack.pop()
            iflines.pop()
            if not stack: # stack has become empty
               msg = ["PreProcessor Error: ENDIF found without a corresponding IF"]
               msg.append("Extra ENDIF found at line %d in %s" % (self.lineno,self.filename))
               msg.append("The real error could be located before the specified line also")
               raise SetupError("\n".join(msg))
         elif m=="USESETUPVARS": # found a use line
            self.checkUse(line)
         else: # regular line
            raise SetupError("Programming Error")
      # finished processing all lines
      if len(stack) != 1: # we have an unbalanced stack
         msg = ["PreProcessor Error: IF found without a corresponding ENDIF"]
         msg.append("IF at line %d in %s does not have a matching ENDIF" % (iflines[-1],self.filename))
         raise SetupError("\n".join(msg))
      if nontrivial: # print the result on stdout
         GVars.out.put("Processed Version of %s" % self.filename,globals.PPDEBUG)
         GVars.out.push()
         for lno,line in ans: GVars.out.put(line.replace("\n",""),globals.PPDEBUG)
         GVars.out.pop()
      return ans

   #Another name of PreProcessLines
   processLines = PreProcessLines
