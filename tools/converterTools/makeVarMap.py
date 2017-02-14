#!/usr/bin/python

#will need, name of flash2 file to deal with,
#           name of target output file (defualt flash2vars, flash2parts)

import sys
from xml.sax import make_parser
from xml.sax import handler

#flash2VarNames = []

class unkNameHandler(handler.ContentHandler):

    inUnkNames = 0

    def startElement(self,name,attrs):
        if name == 'hdf5:Dataset':
            if attrs.getValue('Name') == 'unknown names' :
               self.inUnkNames = 1;

            
    def endElement(self, name): 
        if name == 'hdf5:Dataset':
            self.inUnkNames = 0;
            
    def characters(self, data):
        if self.inUnkNames == 1:
            if not data.isspace() :
                print data.strip().strip('"').encode('ascii')
                

class SimpleDTDHandler(handler.DTDHandler):
    def notationDecl(self,name,publicid,systemid):
        print "Notation: ", name,publicid, systemid
    def unparsedEntityDecl(self,name,publicid,systemid,ndata):
        print "unparsedEntity: ", name,publicid,systemid, ndata


#class handleUnkNames(handler.ContentHandler):
#    def startElement(self, name,attrs):
#        if name == 'hdf5:Dataset' :
#            if attrs.getValue('Name') == "unknown names" 

p=make_parser()
p.setContentHandler(unkNameHandler())
p.setDTDHandler(SimpleDTDHandler())

p.parse(sys.stdin)
#print flash2VarNames



