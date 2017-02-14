#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import time
import array
#import modeldata

print "************ Python Script started ************"

starttime = time.localtime()
year, month, day, hour, minute, second = starttime[0:6]
print "Start time:    %4i-%02i-%02i  %02i:%02i.%02i" % (year, month, day, hour, minute, second)


def read_table(file,table,table_str,readheader,tofloat):
	f = open(file, 'r')
	header = ""
	if readheader: header = f.readline()
	for line in f:
		lstr = line.split()
		larr = []
		append_line = True
		for i in range(0,len(lstr)):
			try:
				if tofloat:
					larr.append(float(lstr[i]))
				else:
					larr.append(lstr[i])
			except ValueError:
				append_line = False
				break
		if append_line:
			table.append(larr)
			table_str.append(line)
	print 'columns in table : '+str(len(table[0]))
	print 'lines in table   : '+str(len(table))
	f.close()
	return header

def write_table(file,table_str,header):
	nl = len(table_str)
	f = open(file, 'w')
	f.write(header)
	for i in range(0,nl):
		f.write(table_str[i])
	f.close()
	print "'"+file+"' written. "+str(nl)+" lines."

def find_line(table,col,val):
	nl = len(table)
	for i in range(0,nl):
		if table[i][col] < val:
			continue
		else:
			if i > 0:
				if (table[i][col]-val) < (val-table[i-1][col]):
					return i
				else:
					return i-1
			else:
				return -1

def find_first(table,col,val):
        nl = len(table)
        for i in range(0,nl):
                if table[i][col] < val:
                        continue
                else:
			return i

def tofloat(tab):
	for i in range(0,len(tab)):
		try:
			tab[i] = float(tab[i])
		except ValueError:
			print tab[i]
			break
	return tab

def clean_recursive(table,table_str,col):
	nl = len(table)
	for i in range(0,nl-1):
		if table[i+1][col] <= table[i][col]:
			iend = i+1
			istart = find_first(table,col,table[iend][col])
			print "removing lines: ", istart, iend, table[istart][col], table[iend-1][col], table[iend][col]
			for j in range(istart,iend):
				# print 'popping line: '+str(table[istart][col])
				table.pop(istart)
				table_str.pop(istart)
			clean_recursive(table,table_str,col)
			break

def copy_file(infile,outfile):
        shellcmd = "cp "+infile+" "+outfile
        print shellcmd
        status = os.system(shellcmd)

#models = []
#models_str = []
#models_header = modeldata.read("models.txt",models,models_str,True)
#modelname = ""
#for l in range(0,len(models)):
#	if modelname != models[l][0]: # new model
#		modelname = models[l][0]
#                print "CLEANING model: "+modelname
filename = "SMT.dat"
outfile = filename+'_cleaned'
copy_file(filename,filename+'_sav')
table = []
table_str = []
header = read_table(filename,table,table_str,True,True)
clean_recursive(table,table_str,0)
write_table(outfile,table_str,header)


print "************ Python Pipeline finished ************"
print "Start time:    %4i-%02i-%02i  %02i:%02i.%02i" % (year, month, day, hour, minute, second)
endtime = time.localtime()
year, month, day, hour, minute, second = endtime[0:6]
print "Endtime time:  %4i-%02i-%02i  %02i:%02i.%02i" % (year, month, day, hour, minute, second)
start_time = time.mktime(starttime)
end_time   = time.mktime(endtime)
operation_time_in_seconds = end_time - start_time
operation_time_in_hours   = operation_time_in_seconds/3600
print "Operation time in seconds: " + str(operation_time_in_seconds) + "  seconds"
print "Operation time in hours  : " + str(operation_time_in_hours)   + "  hours"
print "************ Python Script finished ************"

