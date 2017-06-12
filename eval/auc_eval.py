#!/usr/bin/env python

import glob
import sys 


if len(sys.argv) != 3:
    print 'Usage: # python %s [file*.ev] 0or1' % sys.argv[0]
    quit() 
filename = sys.argv[1]
csv_files = filename#glob.glob(filename)

gcount=0
pos=0
neg=0
auc=0.0
acc=0.0
pnum=0
deval  = {}
dvalue = {}
#print "Read files :",
####for fname in csv_files:
for i in range(1):
    fname = filename
    #print fname,
    f = open(fname)
    line = f.readline()
    while line:
	itemlist = line[:-1].split('\t')
        if itemlist[0]=='1':
            pos=pos+1
            if float(itemlist[1]) >= 0:
                acc = acc + 1
        else:
            neg=neg+1
            if float(itemlist[1]) < 0:
                acc = acc + 1
        #print gcount,itemlist
        dvalue[gcount] = itemlist[0]
        deval[gcount] = itemlist[1]
        gcount = gcount + 1
        line = f.readline()
    f.close()
#print
#print dvalue

for k, v in sorted(deval.items(), key=lambda x:x[1], reverse=True):
    if dvalue[k]=='1':
        pnum = pnum+1
    else:
        auc = auc + pnum
if sys.argv[2] =='0':
    print acc/gcount
else:
    print auc/(pos*neg)

"""
print "All graph :",gcount
print "positive number :",pos
print "negative number :",neg
print "ACC       :",acc/gcount
print "AUC value :",auc/(pos*neg)
"""
