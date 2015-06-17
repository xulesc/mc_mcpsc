#!/usr/bin/python

INFILE = 'pdb40d_j.fa'
OUTFILE = 'dname2'

of = open(OUTFILE,'w')
for line in open(INFILE,'r'):
  if line[0] != '>':
    continue
  data = line.replace('>','').replace('\n','').split(' ')
#  print data
  of.write('%s\t%s\n' %(data[0][1:5], data[1][0]))
of.close()

