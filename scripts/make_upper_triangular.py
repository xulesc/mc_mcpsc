#!/usr/bin/python

import sys, os

def process(dirname, flistFile, jlistFile):
  files = sorted(os.listdir(dirname))
  pairs = []
  for x in range(0, len(files)):
    for y  in range(x + 1, len(files)):
      pairs.append('%d %d' %(x, y))
  flistFile.write('%d\n' %len(files))      
  flistFile.write('%s\n' %'\n'.join(files))
  jlistFile.write('%d\n' %len(pairs))
  jlistFile.write('%s\n' %'\n'.join(pairs))


def process2(pdb40dfname,flistFile,jlistFile):
  files = []; pairs = []
  afiles = os.listdir('/home/anuj/Downloads/astral_40p/')
  afilesm = {}
  for f in afiles: afilesm[f] = 1
  for line in open(pdb40dfname,'r'):
    if line[0] != '>':
      continue
    fname = line.split(' ')[0].replace('>','')
    if afilesm.get('%s.ent' %fname) == None:
      continue
    files.append('%s.ent' %fname)
  for x in range(0,len(files)):
    for y in range(x+1, len(files)):
      pairs.append('%d %d' %(x,y))
  flistFile.write('%d\n' %len(files))      
  flistFile.write('%s' %'\n'.join(files))
  jlistFile.write('%d\n' %len(pairs))
  jlistFile.write('%s' %'\n'.join(pairs))
  

if __name__ == '__main__':
  dirname = sys.argv[1]
  flistfname = sys.argv[2]
  jlistfname = sys.argv[3]
  
  print [dirname, flistfname, jlistfname]
  process(dirname, open(flistfname, 'w'), open(jlistfname, 'w'))

