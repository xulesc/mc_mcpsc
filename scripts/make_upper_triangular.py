#!/usr/bin/python

import sys, os

def process(dirname, flistFile, jlistFile):
  files = os.listdir(dirname)
  pairs = []
  for x in range(0, len(files)):
    for y  in range(x, len(files)):
      pairs.append('%d %d' %(x, y))
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

