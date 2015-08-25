#!/usr/bin/python

scop_file = "/home/anuj/Downloads/dir.cla.scope.2.05-stable.txt"
dom_files = "dom_list"

## read the SCOP file
scop = {}
for line in open(scop_file, 'r'):
  if line[0] == '#':
    continue
  data = line.replace('\n','').split('\t')
  dom_name = data[1].lower()
  chain_name = data[2].lower().replace(':','')
  scop_family = data[3].split('.')[0]
  scop[(dom_name, chain_name)] = scop_family
#print scop

## read the domain files
lcount = 0; doms = {}
for line in open(dom_files, 'r'):
  lcount += 1
  if lcount == 1:
    continue
  doms[line.replace('\n','').replace('pdb','').replace('.ent','')] = 1

seen = {}
for k,v in scop.iteritems():
  if doms.get(k[0]) != None and seen.get(k[0]) == None:
    seen[k[0]] = 1
    print '%s\t%s' %(k[0], v)
