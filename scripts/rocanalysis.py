#!/usr/bin/env python

SCOP_FILE_NAME = 'dir.cla.scope.2.05-stable.txt' #'dir.cla.scop.txt_1.75'
#CE_LOG_FILE_NAME = 'pdb40d.psc.log.ce'
#TM_LOG_FILE_NAME = 'pdb40d.psc.log.tm'
#USM_LOG_FILE_NAME = 'pdb40d.psc.log.usm'
CE_LOG_FILE_NAME = 'apdb40.psc.log.ce'
TM_LOG_FILE_NAME = 'apdb40.psc.log.tm'
USM_LOG_FILE_NAME = 'apdb40.psc.log.usm'
DOM_LEN_FILE_NAME = 'dom.len.dat'
LEN_THRESH = 0

def scop_line_data(line):
	d = line.split('\t')
	return (d[0], (d[2], d[3]))

def get_data_ce_line(line):
	d = line.split(' ')
	cl = lambda x: x.replace('pdb','').split('.')[0]
	return ((cl(d[3]), cl(d[5])), float(d[11][:-2]))

def get_data_tm_line(line):
	d = line.split(' ')
	cl = lambda x: x.replace('pdb','').split('.')[0]
	return ((cl(d[3]), cl(d[5])), 1-float(d[13][:-1]))

def get_data_usm_line(line):
	d = line.split(' ')
	cl = lambda x: x.replace('pdb','').split('.')[0]
	return ((cl(d[3]), cl(d[5])), float(d[11][:-1]))

def is_of_class(k, sc, scop_map):
	return scop_map[k[0]][1][0] == sc and scop_map[k[1]][1][0] == sc

def normalize(data):
	vals = map(lambda x : x[1], data)
	minv = min(vals); maxv = max(vals); diffv = maxv-minv
	return map(lambda x : (x[0], (x[1] - minv)/diffv), data)
	
def is_long_enough(k, dom_len_map):
	return dom_len_map[k[0]] >= LEN_THRESH and dom_len_map[k[1]] >= LEN_THRESH

def fix_psc_data(data_map, sc, scop_map, dom_len_map):
	d = filter(lambda x: is_of_class(x[0], sc, scop_map), data_map.iteritems())
	#d2 = filter(lambda x: is_long_enough(x[0], dom_len_map), d)
	nd = normalize(d)
	return sorted(nd, key=lambda x: x[1])

def _get_val(k, dmap):
	v = dmap.get(k)
	if v == None:
		v = dmap.get((k[1],k[0]))
	return v

def make_mcpsc_data(ce_data,tm_data,usm_data):
	ret = {}
	ce_map = dict(ce_data); tm_map = dict(tm_data); usm_map = dict(usm_data)
	for k, v in ce_map.iteritems():
		ret[k] = sum([v,_get_val(k, tm_map),_get_val(k, usm_map)])/3
	return ret.iteritems()

def get_tpfp(data,scop_map):
	tp = 0; fp = 0
	for d in data:
		sf1 = scop_map[d[0][0]][1][2]
		sf2 = scop_map[d[0][1]][1][2]
	#	print [sf1,sf2]
		tp += sf1 == sf2
		fp += sf1 != sf2
	return (tp,fp)
	
def get_tpfp_curve(data,scop_map):
	l = len(data)
	c = map(lambda x : get_tpfp(data[0:int(x*l*1.0/100.0)],scop_map), range(0,101,10))
	P, N = c[-1]
	return map(lambda x : (x[0]*1.0/P,x[1]*1.0/N), c)

def calc_auc(roc):
	f_auc = lambda x : (roc[x][0]+roc[x-1][0]) * (roc[x][1]-roc[x-1][1]) / 2.0
	return sum(map(f_auc, range(1,len(roc))))

def dom_len_data(line):
	d = line.replace('\n','').split(' ')
	return (d[0],int(d[1]))

dom_len_map = {} #dict(dom_len_data(x) for x in open(DOM_LEN_FILE_NAME))

scop_map = dict(scop_line_data(x) for x in open(SCOP_FILE_NAME) if x[0] != '#')
##print scop_map

ce_map = dict(get_data_ce_line(x) for x in open(CE_LOG_FILE_NAME) if x[0:5] == 'Algo:')
##print ce_map

tm_map = dict(get_data_tm_line(x) for x in open(TM_LOG_FILE_NAME) if x[0:5] == 'Algo:')
##print tm_map

usm_map = dict(get_data_usm_line(x) for x in open(USM_LOG_FILE_NAME) if x[0:5] == 'Algo:')
##print usm_map

for SC in ['a','b','c','d']: 
	ce_data = fix_psc_data(ce_map, SC, scop_map,dom_len_map)
	print '%s:%d' %(SC,len(ce_data))
	ce_tpfp_data = get_tpfp_curve(ce_data,scop_map)
	print 'ce:%f' %calc_auc(ce_tpfp_data)
	tm_data = fix_psc_data(tm_map, SC, scop_map,dom_len_map)
	tm_tpfp_data = get_tpfp_curve(tm_data,scop_map)
	print 'tm:%f' %calc_auc(tm_tpfp_data)
	usm_data = fix_psc_data(usm_map, SC, scop_map,dom_len_map)
	usm_tpfp_data = get_tpfp_curve(usm_data,scop_map)
	print 'usm:%f' %calc_auc(usm_tpfp_data)
	mcpsc_data = sorted(make_mcpsc_data(ce_data,tm_data,usm_data), key=lambda x: x[1])	
	mcpsc_tpfp_data = get_tpfp_curve(mcpsc_data,scop_map)
	print 'mcpsc:%f' %calc_auc(mcpsc_tpfp_data)
