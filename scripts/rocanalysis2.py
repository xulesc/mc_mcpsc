#!/usr/bin/env python

APDB40_TO_ASTRAL_FILE_NAME = 'apdb40_astral_mapping.clean'
ADVTOPS_DATASET_FILE_NAME = 'advTOPSplus-PDB40_subset-results.txt'
CE_LOG_FILE_NAME = 'apdb40.psc.log.ce'
TM_LOG_FILE_NAME = 'apdb40.psc.log.tm'
USM_LOG_FILE_NAME = 'apdb40.psc.log.usm'

def get_data_ce_line(line):
        d = line.split(' ')
        cl = lambda x: x.split('.')[0]
        return ((cl(d[3]), cl(d[5])), float(d[11][:-2]))

def get_data_tm_line(line):
        d = line.split(' ')
        cl = lambda x: x.split('.')[0]
        return ((cl(d[3]), cl(d[5])), 1-float(d[13][:-1]))

def get_data_usm_line(line):
        d = line.split(' ')
        cl = lambda x: x.split('.')[0]
        return ((cl(d[3]), cl(d[5])), float(d[11][:-1]))
        
def get_astral_to_apdb_line(line):
        dom, fname = line.split(' ')
        return (fname[fname.rfind('/')+1:fname.rfind('.')],dom)

def normalize(data):
	vals = map(lambda x : float(x[1]), data.iteritems())
	minv = min(vals); maxv = max(vals); diffv = maxv-minv
	return dict(map(lambda x : (x[0], (float(x[1]) - minv)/diffv), data.iteritems()))

def filter1(d1,d2,astral2pdb40,klassification,klass):
	ad1 = astral2pdb40[d1]; ad2 = astral2pdb40[d2]
	v = klassification.get((ad1,ad2))
	if v == None:
		v = klassification.get((ad2,ad1))
	if v == None:
		return False
	return v[0].split('_')[1][0] == klass and v[1].split('_')[1][0] == klass

def is_match(v1,v2):
	return v1.split('_')[1][2] == v2.split('_')[1][2]

def get_tpfp(psc,astral2pdb40,klassification):
	tp = 0; fp = 0
	for (d1,d2),s in psc:
		ad1 = astral2pdb40[d1]; ad2 = astral2pdb40[d2]
		v = klassification.get((ad1,ad2))
		if v == None:
			v = klassification[(ad2,ad1)]
		v1,v2 = v
		m = is_match(v1,v2)
		tp += is_match(v1,v2)
		fp += abs(1 - is_match(v1,v2))
	return (tp, fp)

def _get_val(k, dmap):
	v = dmap.get(k)
	if v == None:
		v = dmap.get((k[1],k[0]))
	return v

def make_mcpsc_data(ce_map,tm_map,usm_map):
	ret = {}
	for k, v in ce_map.iteritems():
		ret[k] = sum([v,_get_val(k, tm_map)])/2#,_get_val(k, usm_map)])/3
	return ret
	
def get_tpfp_curve(data,astral2pdb40,klassification):
	l = len(data)
	c = map(lambda x : get_tpfp(data[0:int(x*l*1.0/100.0)],astral2pdb40,klassification), range(0,101,10))
	P, N = c[-1]
	return map(lambda x : (x[0]*1.0/P,x[1]*1.0/N), c)

def calc_auc(roc):
	f_auc = lambda x : (roc[x][0]+roc[x-1][0]) * (roc[x][1]-roc[x-1][1]) / 2.0
	return sum(map(f_auc, range(1,len(roc))))

##
APDB40_TO_ASTRAL = dict(map(lambda x : get_astral_to_apdb_line(x), open(APDB40_TO_ASTRAL_FILE_NAME)))
##print APDB40_TO_ASTRAL

ADVTOPS_DATASET = dict(map(lambda x : ((x.split(' ')[0].lower(),x.split(' ')[1].lower()),
  (x.split(' ')[-3].lower(),x.split(' ')[-2].lower())), open(ADVTOPS_DATASET_FILE_NAME)))
##print ADVTOPS_DATASET

ce_map = normalize(dict(get_data_ce_line(x) for x in open(CE_LOG_FILE_NAME) if x[0:5] == 'Algo:'))
##print ce_map

tm_map = normalize(dict(get_data_tm_line(x) for x in open(TM_LOG_FILE_NAME) if x[0:5] == 'Algo:'))
##print tm_map

usm_map = normalize(dict(get_data_usm_line(x) for x in open(USM_LOG_FILE_NAME) if x[0:5] == 'Algo:'))
##print usm_map

mcpsc_map = make_mcpsc_data(ce_map,tm_map,usm_map)
##print mcpsc_map

##
methods = ['ce','tm','usm','mcpsc']
for psc, method in zip([ce_map,tm_map,usm_map,mcpsc_map],methods):
	for klass in ['1','2','3','4']:
		psc_sorted = sorted(filter(lambda x : filter1(x[0][0],x[0][1],APDB40_TO_ASTRAL,ADVTOPS_DATASET,klass), 
			psc.iteritems()),key=lambda x : x[1])
		# print len(psc_sorted)
		psc_tpfp = get_tpfp_curve(psc_sorted,APDB40_TO_ASTRAL,ADVTOPS_DATASET)
		print [method, klass, calc_auc(psc_tpfp)]

