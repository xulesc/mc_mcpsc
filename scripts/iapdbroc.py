#!/usr/bin/env python

APDB40_TO_ASTRAL_FILE_NAME = 'apdb40_astral_mapping.clean'
ADVTOPS_DATASET_FILE_NAME = 'advTOPSplus-PDB40_subset-results.txt'
CE_LOG_FILE_NAME = 'apdb40.psc.log.ce'
TM_LOG_FILE_NAME = 'apdb40.psc.log.tm'
USM_LOG_FILE_NAME = 'apdb40.psc.log.usm'

def get_astral_to_apdb_line(line):
        dom, fname = line.split(' ')
        return (dom, fname[fname.rfind('/')+1:fname.rfind('.')])

APDB40_TO_ASTRAL = dict(map(lambda x : get_astral_to_apdb_line(x), open(APDB40_TO_ASTRAL_FILE_NAME)))
##print APDB40_TO_ASTRAL

ADVTOPS_DATASET = dict(map(lambda x : ((x.split(' ')[0].lower(),x.split(' ')[1].lower()),
  (x.split(' ')[-3].lower(),x.split(' ')[-2].lower())), open(ADVTOPS_DATASET_FILE_NAME)))

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

def is_match(v1,v2):
	return v1.split('_')[1][2] == v2.split('_')[1][2]

def normalize(data):
	vals = map(lambda x : float(x[1]), data.iteritems())
	minv = min(vals); maxv = max(vals); diffv = maxv-minv
	return dict(map(lambda x : (x[0], (float(x[1]) - minv)/diffv), data.iteritems()))

def _get_val(k, dmap):
	v = dmap.get(k)
	if v == None:
		v = dmap.get((k[1],k[0]))
	return v

def filter1(d1,d2,astral2pdb40,klassification,klass):
	ad1 = astral2pdb40[d1]; ad2 = astral2pdb40[d2]
	v = klassification.get((ad1,ad2))
	if v == None:
		v = klassification.get((ad2,ad1))
	if v == None:
		raise ValueError("mismatch %s %s" %(d1,d2))
	return v[0].split('_')[1][0] == klass and v[1].split('_')[1][0] == klass

def make_mcpsc_data(ce_map,tm_map,usm_map):
	ret = {}
	for k, v in ce_map.iteritems():
		ret[k] = sum([v,_get_val(k, tm_map)])/2#,_get_val(k, usm_map)])/3
	return ret
	        
def advastral():
	ret = {}
	for k, v  in ADVTOPS_DATASET.iteritems():
		if APDB40_TO_ASTRAL.get(k[0]) != None and APDB40_TO_ASTRAL.get(k[1]) != None:
			ret[(APDB40_TO_ASTRAL[k[0]], APDB40_TO_ASTRAL[k[1]])] = (v[0], v[1])
	return ret

ce_map = normalize(dict(get_data_ce_line(x) for x in open(CE_LOG_FILE_NAME) if x[0:5] == 'Algo:'))
##print ce_map

tm_map = normalize(dict(get_data_tm_line(x) for x in open(TM_LOG_FILE_NAME) if x[0:5] == 'Algo:'))
##print tm_map

usm_map = normalize(dict(get_data_usm_line(x) for x in open(USM_LOG_FILE_NAME) if x[0:5] == 'Algo:'))
##print usm_map

mcpsc_map = make_mcpsc_data(ce_map,tm_map,usm_map)
##print mcpsc_map

pairs_of_interest = advastral()
klass = 1
psc = ce_map
#psc_sorted = sorted(filter(lambda x : filter1(x[0][0],x[0][1],APDB40_TO_ASTRAL,ADVTOPS_DATASET,klass), 
#	psc.iteritems()),key=lambda x : x[1])
#print psc
#for (d1,d2), v in psc.iteritems():
for (d1,d2), (klass1,klass2) in ADVTOPS_DATASET.iteritems():
	if APDB40_TO_ASTRAL.get(d1) == None or APDB40_TO_ASTRAL.get(d2) == None:
		continue
	f1 = APDB40_TO_ASTRAL[d1]
	f2 = APDB40_TO_ASTRAL[d2]		
	k = (f1,f2)
	v = psc.get(k)
	if v == None:
		k = (f2,f1)
	if psc.get(k) == None:
		continue
	v_ce = ce_map[k]
	v_tm = tm_map[k]
	v_usm = usm_map[k]
	print '%s %s %s %s %s %s %f %f %f' %(d1,d2,f1,f2,klass1.split('_')[1],
		klass2.split('_')[1],v_ce,v_tm,v_usm)
	
