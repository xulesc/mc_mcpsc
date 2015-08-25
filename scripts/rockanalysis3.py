#!/usr/bin/env python

INFILE = 'apdb40.psc.advTops.log'

def get_tpfp(psc):
	tp = 0; fp = 0
	for (d1,d2,k1,k2),s in psc:
		tp += k1 == k2
		fp += k1 != k2
	return (tp, fp)

def get_tpfp_curve(data):
	l = len(data)
	c = map(lambda x : get_tpfp(data[0:int(x*l*1.0/100.0)]), range(0,101,10))
	P, N = c[-1]
	return map(lambda x : (x[0]*1.0/P,x[1]*1.0/N), c)

def calc_auc(roc):
	f_auc = lambda x : (roc[x][0]+roc[x-1][0]) * (roc[x][1]-roc[x-1][1]) / 2.0
	return sum(map(f_auc, range(1,len(roc))))

for x in ['1','2','3','4']:
	ces = []; tms = []; usms = []
	for line in open(INFILE):
		data = line.replace('\n','').split(' ')
		if data[4][0] != x or data[5][0] != x:
			continue
		sk1 = data[4].split('.'); sk2 = data[5].split('.')
		k1 = int(sk1[1]); k2 = int(sk2[1])
		ces.append(((data[0],data[1],k1,k2),float(data[-3])))
		tms.append(((data[0],data[1],k1,k2),float(data[-2])))
		usms.append(((data[0],data[1],k1,k2),float(data[-1])))
	print len(ces)
	## normalize
	ce_v = map(lambda x : x[1], ces)
	tm_v = map(lambda x : x[1], tms)
	usm_v = map(lambda x : x[1], usms)
	min_ce_v = min(ce_v); max_ce_v = max(ce_v); dce = max_ce_v - min_ce_v
	min_tm_v = min(tm_v); max_tm_v = max(tm_v); dtm = max_tm_v - min_tm_v
	min_usm_v = min(usm_v); max_usm_v = max(usm_v); dusm = max_usm_v - min_usm_v
	ce_v = dict(map(lambda x : (x[0],(x[1]-min_ce_v)/dce), ces))
	tm_v = dict(map(lambda x : (x[0],(x[1]-min_tm_v)/dtm), tms))
	usm_v = dict(map(lambda x : (x[0],(x[1]-min_usm_v)/dusm), usms))
	##
	mcpsc_v = {}
	for k,v in ce_v.iteritems():
		mcpsc_v[k] = sum([ce_v[k],tm_v[k]])/2
	##
	print ['ce',x,calc_auc(get_tpfp_curve(sorted(ce_v.iteritems(),key=lambda x: x[1])))]
	print ['tm',x,calc_auc(get_tpfp_curve(sorted(tm_v.iteritems(),key=lambda x: x[1])))]
	print ['usm',x,calc_auc(get_tpfp_curve(sorted(usm_v.iteritems(),key=lambda x: x[1])))]
	print ['mcpsc',x,calc_auc(get_tpfp_curve(sorted(mcpsc_v.iteritems(),key=lambda x: x[1])))]

