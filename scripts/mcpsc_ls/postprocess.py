#!/usr/bin/python

from timeit import default_timer as timer
import numpy as np

###############################################################################
CE_INFILE='batch_run/results_save/ce_results_1.txt'
FS_INFILE='batch_run/results_save/fast_results_1.txt'
GR_INFILE='batch_run/results_save/results.txt.sim'
TM_INFILE='batch_run/results_save/tm_results_1.txt'
UM_INFILE='batch_run/results_save/usm_results.txt'
GT_INFILE='workspace/data/ground_truth'

PP_OUTFILE='batch_run/postprocessed.dat'

PREP_PP=False
DO_NNI=True
DO_AUC=True
DO_SUB=False

###############################################################################
GT_SF_EXTRACT 	= lambda x : x.split('.')[2]
GT_SF_EXTRACT2 	= lambda x : '.'.join(x.split('.')[0:2])
EXCLUDE_HASH  	= lambda x : not x.startswith('#')

###############################################################################
def read_psc_data(outfile, pm, fname, idx, inv=0, sep=' '):
  start = timer()
  retlist = []
  for line in open(fname):
    data = line.replace('\n','').split(sep)
    try:
      if data[idx] == '-0': ## this is for gr-align failed cases
        continue
      retlist.append(((data[0],data[1]), abs(inv-float(data[idx]))))
    except:
      pass
  anp = np.array(map(lambda x : x[1], retlist))
  normanp = list((anp-min(anp))/(max(anp)-min(anp)))
  ret = dict(map(lambda x : (retlist[x][0], normanp[x]), xrange(len(normanp))))
  outfile.write('# %s data read count %d in %d seconds\n' 
    %(pm,len(ret),timer() - start))
  return ret

def read_gt(fname, outfile):
  start = timer()
  ret = {}
  for line in open(fname):
    data = line.replace('\n','').split('\t')
    sf1 = GT_SF_EXTRACT(data[2]); sf2 = GT_SF_EXTRACT(data[3])
    ret[(data[0],data[1])] = [data[2],data[3],data[6],data[7],sf1 == sf2]
  outfile.write('# gt data read count %d in %d seconds\n'
    %(len(ret),timer() - start))
  return ret
  
def read_post_data(fname, flambda=None):
  ret = []
  for line in filter(EXCLUDE_HASH, open(fname)):
    data = line.replace('\n','').split(' ')
    psc_vals = map(lambda x : float(x), data[7:])
    ##
    mcpsc_val = -1.0
    # we use only ce,fast,tmalign because the others are below 50% coverage
    upsc_vals = [psc_vals[0], psc_vals[1], psc_vals[3]]
    valid_psc_vals = filter(lambda x : x >= 0, upsc_vals)
    if valid_psc_vals != None and len(valid_psc_vals) > 0:
      mcpsc_val = sum(valid_psc_vals) / len(valid_psc_vals)
      #mcpsc_val = min(valid_psc_vals)
    ##
    ret.append([data[0], data[1], GT_SF_EXTRACT2(data[4]), 
      GT_SF_EXTRACT2(data[5]), int(data[6])] + psc_vals + [mcpsc_val])
    #print ret[len(ret)-1]      
  print len(ret)
  if filter != None:
    ret=filter(flambda,ret)
  print len(ret)
  return ret
  
def homologous_domains(a, b): return GT_SF_EXTRACT2(a) == GT_SF_EXTRACT2(b)

def get_micro_stats(pscdata, thresh, total_homologs, total_non_homologs, midx):
  lpscdata = filter(lambda x : x[midx] >= 0 and x[midx] <= thresh, pscdata)
  tp = len(filter(lambda x : x[2] == x[3], lpscdata))
  fp = len(lpscdata) - tp
  fn = total_homologs - tp
  tn = total_non_homologs - fp
  return [tp, fp, fn, tn, float(tp)/max(.0001,float(tp+fp)), 
    float(tp)/max(.0001,float(tp+fn)), float(tp)/total_homologs, 
    float(fp)/total_non_homologs]
    
def calc_auc(roc):
  f_auc = lambda x : (roc[x][0]+roc[x-1][0]) * (roc[x][1]-roc[x-1][1]) / 2.0
  return sum(map(f_auc, range(1,len(roc))))                

###############################################################################
if PREP_PP:
  pp_outfile 		= open(PP_OUTFILE,'w')
  cedata 		= read_psc_data(pp_outfile,'ce',CE_INFILE,7)
  fastdata	 	= read_psc_data(pp_outfile,'fs',FS_INFILE,5)
  grdata 		= read_psc_data(pp_outfile,'gr',GR_INFILE,6,sep='\t')
  tmdata 		= read_psc_data(pp_outfile,'tm',TM_INFILE,8,inv=1)
  usmdata 		= read_psc_data(pp_outfile,'um',UM_INFILE,2)
  gtdata		= read_gt(GT_INFILE,pp_outfile)
  
  for k, v in gtdata.iteritems():
    cev = cedata.get(k); fastv = fastdata.get(k)
    grv = grdata.get(k); tmv = tmdata.get(k)
    usmv = usmdata.get(k)
    pp_outfile.write('%s %s %s %s %s %s %d %f %f %f %f %f\n' %(k[0], k[1], 
      v[0], v[1], v[2], v[3], v[-1], max(-1,cev), max(-1,fastv), max(-1,grv), 
      max(-1,tmv), max(-1,usmv)))
  pp_outfile.close()

###############################################################################
def do_nni_core(flambda=None):
  pscdata = read_post_data(PP_OUTFILE, flambda)
  methods = ['ce', 'fast', 'gralign', 'tmalign', 'usm', 'mcpsc']
  method_index = [5, 6, 7, 8, 9, 10]
  for m,midx in zip(methods,method_index):
    print 'processing %s' %m
    fpscdata = filter(lambda x : x[midx] >= 0, pscdata)
    total_homologs = len(filter(lambda x : x[2] == x[3], fpscdata))
    total_non_homologs = len(fpscdata) - total_homologs
    of = open('batch_run/perf/%s_perf.dat' %m, 'w')
    for pthresh in xrange(-1,101):
      thresh = pthresh * 1.0 / 100
      stats = get_micro_stats(fpscdata, thresh, total_homologs, 
        total_non_homologs, midx)
      of.write('%s\n' %' '.join(map(lambda x : str(x), [thresh] + stats)))
    of.close()
  print 'done'

#flambda=lambda x : x[2][0] == sclass or x[3][0] == sclass
if DO_NNI:
  do_nni_core()

###############################################################################
def do_auc_core():
  #roc_data_dirs=['batch_run/perf-cath-avg', 'batch_run/perf-cath-min',
  #  'batch_run/perf-scop-avg', 'batch_run/perf-scop-min']
  roc_data_dirs=['batch_run/perf']
  roc_data_meth=['ce_perf.dat', 'fast_perf.dat', 'gralign_perf.dat',
    'mcpsc_perf.dat', 'tmalign_perf.dat', 'usm_perf.dat']
  roc_data_files=[]
  for rd in roc_data_dirs:
    for rf in roc_data_meth:
      roc_data_files.append('%s/%s' %(rd,rf))
  for roc_data_file in roc_data_files:
    roc_data=map(lambda x : x.split(' ')[-2:], open(roc_data_file))
    roc_data=map(lambda x : (float(x[0]), float(x[1])), roc_data )
    print '%s %f' %(roc_data_file, calc_auc(roc_data))

if DO_AUC:
  do_auc_core()        

###############################################################################
if DO_SUB:
  for sclass in ['1','2','3','4']:
    print sclass
    do_nni_core(flambda=lambda x : x[2][0] == sclass or x[3][0] == sclass)
    do_auc_core()

###############################################################################

