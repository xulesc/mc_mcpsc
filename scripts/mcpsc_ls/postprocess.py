#!/usr/bin/python

from timeit import default_timer as timer
import numpy as np

###############################################################################
CE_INFILE='batch_run/results_save/ce_results_1.txt'
FS_INFILE='batch_run/results_save/fast_results_1.txt'
GR_INFILE='batch_run/results_save/results.txt.sim'
TM_INFILE='batch_run/results_save/tm_results_1.txt'
UM_INFILE='batch_run/results_save/usm_results.txt'
GT_INFILE='workspace/data/tsv.3.2.0.VS.1.75.outpairs'

PP_OUTFILE='batch_run/postprocessed.dat'

PREP_PP=False

###############################################################################
GT_SF_EXTRACT = lambda x : x.split('.')[2]

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
  normanp = list((anp-anp.min())/(anp.max()-anp.min()))
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

###############################################################################
if PREP_PP:
  pp_outfile 		= open(PP_OUTFILE,'w')
  cedata 		= read_psc_data(pp_outfile,'ce',CE_INFILE,7)
  fastdata	 	= read_psc_data(pp_outfile,'fs',FS_INFILE,5)
  grdata 		= read_psc_data(pp_outfile,'gr',GR_INFILE,6,sep='\t')
  tmdata 		= read_psc_data(pp_outfile,'tm',TM_INFILE,8,inv=1)
  usmdata 		= read_psc_data(pp_outfile,'um',UM_INFILE,2)
  gtdata		= read_gt(pp_outfile,GT_INFILE)

  for k, v in gtdata.iteritems():
    cev = cedata.get(k); fastv = fastdata.get(k)
    grv = grdata.get(k); tmv = tmdata.get(k)
    usmv = usmdata.get(k)
    pp_outfile.write('%s %s %s %s %d %f %f %f %f %f\n' %(k[0], k[1], v[0], v[1],
      v[-1], max(-1,cev), max(-1,fastv), max(-1,grv), max(-1,tmv), max(-1,usmv)))
  pp_outfile.close()

###############################################################################
