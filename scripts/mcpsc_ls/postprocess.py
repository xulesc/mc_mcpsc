#!/usr/bin/python

from timeit import default_timer as timer
import numpy as np
import sys
from collections import defaultdict
import operator
import math, random
import itertools

###############################################################################
CE_INFILE='batch_run/results_save/ce_results_1.txt'
FS_INFILE='batch_run/results_save/fast_results_1.txt'
GR_INFILE='batch_run/results_save/results.txt.sim'
TM_INFILE='batch_run/results_save/tm_results_1.txt'
UM_INFILE='batch_run/results_save/usm_results.txt'
GT_INFILE='workspace/data/ground_truth'

#PP_OUTFILE='batch_run/postprocessed.dat'
#PP_OUTFILE='batch_run/postprocessed.allpsc.dat'
PP_OUTFILE='batch_run/postprocessed.dat.fill_avg'
UF_OUTFILE='batch_run/unified.dat'
FV_OUTFILE='batch_run/fv.dat'

PREP_PP=False
DO_NNI=True
DO_AUC=True
DO_SUB=False
DO_UNF=False
DO_FVS=False
DO_MBD=False
DO_CON=False
DO_FILL=False
DO_FILL2=False

###############################################################################
GT_SF_EXTRACT 	= lambda x : x.split('.')[2]
GT_SF_EXTRACT2 	= lambda x : '.'.join(x.split('.')[:3])
EXCLUDE_HASH  	= lambda x : not x.startswith('#')
conv 		= {'a':'1','b':'2','c':'3','d':'4'}
MOST_COMMON	= lambda lst : max(set(lst), key=lst.count)


###############################################################################
def read_psc_data(outfile, gtdata, pm, fname, idx, inv=0, sep=' '):
  start = timer()
  retlist = []
  for line in open(fname):
    data = line.replace('\n','').split(sep)
    try:
      if data[idx] == '-0': continue ## this is for gr-align failed cases
      if line.find('(0.0%)') != -1: continue ## this is for ce failed cases
      if gtdata.get((data[0],data[1])) == None: continue ## drop extra pairs
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
    upsc_vals = psc_vals # [psc_vals[0], psc_vals[1], psc_vals[3]]
    valid_psc_vals = filter(lambda x : x >= 0, upsc_vals)
    if valid_psc_vals != None and len(valid_psc_vals) > 0:
      mcpsc_val = sum(valid_psc_vals) / len(valid_psc_vals)
    ##
    ret.append([data[0], data[1], GT_SF_EXTRACT2(data[4]), 
      GT_SF_EXTRACT2(data[5]), int(data[6])] + psc_vals + [mcpsc_val])
  print len(ret)
  if filter != None:
    ret=filter(flambda,ret)
  print len(ret)
  return ret
  
def homologous_domains(a, b): return GT_SF_EXTRACT2(a) == GT_SF_EXTRACT2(b)

def get_micro_stats(pscdata, thresh, total_homologs, total_non_homologs):
  lpscdata = filter(lambda x : x[2] >= 0 and x[2] <= thresh, pscdata)
  tp = len(filter(lambda x : x[0] == x[1], lpscdata))
  fp = len(lpscdata) - tp
  fn = total_homologs - tp
  tn = total_non_homologs - fp
  return [tp, fp, fn, tn, float(tp)/max(.0001,float(tp+fp)), 
    float(tp)/max(.0001,float(tp+fn)), float(tp)/total_homologs, 
    float(fp)/total_non_homologs]
    
def calc_auc(roc):
  f_auc = lambda x : (roc[x][0]+roc[x-1][0]) * (roc[x][1]-roc[x-1][1]) / 2.0
  return sum(map(f_auc, range(1,len(roc))))                

def read_psc2(fname,colcount,sep=' '):
  print fname
  ret = {}; retlist = []
  for line in open(fname):
    if line.find('-0') != -1: continue
    if line.find('-nan') != -1: continue
    if line[0]=='Q': continue
    data = line.replace('\n','').replace('%','').replace(',','').split(sep)
    if len(data) < colcount: continue
    fdata = []
    for d in data:
      if d.find('(') == -1:
        fdata.append(d)
      else:
        fdata.append(d.split('(')[1].split(')')[0])
    ret[(fdata[0],fdata[1])] = fdata[2:]
  print len(ret)
  return ret  

###############################################################################
if PREP_PP:
  pp_outfile 		= open(PP_OUTFILE,'w')
  gtdata		= read_gt(GT_INFILE,pp_outfile)
  cedata 		= read_psc_data(pp_outfile,gtdata,'ce',CE_INFILE,7)
  fastdata	 	= read_psc_data(pp_outfile,gtdata,'fs',FS_INFILE,5)
  grdata 		= read_psc_data(pp_outfile,gtdata,'gr',GR_INFILE,6,sep='\t',inv=1)
  tmdata 		= read_psc_data(pp_outfile,gtdata,'tm',TM_INFILE,8,inv=1)
  usmdata 		= read_psc_data(pp_outfile,gtdata,'um',UM_INFILE,2)
  
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
    fpscdata = map(lambda x: [x[2],x[3],x[midx]], 
      filter(lambda x : x[midx] >= 0, pscdata))
    total_homologs = len(filter(lambda x : x[0] == x[1], fpscdata))
    total_non_homologs = len(fpscdata) - total_homologs
    of = open('batch_run/perf/%s_perf.dat' %m, 'w')
    for pthresh in xrange(-1,101):
      thresh = pthresh * 1.0 / 100
      stats = get_micro_stats(fpscdata, thresh, total_homologs, 
        total_non_homologs)
      of.write('%s\n' %' '.join(map(lambda x : str(x), [thresh] + stats)))
    of.close()
  print 'done'

if DO_NNI:
  do_nni_core()

###############################################################################
def do_auc_core():
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
if DO_UNF:
  reader1 = lambda x : (x.split(' ')[:2], x.split(' ')[2:])
  reader2 = lambda x : (x.split('\t')[:2], x.split('\t')[2:])
  gtdata = read_gt(GT_INFILE, sys.stdout)
  cedata = read_psc2(CE_INFILE,11)
  fsdata = read_psc2(FS_INFILE,6)
  grdata = read_psc2(GR_INFILE,7,sep='\t')
  tmdata = read_psc2(TM_INFILE,9)
  umdata = read_psc2(UM_INFILE,3)
  
  ofile = open(UF_OUTFILE,'w')
  for k, v in gtdata.iteritems():
    cev = cedata.get(k)
    if cev == None:
      cev = ['-1']*(11-2)
    cev = [cev[1]] + cev[3:]
    fsv = fsdata.get(k)
    if fsv == None:
      fsv = ['-1']*(6-2)
    grv = grdata.get(k)
    if grv == None:
      grv = ['-1']*(7-2)
    tmv = tmdata.get(k)
    if tmv == None:
      tmv = ['-1']*(9-2)
    tmv = tmv[2:]
    umv = umdata.get(k)
    if umv == None:
      umv = ['-1']*(3-2) 
    ofile.write('%s %s %s %s %d %s %s %s %s %s\n' %(k[0],k[1],
      v[2], v[3], (GT_SF_EXTRACT2(v[2]) == GT_SF_EXTRACT2(v[3])),
      ' '.join(cev), ' '.join(fsv), ' '.join(grv),
      ' '.join(tmv), ' '.join(umv)))
  ofile.close()
    
###############################################################################
if DO_FVS:
  ofile = open('%s.csv' %FV_OUTFILE,'w')
  klasses = []; fvs = []
  for line in open(UF_OUTFILE):
    if line.find('-1') != -1: continue # skips lines with incomplete dimensions
    data = line.replace('\n','').split(' ')
    fv = map(lambda x: float(x), data[5:])
    klasses.append(data[4])
    fvs.append(fv)
    
  for dim in range(len(fvs[0])):
    print 'normalizing: %d' %dim
    dim_data = map(lambda x : x[dim], fvs)
    maxv = max(dim_data)
    minv = min(dim_data)
    diff = maxv - minv
    for fv in fvs:
      fv[dim] = (fv[dim]-minv)/(diff)
        
  for k, v in zip(klasses, fvs):    
    fv = ' '.join(map(lambda x : '%d:%f' %(x,v[x]), range(len(v))))
    ofile.write('%s %s\n' %(k,fv))
          
  ofile.close()
  
###############################################################################
if DO_MBD:
  do_consensus = True
  cath_level = 0
  do_top = False
  ## read domains and cath/scop classification
  gtdata = defaultdict(str)
  for line in open(GT_INFILE):
    data = line.replace('\n','').split('\t')
    d1 = data[0]; d2 = data[1]; 
    cath1 = data[-2]; cath2 = data[-1]; iscop = False ## CATH
    #cath1 = data[2]; cath2 = data[3]; iscop = True;  ## SCOP
    gtdata[d1] = cath1
    gtdata[d2] = cath2
  print '%d domains' %len(gtdata)
  ##
  pscdata = read_post_data(PP_OUTFILE)
  methods = ['CE', 'FAST', 'GRALIGN', 'TMALIGN', 'USM', 'MCPSC']
  method_index = [5, 6, 7, 8, 9, 10]
  ## write fv
  fvo = open('fv.dat','w')
  ##
  fpscdatas = []
  for midx in method_index:
    fpscdata = map(lambda x: [x[0],x[1],x[midx]], filter(lambda x : x[midx] >= 0, pscdata))
    datamap = defaultdict(lambda: defaultdict(int))
    for d in fpscdata:
      datamap[d[0]][d[1]] = d[2]
      datamap[d[1]][d[0]] = d[2]
    fpscdatas.append(datamap)
  print 'Method "ignore" "Class 2" "Class 3" "Class 1" "ignore" "Class 3" "Class 1" "Class 2" "ignore"'
  stash = defaultdict(lambda: defaultdict(str))
  for m,midx,ctr in zip(methods,method_index,range(midx)):
    #if m != 'MCPSC': continue
    fout = open('%s.detail'%m,'w')
    datamap = fpscdatas[ctr]
    blockdata = []
    pos = 0; tot = 0;
    for x in range(4):
      blockdata.append([])
      for y in range(4):
        blockdata[x].append(0)
    for k, v in datamap.iteritems():
      if gtdata[k] == None: continue
      cath_query = gtdata[k].split('.')[cath_level]
      if m != 'MCPSC' or do_consensus == False:
        sorted_x = sorted(v.items(), key=operator.itemgetter(1))
        if gtdata[sorted_x[0][0]] == None: continue
        cath_nn = gtdata[sorted_x[0][0]].split('.')[cath_level]
        stash[m][k] = sorted_x[0][0]
      if m == 'MCPSC' and do_consensus:
        cath_nns = map(lambda x : gtdata[x].split('.')[cath_level], map(lambda x : stash[x][k], methods[:-1]))
        cath_nn  = MOST_COMMON(filter(lambda x : x != '', cath_nns))
      if do_top:
        pos += cath_query == cath_nn
        tot += 1
      cath_query = int(cath_query); cath_nn = int(cath_nn)
      fout.write('%d_%d:%s_%s\n' %(cath_query-1,cath_nn-1,k,sorted_x[0]))
      blockdata[cath_query-1][cath_nn-1] += 1
    fout.close()
    if do_top:
      print '# %s accuracy = %f' %(m,pos*1.0/tot)
      continue
    # print blockdata
    cor = sum(map(lambda x : blockdata[x][x], range(len(blockdata)-1)))
    total = sum(map(lambda x : sum(x), blockdata))
    print '# %s correct allocations = %d from total = %d giving %0.2f' %(m,
      cor,total,cor*1.0/total)
    # continue
    for b in blockdata:
      s = sum(b[:3])
      for x in range(3):
        b[x] = b[x]*1.0/s
    op = []
    for x in range(3): 
      for y in range(3): 
        op.append(blockdata[x][y])
    print '%s %s' %(m,' '.join(map(lambda x : str(x), op)))

###############################################################################
if DO_CON: ## should be run only for the smaller set
  cath_level = 0
  ## read domains and cath/scop classification
  gtdata = defaultdict(str)
  for line in open(GT_INFILE):
    data = line.replace('\n','').split('\t')
    d1 = data[0]; d2 = data[1]; 
    cath1 = data[-2]; cath2 = data[-1]; iscop = False ## CATH
    gtdata[d1] = cath1
    gtdata[d2] = cath2
  print '%d domains' %len(gtdata)
  ##
  pscdata = read_post_data(PP_OUTFILE)
  methods = ['CE', 'FAST', 'GRALIGN', 'TMALIGN', 'USM']
  method_index = [5, 6, 7, 8, 9]
  ## write fv
  fvo = open('fv.dat','w')
  ## read all data
  fpscdatas = []
  for midx in method_index:
    fpscdata = map(lambda x: [x[0],x[1],x[midx]], filter(lambda x : x[midx] >= 0, pscdata))
    datamap = defaultdict(lambda: defaultdict(int))
    for d in fpscdata:
      datamap[d[0]][d[1]] = d[2]
      datamap[d[1]][d[0]] = d[2]
    fpscdatas.append(datamap)    
  ## sort
  for fpscdata in fpscdatas:
    for k,v in fpscdata.iteritems():
      fpscdata[k] = sorted(v.items(), key=operator.itemgetter(1))
  dlist = fpscdatas[0].keys()  
  ## 
  # print 'Method "ignore" "Class 2" "Class 3" "Class 1" "ignore" "Class 3" "Class 1" "Class 2" "ignore"'
  combs = []
  lst = range(len(methods))
  for i in xrange(0, len(methods)+1):
    els = [list(x) for x in itertools.combinations(lst, i)]
    combs.extend(els)
  combs = combs[1:]
  for perm in combs:
    # print perm
    perm_methods = map(lambda x : methods[x], perm)
    perm_fpscdatas = map(lambda x : fpscdatas[x], perm)
    pos = 0; top = 0
    for query in dlist:
      # if gtdata.get(query) == None: continue
      cath_query = gtdata[query].split('.')[cath_level]      
      # print perm_fpscdatas[0][query][0][0]
      cath_nns = map(lambda x : gtdata[x[query][0][0]].split('.')[cath_level], perm_fpscdatas)
      cath_nn = MOST_COMMON(filter(lambda x : x != '', cath_nns))
      pos += cath_query == cath_nn
      top += 1
    print '%s : %.2f' %('-'.join(perm_methods), pos*1.0/top)      
  ## clever
  #for query in dlist:
  #  cath_query = gtdata[query].split('.')[cath_level]
  #  cath_nns = map(lambda x : gtdata[x[query][0][0]].split('.')[cath_level], perm_fpscdatas)
  #  if cath_query in ['1','2']:
  #    cath_nn = cath_nns[2]
  #  else:
  #    cath_nn = MOST_COMMON(filter(lambda x : x != '', cath_nns))
  #  pos += cath_query == cath_nn
  #  top += 1
  #print 'Clever: %f' %(pos*1.0/top)
          
###############################################################################
def coverage(gt,data):
  r = lambda x : data.get(x[0]) != None and data.get(x[0]).get(x[1]) != None
  return len(filter(r,gt)) * 1.0 / len(gt)
  
def retrieve(data,q,t):
  qnn = None; tnn = None
  try:
    qnn = data.get(q)[0][0]
    tnn = data.get(t)[0][0]
  except:
    pass
  return [qnn,tnn]            
    
if DO_FILL:
  gtdata = defaultdict(str)
  gtpairs = []
  for line in open(GT_INFILE):
    data = line.replace('\n','').split('\t')
    d1 = data[0]; d2 = data[1]; 
    cath1 = data[-2]; cath2 = data[-1]; iscop = False ## CATH
    gtdata[d1] = cath1
    gtdata[d2] = cath2
    gtpairs.append([d1,d2])
  print '%d domains' %len(gtdata)
  ##
  pscdata = read_post_data(PP_OUTFILE)
  methods = ['CE', 'FAST', 'GRALIGN', 'TMALIGN', 'USM']
  method_index = [5, 6, 7, 8, 9]
  ## read all data
  fpscdatas = []; ofpscdatas = []; averages = []
  for midx in method_index:
    fpscdata = map(lambda x: [x[0],x[1],x[midx]], filter(lambda x : x[midx] >= 0, pscdata))
    datamap = defaultdict(lambda: defaultdict(int)); ofpscdata = {}
    sum = 0
    for d in fpscdata:
      ofpscdata[(d[0],d[1])] = d[2]
      sum += d[2]
    averages.append(sum*1.0/len(fpscdata))
    fpscdatas.append(datamap)
    ofpscdatas.append(ofpscdata)
    print '%s coverage %f' %(methods[midx-5], coverage(gtpairs,datamap))    
  ##
  norms = []
  for ofpscdata in ofpscdatas:
    v = np.array(ofpscdata.values())
    norms.append([np.mean(v),np.std(v)])
  ##
  for pair in gtpairs:
    p1,p2 = pair
    for average,ofpscdata,norm in zip(averages,ofpscdatas,norms):
      if ofpscdata.get((p1,p2)) != None: continue
      ofpscdata[(p1,p2)] = average # random.normalvariate(norm[0],norm[1]) #random.uniform(0,1) #average
      
  fout = open('%s.fill' %PP_OUTFILE, 'w')
  for pair in gtpairs:
    p1,p2 = pair; c1,c2 = [gtdata[p1],gtdata[p2]]
    dummy = '-1'
    psc_scores = map(lambda x : '%f' %x[(p1,p2)], ofpscdatas)
    fout.write('%s %s sc1 sc2 %s %s %s %s\n' %(p1,p2,c1,c2,dummy,' '.join(psc_scores)))
  fout.close()
  
###############################################################################
def getnn(datamap,dom):
  if datamap.get(dom) == None: return None
  return datamap.get(dom)[0][0]
  
if DO_FILL2:
  gtdata = defaultdict(str)
  gtpairs = []
  for line in open(GT_INFILE):
    data = line.replace('\n','').split('\t')
    d1 = data[0]; d2 = data[1]; 
    cath1 = data[-2]; cath2 = data[-1]; iscop = False ## CATH
    gtdata[d1] = cath1
    gtdata[d2] = cath2
    gtpairs.append([d1,d2])
  print '%d domains' %len(gtdata)
  ##
  pscdata = read_post_data(PP_OUTFILE)
  methods = ['CE', 'FAST', 'GRALIGN', 'TMALIGN', 'USM']
  method_index = [5, 6, 7, 8, 9]
  ## read all data
  fpscdatas = []; ofpscdatas = []; averages = []
  for midx in method_index:
    fpscdata = map(lambda x: [x[0],x[1],x[midx]], filter(lambda x : x[midx] >= 0, pscdata))
    datamap = defaultdict(lambda: defaultdict(int)); ofpscdata = {}; sum = 0
    for d in fpscdata:
      datamap[d[0]][d[1]] = d[2]
      datamap[d[1]][d[0]] = d[2]
      ofpscdata[(d[0],d[1])] = d[2]
      sum += d[2]
    averages.append(sum*1.0/len(fpscdata))    
    ofpscdatas.append(ofpscdata)
    fpscdatas.append(datamap)    
    print '%s coverage %f' %(methods[midx-5], coverage(gtpairs,datamap))
  ## sort
  for fpscdata in fpscdatas:
    for k,v in fpscdata.iteritems():
      fpscdata[k] = sorted(v.items(), key=operator.itemgetter(1))
  dlist = fpscdatas[0].keys()  
  ## 
  for pair in gtpairs:
    p1,p2 = pair
    for average,ofpscdata in zip(averages,ofpscdatas):
      if ofpscdata.get((p1,p2)) != None: continue      
      qnns = filter(lambda x : x != None, map(lambda x : getnn(fpscdatas[x],p1), range(5)))
      tnns = filter(lambda x : x != None, map(lambda x : getnn(fpscdatas[x],p2), range(5)))
      if len(qnns) == 0 or len(tnns) == 0:
        ofpscdata[(p1,p2)] = average
        continue
      # print 'here'
      qnn = MOST_COMMON(qnns); tnn = MOST_COMMON(tnns)
      if ofpscdata.get((p1,tnn)) != None:
        s = ofpscdata.get((p1,tnn))
      elif ofpscdata.get((qnn,p2)) != None:
        s = ofpscdata.get((qnn,p2))
      elif ofpscdata.get((qnn,tnn)) != None:
        s = ofpscdata.get((qnn,tnn))
      else:
        s = average
      ofpscdata[(p1,p2)] = s
  ##
  fout = open('%s.fill' %PP_OUTFILE, 'w')
  for pair in gtpairs:
    p1,p2 = pair; c1,c2 = [gtdata[p1],gtdata[p2]]
    dummy = '-1'
    psc_scores = map(lambda x : '%f' %x[(p1,p2)], ofpscdatas)
    fout.write('%s %s sc1 sc2 %s %s %s %s\n' %(p1,p2,c1,c2,dummy,' '.join(psc_scores)))
  fout.close()
  
      
###############################################################################

