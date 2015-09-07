#!/usr/bin/env python

import sys
import os
import subprocess
import re
import zlib
from multiprocessing import Pool
from timeit import default_timer as timer
import threading
import shutil

# CONFIGURATION
PROGRAMS = ['ce', 'tmalign', 'fast', 'gralign', 'usm']
PROGDIR = 'programs'
DATADIR = 'data'
WORKDIR = 'work'
PDBEXTN = '.ent'
THREADS = 16
# END OF CONFIGURATION

FASTOUT = open('work/fast_results.txt', 'w')
MCOUT = open('work/mc_results.txt', 'w')
TMOUT = open('work/tm_results.txt', 'w')
CEOUT = open('work/ce_results.txt', 'w')

lock = threading.Lock()

# PRE-PROCESS


class GRALIGN_PRE_PROCESSOR:

    """
    Provides pre-processing functionality needed by GRALIGN.

    The class provides functionality to run the two step pre-processing
    required by the GRALIGN program. The first step is provided by the
    GRALIGN program CMAP and the second step is provided by the GRALIGN
    program DCOUNT. This class wraps the process of running these two
    programs.

    At the end of the process basic number statistics of how many PDB files
    were processed, how many contact maps were generated and how many
    signature files were generated is output.
    """

    def __init__(self, binPath1, binPath2, tmpDir):
        self._binPath1 = binPath1
        self._binPath2 = binPath2
        self._tmpDir = '%s/gralign' % tmpDir
        if not os.path.exists(self._tmpDir):
            os.makedirs(self._tmpDir)

    def pre_process_all_to_all(self, pdbDir):
        for file in os.listdir(pdbDir):
            # ./CMap -i 1amk.pdb -c A -o 1amkA.gw -d 12.0
            # ./DCount -i 1amkA.gw -o 1amkA.ndump
            f1 = file.replace(PDBEXTN, '.gw')
            f2 = file.replace(PDBEXTN, '.ndump')
            cmd = './%s -i %s/%s -c A -o %s/%s -d 12.0' % (
                self._binPath1, pdbDir, file, self._tmpDir, f1)
            proc = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
            # print proc.stdout.readlines()
            cmd = './%s -i %s/%s -o %s/%s' % (
                self._binPath2, self._tmpDir, f1, self._tmpDir, f2)
            proc = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
            # print proc.stdout.readlines()
        l1 = len(os.listdir(pdbDir))
        l2 = len(
            filter(lambda x: x.endswith(".gw"), os.listdir(self._tmpDir)))
        l3 = len(
            filter(lambda x: x.endswith(".ndump"), os.listdir(self._tmpDir)))
        print 'pdb files: %d, contact maps: %d, signature files: %d'
        %(l1, l2, l3)


# END OF PRE-PROCESS

# HANDLERS
class CE_HANDLER:

    def __init__(self, binPath, tmpDir):
        self._binPath = binPath
        self._tmpDir = tmpDir
        if not os.path.exists(self._tmpDir):
            os.makedirs(self._tmpDir)

    def process_pair(self, fname1, fname2):
        global lock
        # note special requirement for pom/mkDB from exec location
        # ./CE - $PDB_DIR/$1 - $PDB_DIR/$2 - $TMP_DIR
        f1 = fname1.split('/')[-1].replace(PDBEXTN, '')
        f2 = fname2.split('/')[-1].replace(PDBEXTN, '')
        wdir = '%s/%s_%s' % (self._tmpDir, f1, f2)
        os.makedirs(wdir)
        cmd = './%s - %s - %s - %s' % (self._binPath, fname1, fname2, wdir)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        new_pair = True
        res = []
        for line in proc.stdout.readlines():
            line = line.replace('\n', '')
            if line.find('Size=') != -1:
                if new_pair:
                    new_pair = False
                    res.append([f1, f2])
                data = line.split(' ')
                dom_chain = data[2]
                chain = dom_chain.split(':')[1]
                dlen = data[3].replace('(', '').replace(')', '').replace(
                    '=', '').replace('Size', '')
                res[len(res) - 1] += [chain, dlen]

            if line.find('Alignment length') != -1:
                new_pair = True
                data = line.split(' ')
                res[len(res) - 1] += [
                    data[3], data[6].replace('A', ''), data[9], data[12],
                    data[19]]
        if len(res) > 1:
            print 'pair: %s-%s has multi chain domain' % (fname1, fname2)
        # each entry in res: fname1, fname2, chain1, len1, chain2, len2,
        # Alignment length, Rmsd, Z-Score, Gaps, Sequence identities
        if len(res) == 0:
            res.append([f1, f2])
        with lock:
            CEOUT.write('%s\n' % ' '.join(res[0]))
        shutil.rmtree(wdir, ignore_errors=True)
        return res


class TM_HANDLER:

    def __init__(self, binPath):
        self._binPath = binPath

    def process_pair(self, fname1, fname2):
        global lock
        # ./$prg $dir/$x $dir/$y
        cmd = './%s %s %s' % (self._binPath, fname1, fname2)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        f1 = fname1.split('/')[-1].replace(PDBEXTN, '')
        f2 = fname2.split('/')[-1].replace(PDBEXTN, '')
        ret = [[f1, f2]]
        for line in proc.stdout.readlines():
            line = line.replace('\n', '')
            if line.find('Length of Chain_') != -1:
                ret[0].append(
                    line.split(':')[1].replace('residues', '').strip())
            if line.find('Aligned length') != -1:
                data = re.split('\s+', line)
                ret[0] += [data[2], data[4], data[6]]
            if line.find('TM-score=') != -1:
                ret[0].append(line.split(' ')[1])
        # each entry in res: fname1, fname2, len1, len2, aligned length, RMSD,
        # sequence id, tm1, tm2
        if len(ret) == 0:
            ret.append([f1, f2])
        with lock:
            TMOUT.write('%s\n' % ' '.join(ret[0]))
        return ret


class FAST_HANDLER:

    def __init__(self, binPath):
        self._binPath = binPath

    def process_pair(self, fname1, fname2):
        global lock
        cmd = './%s %s %s' % (self._binPath, fname1, fname2)  # , pdb1, pdb2)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        ret = [
            [fname1.split('/')[-1].replace(PDBEXTN, ''), fname2.split('/')[-1].
                replace(PDBEXTN, '')]]
        for line in proc.stdout.readlines():
            if line.find('RMSD=') != -1:
                data = line.replace('\n', '').split(' ')
                ret[0] += [data[0].split('=')[1], data[3].split('=')[
                    1], data[4].split('=')[1], data[5].split('=')[1]]
        with lock:
            FASTOUT.write('%s\n' % ' '.join(ret[0]))
        # each entry of res: fname1, fname2, len1, len2, rmsd
        return ret

    def process_alltoall(self, pdb_pairs):
        fast_out_dir = '%s/fast' % WORKDIR
        os.makedirs(fast_out_dir)
        fast_out = open('%s/results.txt' % fast_out_dir, 'w')
        for pdb_pair in pdb_pairs:
            score = self.process_pair(pdb_pair[0], pdb_pair[1])
            fast_out.write('%s\n' % (' '.join(score[0])))
        fast_out.close()


class MC_HANDLER:

    def __init__(self, binPath):
        self._binPath = binPath

    def process_pair(self, fname1, fname2):
        global lock
        cmd = './%s %s %s' % (self._binPath, fname1, fname2)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        l1 = '-1'
        l2 = '-1'
        pairs = '-1'
        rmsd = '-1'
        maxsub = '-1'
        dlen = '-1'
        grmsd = '-1'
        tm = '-1'
        sid = '-1'
        f1 = fname1.split('/')[-1].replace(PDBEXTN, '')
        f2 = fname2.split('/')[-1].replace(PDBEXTN, '')

        for line in proc.stdout.readlines():
            try:
                line = line.replace('\n', '')
                if line.find('Prediction size=') != -1:
                    l1 = line.split(' ')[4]
                if line.find('Experiment size=') != -1:
                    l2 = line.split(' ')[4]
                if line.find('Pairs=') != -1:
                    data = line.split(' ')
                    pairs = data[3].replace(
                        ',', '')
                    rmsd = data[5].replace(',', '')
                    maxsub = data[6].split('=')[1]
                    dlen = data[8]
                    grmsd = data[10]
                    tm = data[
                        11].split('=')[1].replace(',', '')
                    sid = data[12].split('=')[1]
            except:
                pass
        ret = [[f1, f2, l1, l2, pairs, rmsd, maxsub, dlen, grmsd, tm, sid]]
        with lock:
            MCOUT.write('%s\n' % ' '.join(ret[0]))
        return ret


class GR_HANDLER:

    def __init__(self, binPath):
        self._binPath = binPath

    def process_alltoall(self, dirName):
        #
        q = '%s/pairs.lst' % dirName
        f = open(q, 'w')
        files = filter(lambda x: x.endswith(".ndump"), os.listdir(dirName))
        for file in files:
            f.write('%s\n' % file.replace('.ndump', ''))
        f.close()
        # ./GR-Align -q skolnick.lst -r ./skolnick -o results.txt
        cmd = '%s -q %s -r %s -o %s' % (
            self._binPath, q, dirName, '%s/results.txt' % dirName)
        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in proc.stdout.readlines():
            print line


class USM_HANDLER:

    def __init__(self):
        self._MAX_CONTACTS = 1000

    def _trunc(self, x):
        if len(x) < self._MAX_CONTACTS:
            return x
        return x[:self._MAX_CONTACTS]

    def process_pair(self, fname1, fname2):
        reader = lambda x: not x.startswith('|')
        cm1 = ''.join(self._trunc(filter(reader, open(fname1).readlines())))
        cm2 = ''.join(self._trunc(filter(reader, open(fname2).readlines())))
        comp = lambda x: float(len(zlib.compress(x)))
        x = comp(cm1)
        y = comp(cm2)
        xy = comp(
            cm1 + cm2)
        yx = comp(cm2 + cm1)
        f1 = fname1.split('/')[-1].replace('.gw', '')
        f2 = fname2.split('/')[-1].replace('.gw', '')
        return [[f1, f2, '%f' % (max(yx - y, xy - x) / max(x, y))]]

    def process_alltoall(self, cm_pairs):
        usm_out_dir = '%s/usm' % WORKDIR
        os.makedirs(usm_out_dir)
        usm_out = open('%s/results.txt' % usm_out_dir, 'w')
        for cm_pair in cm_pairs:
            score = self.process_pair(cm_pair[0], cm_pair[1])
            dom1 = cm_pair[0].split('/')[-1].replace('.gw', '')
            dom2 = cm_pair[1].split('/')[-1].replace('.gw', '')
            usm_out.write('%s %s %f\n' % (dom1, dom2, score))
        usm_out.close()

# END OF HANDLERS

# PROCESS
start = timer()
gralign_pre_processor = GRALIGN_PRE_PROCESSOR(
    'programs/CMap', 'programs/DCount', 'work')
gralign_pre_processor.pre_process_all_to_all('data')
end = timer()
print 'gralign preprocessing took %d seconds' % (end - start)

ce = CE_HANDLER('%s/%s' % (PROGDIR, PROGRAMS[0]), '%s/ce' % WORKDIR)
# print ce.process_pair('%s/1aa9.pdb' %DATADIR, '%s/1ash.pdb' %DATADIR)

tm = TM_HANDLER('%s/%s' % (PROGDIR, PROGRAMS[1]))
# print tm.process_pair('%s/1aa9.pdb' %DATADIR, '%s/1ash.pdb' %DATADIR)

fast = FAST_HANDLER('%s/%s' % (PROGDIR, PROGRAMS[2]))
# print fast.process_pair('%s/1aa9.pdb' %DATADIR, '%s/1ash.pdb' %DATADIR)

mc = MC_HANDLER('%s/%s' % (PROGDIR, PROGRAMS[3]))
# print mc.process_pair('%s/1aa9.pdb' %DATADIR, '%s/1ash.pdb' %DATADIR)

usm = USM_HANDLER()
# print usm.process_pair('work/gralign/1aa9.gw', 'work/gralign/1ash.gw')

gr = GR_HANDLER('%s/%s' % (PROGDIR, PROGRAMS[4]))

pdb_files = os.listdir(DATADIR)
psc_pairs = []
for line in open('../workspace/data/ground_truth'):
    ldata = line.split('\t')
    psc_pairs.append(
        ('%s/%s.ent' % (DATADIR, ldata[0]), '%s/%s.ent' % (DATADIR, ldata[1])))
cm_files = filter(lambda x: x.endswith(".gw"), os.listdir('work/gralign'))
cm_pairs = []
for x in range(0, len(cm_files)):
    for y in range(x, len(cm_files)):
        cm_pairs.append(('work/gralign/%s' % cm_files[x], 'work/gralign/%s'
                        % cm_files[y]))

# Do Serial
start = timer()
gr.process_alltoall('work/gralign')
end = timer()
print 'gralign processed in %d seconds' % (end - start)
# usm.process_alltoall(cm_pairs)


def fast_process_pair(ps): return fast.process_pair(ps[0], ps[1])


def usm_process_pair(ps): return usm.process_pair(ps[0], ps[1])


def mc_process_pair(ps): return mc.process_pair(ps[0], ps[1])


def tm_process_pair(ps): return tm.process_pair(ps[0], ps[1])


def ce_process_pair(ps): return ce.process_pair(ps[0], ps[1])


def doMulti(outfilename, procmethod, pairs):
    start = timer()
    p = Pool(THREADS)
    reduced = p.map(procmethod, pairs)
    out = open(outfilename, 'w')
    for res in reduced:
        out.write('%s\n' % ' '.join(res[0]))
    out.close()
    return timer() - start

usm_seconds = doMulti('%s/usm_results.txt' % WORKDIR, usm_process_pair,
                      cm_pairs)
print 'usm processed in %d seconds' % usm_seconds
fast_seconds = doMulti('%s/fast_results_1.txt' % WORKDIR, fast_process_pair,
                       psc_pairs)
print 'fast processed in %d seconds' % fast_seconds
FASTOUT.close()

tm_seconds = doMulti('%s/tm_results_1.txt' % WORKDIR, tm_process_pair,
                     psc_pairs)
print 'tmalign processed in %d seconds' % tm_seconds
TMOUT.close()
ce_seconds = doMulti('%s/ce_results_1.txt' % WORKDIR, ce_process_pair,
                     psc_pairs)
print 'ce processed in %d seconds' % ce_seconds
CEOUT.close()


# END OF PROCESS
