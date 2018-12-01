#!/usr/bin/env python

import sys, time
from tqdm import tqdm

k = int(sys.argv[1])
seqs = []
seq = ''
c = 0
ks = set()
start_time = time.time()
with(open(sys.argv[2], 'r')) as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                #seqs.append(seq)
                tstart_time = time.time()
                for i in xrange(len(seq) - k + 1):
                    kmer = seq[i : i + k]
                    ks.add(kmer)
                telapsed_time = time.time() - tstart_time
                print c, '\tadded seq of len', len(seq), telapsed_time
            c += 1
            seq = ''
        else:
            seq += line.strip()
    if len(seq) > 0:
        #seqs.append(seq)
        tstart_time = time.time()
        for i in xrange(len(seq) - k + 1):
            kmer = seq[i : i + k]
            ks.add(kmer)
        telapsed_time = time.time() - tstart_time
        print c, '\tadded seq of len', len(seq), telapsed_time

#print 'read', len(seqs), 'seqs, adding to ks...'

elapsed_time = time.time() - start_time
print 'elapsed time:', elapsed_time
#for seq in seqs:
#    print '\tadding seq...'
#    sys.stdout.flush()
    #for i in tqdm(xrange(len(seq) - k + 1)):

print len(ks), 'kmers'
print 'done!'
print 'checking correctness'

c = 0
if len(sys.argv) > 3:
    with open(sys.argv[2], 'r') as fh:
        seq = ''
        for line in fh:
            if line[0] == '>':
                if len(seq) > 0:
                    for i in range(len(seq) - k + 1):
                        kmer = seq[i : i + k]
                        assert kmer in ks
                        c += 1
                seq = ''
            else:
                seq += line.strip()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        assert kmer in ks
        c += 1
print 'checked', c, 'kmers'
del ks
