#!/usr/bin/env python

import sys
import time

start_k = int(sys.argv[1])
end_k = int(sys.argv[2])

seqs = []
seq = ''
with(open(sys.argv[3], 'r')) as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                seqs.append(seq)
            seq = ''
        else:
            seq += line.strip()
    if len(seq) > 0:
        seqs.append(seq)

for k in range( start_k, end_k + 1 ):
    start_time = time.time()
    ks = set()
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            ks.add(kmer)
    elapsed_time = time.time() - start_time
    print k, elapsed_time
    sys.stdout.flush()
    del ks
print 'done!'
#print len(ks)

#for kmer in ks:
#    print kmer
#print len(ks)
