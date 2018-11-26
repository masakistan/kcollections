#!/usr/bin/env python

import kcollections, sys, time
from tqdm import tqdm

k = int(sys.argv[1])
seqs = []
seq = ''
with(open(sys.argv[2], 'r')) as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                seqs.append(seq)
            seq = ''
        else:
            seq += line.strip()
    if len(seq) > 0:
        seqs.append(seq)

print 'read', len(seqs), 'seqs, adding to ks...'

ks = kcollections.Kset(k)
ks.parallel_add_init(4)

for seq in seqs:
    print '\tadding seq...'
    sys.stdout.flush()
    ks.parallel_add_seq(seq, len(seq))
ks.parallel_add_join()

print len(ks), 'kmers'
print 'done!'
