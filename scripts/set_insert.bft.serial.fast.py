#!/usr/bin/env python

import kcollections, sys, time

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

ks = kcollections.Kset(k)

print 'read', len(seqs), 'seqs, adding to ks...'

for seq in seqs:
    print '\tadding seq...'
    sys.stdout.flush()
    ks.add_seq(seq, len(seq))
print 'done!'
