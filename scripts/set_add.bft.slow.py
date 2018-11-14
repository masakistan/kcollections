#!/usr/bin/env python

import kcollections, sys

k = int(sys.argv[1])
ks = kcollections.Kset(k)

seq = ''
with(open(sys.argv[2], 'r')) as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                for i in range(len(seq) - k + 1):
                    kmer = seq[i : i + k]
                    ks.add(kmer)
            seq = ''
        else:
            seq += line.strip()
    if len(seq) > 0:
        for i in range(len(seq) - k + 1):
                kmer = seq[i : i + k]
                ks.add(kmer)
print 'done!'
