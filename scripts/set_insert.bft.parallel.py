#!/usr/bin/env python

import kcollections, sys, time

k = int(sys.argv[1])
ks = kcollections.Kset(k)
ks.parallel_insert_init(4)

with open(sys.argv[2], 'r') as fh:
    for line in fh:
        line = line.strip()
        print line
        ks.parallel_insert(line)
print 'done!'
