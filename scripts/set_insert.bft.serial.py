#!/usr/bin/env python

import kcollections, sys, time

k = int(sys.argv[1])
ks = kcollections.Kset(k)

with open(sys.argv[2], 'r') as fh:
    for line in fh:
        line = line.strip()
        ks.add(line)
print 'done!'
