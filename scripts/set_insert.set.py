#!/usr/bin/env python

import sys, time

k = int(sys.argv[1])
ks = set()

with open(sys.argv[2], 'r') as fh:
    for line in fh:
        line = line.strip()
        ks.add(line)
print 'done!'
