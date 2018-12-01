#!/usr/bin/env python

import kcollections, sys, time

k = int(sys.argv[1])
ks = kcollections.Kset(k)
seqs = []
seq = ''
with(open(sys.argv[2], 'r')) as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                #seqs.append(seq)
                print 'adding seq of len', len(seq)
                ks.add_seq(seq, len(seq))
            seq = ''
        else:
            seq += line.strip()
    if len(seq) > 0:
        print 'adding seq of len', len(seq)
        ks.add_seq(seq, len(seq))
        #seqs.append(seq)


#print 'read', len(seqs), 'seqs, adding to ks...'

#for seq in seqs:
#    print '\tadding seq...'
#    sys.stdout.flush()
#    ks.add_seq(seq, len(seq))
print len(ks), 'kmers'
print 'done!'
