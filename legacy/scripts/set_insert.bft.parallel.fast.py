#!/usr/bin/env python

import kcollections, sys, time
from tqdm import tqdm

k = int(sys.argv[1])
threads = int(sys.argv[2])
ks = kcollections.Kset(k)
ks.parallel_add_init(threads)

seqs = []
seq = ''
c = 0

start_time = time.time()

with(open(sys.argv[3], 'r')) as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                #seqs.append(seq)
                tstart_time = time.time()
                ks.parallel_add_seq(seq, len(seq))
                telapsed_time = time.time() - tstart_time
                print c, '\tadded seq of len', len(seq), telapsed_time
            c += 1
            seq = ''
        else:
            seq += line.strip()
    if len(seq) > 0:
        #seqs.append(seq)
        tstart_time = time.time()
        ks.parallel_add_seq(seq, len(seq))
        telapsed_time = time.time() - tstart_time
        print c, '\tadded seq of len', len(seq), telapsed_time

ks.parallel_add_join()

elapsed_time = time.time() - start_time
print 'elapsed time:', elapsed_time
#print 'read', len(seqs), 'seqs, adding to ks...'

#for seq in seqs:
#    print '\tadding seq...'
#    sys.stdout.flush()
#    ks.parallel_add_seq(seq, len(seq))

print len(ks), 'kmers'
print 'done!'
print 'checking correctness'

c = 0
if len(sys.argv) > 4:
    with open(sys.argv[3], 'r') as fh:
        seq = ''
        for line in fh:
            if line[0] == '>':
                if len(seq) > 0:
                    for i in range(len(seq) - k + 1):
                        kmer = seq[i : i + k]
                        assert kmer in ks, "not find: " + kmer
                        c += 1
                seq = ''
            else:
                seq += line.strip()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        assert kmer in ks, "not find: " + kmer
        c += 1
print 'checked', c, 'kmers'

del ks
