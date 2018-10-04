import sys
from tqdm import tqdm
import kc

k = 27
kd = set()

c = 0
with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        kd.add( kmer )
        c ++ 1

assert len( kd ) == c
print 'Done! Processed', kd.size(), 'kmers!'

with open( sys.argv[ 1 ], 'r' ) as fh:
    c = 0
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        assert kmer in kd, 'ERROR: line ' + str( c ) + ' kmer: ' + kmer
        c += 1

print 'Done! Passed contains test!'
