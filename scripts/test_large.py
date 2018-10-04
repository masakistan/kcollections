import sys
from tqdm import tqdm
import kc

k = 27
kd = kc.Kdict( k )
bkmer = kc.Bkmer( k )

c = 0
with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        kd.insert( kmer )
        c += 1

assert kd.size() == c
print 'Done! Processed', kd.size(), 'kmers!'

with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        assert kd.contains( kmer ), 'ERROR: line ' + str( c ) + ' kmer: ' + kmer

print 'Done! Passed contains test!'

