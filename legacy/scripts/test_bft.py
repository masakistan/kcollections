import sys
from tqdm import tqdm
import kcollections

k = int( sys.argv[ 1 ] )
ks = kcollections.Kset( k )

c = 0
with open( sys.argv[ 2 ], 'r' ) as fh:
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        ks.add( kmer )
        c += 1

assert len( ks ) == c
print 'Done! Processed', len( ks ), 'kmers!'

with open( sys.argv[ 2 ], 'r' ) as fh:
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        assert kmer in ks, 'ERROR: line ' + str( c ) + ' kmer: ' + kmer

print 'Done! Passed contains test!'

