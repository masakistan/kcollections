import sys
from tqdm import tqdm
import kc

k = 27
kd = kc.Kdict( k )
bkmer = kc.Bkmer( k )

c = 0
kmers = []
with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        kd.insert( kmer )
        #kmers.append( kmer )
        #bkmer.set_seq( kmer, k )
        #print kmer
        #kd.insert_bkmer( bkmer )
        
        
        c += 1
        #if c % 1000000 == 0:
        #    print '\nstart insert proc'
        #    sys.stdout.flush()
        #    kmers.sort()
        #    for kmer in tqdm( kmers ):
        #        kd.insert( kmer )
        #    kmers = []
        #    print '\tdone'
        #    sys.stdout.flush()
        
        #if c % 10000 == 0:
        #    sys.stdout.write( '.' )
        #    sys.stdout.flush()
    sys.stdout.write( '\n' )
    sys.stdout.flush()

assert kd.size() == c
print 'Done! Processed', kd.size(), 'kmers!'

with open( sys.argv[ 1 ], 'r' ) as fh:
    c = 0
    for line in tqdm( fh ):
        kmer = line.strip().split()[ 0 ]
        
        #bkmer.set_seq( kmer, k )
        #assert kd.contains_bkmer( bkmer ), 'ERROR: line ' + str( c ) + ' kmer: ' + kmer
        
        assert kd.contains( kmer ), 'ERROR: line ' + str( c ) + ' kmer: ' + kmer
        c += 1

print 'Done! Passed contains test!'

