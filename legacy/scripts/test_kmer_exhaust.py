from __future__ import print_function
import sys, random, itertools
from tqdm import tqdm
import kcollections

bases=['A','T','G','C']

good = []
bad = []

k = int( sys.argv[ 1 ] )
ks = kcollections.Kset( k )

for kmer in tqdm( itertools.product(bases, repeat=k) ):
    kmer = ''.join( kmer )
    if( random.random() > .75 ):
        good.append( kmer )
    else:
        bad.append( kmer )

random.shuffle( good )
random.shuffle( bad )

c = 0
for kmer in good:
    ks.add( kmer )
    c += 1

assert len( ks ) == c
print('Done! Processed', len( ks ), 'kmers!')

for kmer in tqdm( good ):
    assert kmer in ks, 'ERROR: kmer not found! ' + kmer

for kmer in tqdm( bad ):
    assert kmer not in ks, 'ERROR: erroneous kmer found! ' + kmer

print('Done! Passed contains test!')
