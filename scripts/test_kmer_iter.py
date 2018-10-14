import sys, kcollections


k = 27
ks = kcollections.Kset( k )
c = 0
with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in fh:
        kmer = line.strip().split()[ 0 ]
        ks.add( kmer )
        assert kmer in ks, str( c ) + ' ' + kmer + ' not in ks'
        c += 1

print 'Done insering kmers into ks!'

c = 0
for idx, kmer in enumerate( ks ):
    print idx, kmer
    c += 1
assert c == len( ks )
print 'Printed a total of', c, 'kmers set!'

kd = kcollections.Kdict( k )
c = 0
with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in fh:
        kmer = line.strip().split()[ 0 ]
        kd[ kmer ] = c
        assert kd[ kmer ] == c, str( c ) + ' ' + kmer + ' value in dictionary does not match!'
        assert kmer in kd, str( c ) + ' ' + kmer + ' not in kd'
        c += 1

#print 'Done inserting kmers into kd!'

#for k in kd.keys():
#    print k, kd[ kmer ]

#print 'Dont iterating over keys in kd!'

#for v in kd.values():
#    print v

c= 0
for idx, ( k, v ) in enumerate( kd.iteritems() ):
    print idx, k, v
    c += 1
assert c == len( kd )
print 'Printed a total of', c, 'kmers from dictionary!'

print kd[ 'AAATTCCAACCGTGTCTTCTCCATTAG' ]


