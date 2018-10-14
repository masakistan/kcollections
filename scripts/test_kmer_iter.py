import sys
from _Kdict import Kdict as KdictParent
from _Kset import Kset as KsetParent


class Kdict( KdictParent ):
    def __init__( self, k ):
        super( Kdict, self ).__init__( k )

    def __iter__( self ):
        return self._get_kmers( self.get_root(), self.k )

    def _get_kmers( self, v, k, prefix = '' ):
        # get from UC
        for i in range( self.get_uc_size( v ) ):
            yield prefix + self.get_uc_kmer( v, k, i )

        for i in range( self.get_cc_size( v ) ):
            for j in range( self.get_cc_child_size( v, i ) ):
                child_prefix = prefix + self.get_cc_child_suffix( v, i, j )
                for kmer in self._get_kmers( self.get_cc_child_vertex( v, i, j ), k - 4, child_prefix ):
                    yield kmer

    def iteritems( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def keys( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield kmer

    def values( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield self[ kmer ]


class Kset( KsetParent ):
    def __init__( self, k ):
        super( Kset, self ).__init__( k )

    def __iter__( self ):
        return self._get_kmers( self.get_root(), self.k )

    def _get_kmers( self, v, k, prefix = '' ):
        # get from UC
        for i in range( self.get_uc_size( v ) ):
            yield prefix + self.get_uc_kmer( v, k, i )

        for i in range( self.get_cc_size( v ) ):
            for j in range( self.get_cc_child_size( v, i ) ):
                child_prefix = prefix + self.get_cc_child_suffix( v, i, j )
                for kmer in self._get_kmers( self.get_cc_child_vertex( v, i, j ), k - 4, child_prefix ):
                    yield kmer

k = 27
ks = Kset( k )
c = 0
'''with open( sys.argv[ 1 ], 'r' ) as fh:
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
'''

kd = Kdict( k )
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


