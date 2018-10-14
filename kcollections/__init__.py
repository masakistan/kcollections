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

