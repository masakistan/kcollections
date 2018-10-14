from _Kdict as KdictS
from _Kset as 


class Kdict( KdictParent ):
    def __init__( self, k ):
        super().__init__( k )

    def __iter__( self ):
        return self._get_kmers( super.kc.v, super.k )

    def _get_kmers( v, k, prefix = '' ):
        uc = v.uc
        luc = uc.size
        for i in range( luc ):
            yield prefix + get_kmer_from_uc( uc, k, i )


class Kset( KsetParent ):
    def __init__( self, k ):
        super().__init__( k )

    def __iter__( self ):
        return self._get_kmers( super.kc.v, super.k )

    def _get_kmers( v, k, prefix = '' ):
        uc = v.uc
        luc = uc.size
        for i in range( luc ):
            yield prefix + get_kmer_from_uc( uc, k, i )

        lcc = v.cc_size
        for i in range( lcc ):
            pass


