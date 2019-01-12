from ._Kset import Kset as KsetParent
from ._Kset import testbit32 as testbit32

class Kset( KsetParent ):
    def __init__( self, k ):
        super( Kset, self ).__init__( k )

    def items( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def iteritems( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def __iter__( self ):
        return self._get_kmers( self.get_root(), self.k )

    def _get_kmers( self, v, k, prefix = '' ):
        # get from UC
        for i in range( self.get_uc_size( v ) ):
            yield prefix + self.get_uc_kmer( v, k, i )

        for i in range( self.get_vs_size( v ) ):
            child_prefix = prefix + self.get_child_suffix( v, i )
            for kmer in self._get_kmers( self.get_child_vertex( v, i ), k - 4, child_prefix ):
                yield kmer

    def __str__( self ):
        return '{' + ','.join( self ) + '}'

    def __repr( self ):
        return self.__str__()

    def copy( self ):
        new_set = Kset( self.k )
        for kmer in self:
            new_set.add( kmer )
        return new_set

    def update( self, *iters ):
        for _iter in iters:
            for item in _iter:
                self.add( item )
        return self

    def discard( self, item ):
        if item in self:
            del self[ item ]

    def pop( self ):
        kmer = next( self._get_kmers( self.get_root(), self.k ) )
        del self[ kmer ]
        return kmer

    def isdisjoint( self, other ):
        for kmer in other:
            if kmer in self:
                return False
        return True

    def issubset( self, other ):
        for kmer in self:
            if kmer not in other:
                return False
        return True

    def issuperset( self, other ):
        for kmer in other:
            if kmer not in self:
                return False
        return True

    def intersection( self, *other_sets ):
        new_set = Kset( self.k )
        for kmer in self:
            good = True
            for other_set in other_sets:
                if kmer not in other_set:
                    good = False
                    break
            if good:
                new_set.add( kmer )
        return new_set

    def intersection_update( self, *other_sets ):
        self = self.intersection( *other_sets )
        return self

    def difference( self, other ):
        new_set = Kset( self.k )
        for kmer in self:
            if kmer not in other:
                new_set.add( kmer )
        return new_set

    def difference_update( self, *other_sets ):
        for other_set in other_sets:
            for kmer in other_set:
                self.discard( kmer )
        return self

    def symmetric_difference( self, other_set ):
        new_set = Kset( self.k )
        for kmer in self.__iter__():
            if kmer not in new_set:
                if (
                        ( kmer not in self and kmer in other_set )
                        or ( kmer in self and kmer not in other_set )
                        ):
                    new_set.add( kmer )
        for kmer in other_set:
            if kmer not in new_set:
                if(
                        ( kmer not in self and kmer in other_set )
                        or (kmer in self and kmer not in other_set )
                        ):
                    new_set.add( kmer )
        return new_set

    def symmetric_difference_update( self, other_set ):
        self = self.symmetric_difference( other_set )
        return self

    def union( self, *other_sets ):
        new_set = Kset( self.k )
        for other_set in other_sets:
            for kmer in other_set:
                new_set.add( kmer )
        return new_set

    def __and__( self, other ):
        return self.intersection( other )

    def __or__( self, other ):
        return self.update( other )

    def __xor__( self, other ):
        return self.symmetric_difference( other )

    def __sub__( self, other ):
        return self.difference( other )

