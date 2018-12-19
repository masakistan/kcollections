from ._Kdict import Kdict as KdictParent
from ._Kset import Kset as KsetParent
from ._Kcounter import Kcounter as KcounterParent
from ._Kcolor import Kcolor as KcolorParent


class Kcolor( KcolorParent ):
    def __init__( self, k ):
        super( Kcolor, self ).__init__( k )

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
        res = []
        for key, val in self.items():
            res.append( key + ':' + str( val ) )
        return '{' + ','.join( res ) + '}'

    def __repr__( self ):
        return self.__str__()

    def items( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def iteritems( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def keys( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield kmer

    def values( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield self[ kmer ]

    def copy( self ):
        new_kcolor = Kcolor( self.k )
        for kmer in self:
            new_kcolor[ kmer ] = self[ kmer ]
        return new_kcolor

    def get( self, key, value = None ):
        if key in self:
            return self[ key ]
        else:
            return value

    def popitem( self ):
        kmer, item = next( self.items() )
        del self[ kmer ]
        return ( kmer, item )

    def setdefault( self, key, value = None ):
        if key not in self:
            self[ key ] = value
            return value
        else:
            return self[ key ]

    def pop( self, key, *default ):
        if key in self:
            value = self[ key ]
            del self[ key ]
            return value
        else:
            if len( default ) > 0:
                return default
            else:
                raise KeyError( key )

    def update( self, *others ):
        for other in others:
            for item in other:
                try:
                    key, val = item
                except:
                    key = item
                    val = other[ key ]
                self[ key ] = val


class Kdict( KdictParent ):
    def __init__( self, k ):
        super( Kdict, self ).__init__( k )

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
        res = []
        for key, val in self.items():
            res.append( key + ':' + str( val ) )
        return '{' + ','.join( res ) + '}'

    def __repr__( self ):
        return self.__str__()

    def items( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def iteritems( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def keys( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield kmer

    def values( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield self[ kmer ]

    def copy( self ):
        new_kdict = Kdict( self.k )
        for kmer in self:
            new_kdict[ kmer ] = self[ kmer ]
        return new_kdict

    def get( self, key, value = None ):
        if key in self:
            return self[ key ]
        else:
            return value

    def popitem( self ):
        kmer, item = next( self.items() )
        del self[ kmer ]
        return ( kmer, item )

    def setdefault( self, key, value = None ):
        if key not in self:
            self[ key ] = value
            return value
        else:
            return self[ key ]

    def pop( self, key, *default ):
        if key in self:
            value = self[ key ]
            del self[ key ]
            return value
        else:
            if len( default ) > 0:
                return default
            else:
                raise KeyError( key )

    def update( self, *others ):
        for other in others:
            for item in other:
                try:
                    key, val = item
                except:
                    key = item
                    val = other[ key ]
                self[ key ] = val



class Kset( KsetParent ):
    def __init__( self, k ):
        super( Kset, self ).__init__( k )

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


class Kcounter( KcounterParent ):
    def __init__( self, k ):
        super( Kcounter, self ).__init__( k )

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
        res = []
        for key, val in self.items():
            res.append( key + ':' + str( val ) )
        return '{' + ','.join( res ) + '}'

    def __repr__( self ):
        return self.__str__()

    def items( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def iteritems( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield ( kmer, self[ kmer ] )

    def keys( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield kmer

    def values( self ):
        for kmer in self._get_kmers( self.get_root(), self.k ):
            yield self[ kmer ]

    def copy( self ):
        new_kcounter = Kcounter( self.k )
        for kmer in self:
            new_kcounter[ kmer ] = self[ kmer ]
        return new_kcounter

    def get( self, key, value = 0 ):
        if key in self:
            return self[ key ]
        else:
            return value

    def popitem( self ):
        kmer, item = next( self.items() )
        del self[ kmer ]
        return ( kmer, item )

    def setdefault( self, key, value = 0 ):
        if key not in self:
            self[ key ] = value
            return value
        else:
            return self[ key ]

    def pop( self, key, *default ):
        if key in self:
            value = self[ key ]
            del self[ key ]
            return value
        else:
            if len( default ) > 0:
                return default
            else:
                raise KeyError( key )

    def update( self, *others ):
        for other in others:
            for item in other:
                try:
                    key, val = item
                except:
                    key = item
                    val = other[ key ]
                self[ key ] = val
