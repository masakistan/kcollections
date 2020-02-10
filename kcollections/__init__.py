from ._Kdict import *
from ._Kset import Kset as KsetParent
from ._Kcounter import Kcounter as KcounterParent


def create_kdict(base):
    class tkdict(base):
        def __init__(self, k=0, caster=None, seq_caster=None, rcaster=None):
            super(tkdict, self).__init__(k)
            self.caster = caster
            self.rcaster = rcaster
            self.seq_caster = seq_caster

        def __str__( self ):
            res = []
            for key, val in self.items():
                res.append( key + ':' + str( val ) )
            return '{' + ','.join( res ) + '}'

        '''
        def __getitem__(self, key):
            if key in self:
                val = super(tkdict, self).__getitem__(key)
                print(self.rcaster)
                try:
                    return self.rcaster(val)
                except:
                    return val
            else:
                 raise KeyError("kmer {} not in kdict".format(key))
        '''
        def __setitem__(self, key, val):
            #print('setitem')
            if self.caster is not None:
                #print('casting', val)
                val = self.caster(val)
                #print(type(val))
            super(tkdict, self).__setitem__(key, val)

        def __repr__( self ):
            return self.__str__()

        def items( self ):
           for kmer, val in self.__iter__():
               try:
                   yield kmer, self.rcaster(val)
               except:
                   yield kmer, val

        def parallel_add_seq(self, seq, values):
            if self.seq_caster:
                values = self.seq_caster(map(self.caster, values))
            super(tkdict, self).parallel_add_seq(seq, values)

        def add_seq(self, seq, values):
            if self.seq_caster:
                values = self.seq_caster(map(self.caster, values))
            super(tkdict, self).add_seq(seq, values)
            
        def iteritems( self ):
           for kmer, val in self.__iter__():
               try:
                   yield kmer, self.rcaster(val)
               except:
                   yield kmer, val

        def keys( self ):
            for kmer, val in self.__iter__():
                yield kmer

        def values( self ):
            for kmer, val in self.__iter__():
                try:
                    yield self.rcaster(val)
                except:
                    yield val

        def copy( self ):
            new_kdict = Kdict( self.k )
            for kmer in self:
                new_kdict[ kmer ] = self[ kmer ]
            return new_kdict

        def get( self, key, value = None ):
            if key in self:
                val = self[key]
                try:
                    return self.rcaster(val)
                except:
                    return val
            else:
                return value

        def popitem( self ):
            kmer, item = next( self.items() )
            del self[ kmer ]
            try:
                item = self.rcaster(item)
            except:
                pass
            return ( kmer, item )

        def setdefault( self, key, value = None ):
            if key not in self:
                self[ key ] = value
                return value
            else:
                value = self[ key ]
                try:
                    return self.rcaster(value)
                except:
                    return value

        def pop( self, key, *default ):
            if key in self:
                value = self[ key ]
                del self[ key ]
                try:
                    return self.rcaster(value)
                except:
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
                    try:
                        self[ key ] = self.caster(val)
                    except:
                        self[key] = val

                        '''
        def add_seq(self, seq, values):
            if self.caster:
                caster = self.caster.__name__
                split = caster.find('_')
                caster = eval(caster[:split] + '_vector' + caster[split:])
                values = caster(values)
            super(tkdict, self).add_seq(seq, values)

        def parallel_add_seq(self, seq, values):
            if self.caster:
                print('casting')
                caster = self.caster.__name__
                split = caster.find('_')
                caster = eval(caster[:split] + '_vector' + caster[split:])
                print(caster)
                values = caster(values)
                print(type(values))
            super(tkdict, self).parallel_add_seq(seq, values)
'''
            

    return tkdict

def Kdict(val_type, k):
    try:
        iter(val_type)
    except:
        type_name = 'Kdict_' + val_type.__name__
        caster = None
        rcaster = None
        seq_caster = None
        kd = create_kdict(eval(type_name))(k)
    else:
        type_name = 'Kdict_' + '_'.join([x.__name__ if x != list else 'vector' for x in val_type])
        # NOTE: for now, we just store everything as a list
        caster = 'o' + '_'.join([x.__name__ if x != list else 'vector' for x in val_type])
        seq_caster = 'ovector_' + '_'.join([x.__name__ if x != list else 'vector' for x in val_type])
        # NOTE: this does not allow for nexted collections
        rcaster = val_type[0]
        kd = create_kdict(eval(type_name))(k, eval(caster), eval(seq_caster), rcaster)
    return kd


class Kset(KsetParent):
    def __init__(self, k=0):
        super(Kset, self).__init__(k)

    def __str__( self ):
        return '{' + ','.join( self ) + '}'

    def __repr__( self ):
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
        kmer = next( self.__iter__())
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


class Kcounter(KcounterParent):
    def __init__(self, k=0):
        super(Kcounter, self).__init__(k)

    def __str__( self ):
        res = []
        for key, val in self.items():
            res.append( key + ':' + str( val ) )
        return '{' + ','.join( res ) + '}'

    def __repr__( self ):
        return self.__str__()

    def keys( self ):
        for kmer, val in self.__iter__():
            yield kmer

    def values( self ):
        for kmer, val in self.__iter__():
            yield val

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
