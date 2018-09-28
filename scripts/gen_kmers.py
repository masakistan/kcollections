import itertools, random, sys

bases=['A','T','G','C']
k = 27
random.seed( 1337 )

for p in itertools.product( bases, repeat = k ):
    prob = random.random()
    if prob < 0.0000005:
        print ''.join( p )
        sys.stdout.flush()
