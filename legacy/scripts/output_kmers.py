import sys

k = 27
seq = ''
with open( sys.argv[ 1 ], 'r' ) as fh:
    for line in fh:
        if line[ 0 ] == '>':
            continue
        seq += line.strip()

for i in range( len( seq ) - k + 1 ):
    print seq[ i : i + k ]
