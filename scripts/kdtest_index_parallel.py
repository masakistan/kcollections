import kcollections, sys
from Bio import SeqIO
k = 27

seqs = []
for record in SeqIO.parse(sys.argv[1], 'fasta'):
    #print len(record.seq)
    seqs.append(str(record.seq))
    break

for i, seq in enumerate(seqs):
    print('seqs', i, len(seq))

def merge_func(o, n):
    o.append(n[0])
    return o
    
kd = kcollections.Kdict((list, int), k)
kd.set_merge_func(merge_func)
kd.parallel_add_init(4)

for seq in seqs:
    print('adding', len(seq))
    kd.parallel_add_seq(seq, iter([[i] for i in range(len(seq))]))
kd.parallel_add_join()

print('kmer set contains:', len(kd), 'kmers')

for kmer, val in kd.iteritems():
    if isinstance(val, list) and len(val) > 7:
        print(kmer, val)
                
