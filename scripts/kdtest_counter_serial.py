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
    
kd = kcollections.Kdict(int, k)
kd.set_merge_func(lambda o, n: o + n)

for i, seq in enumerate(seqs):
    print('adding', i, len(seq))
    kd.add_seq(seq, [1 for _ in range(len(seq))])

print('\t', len(kd), None)

for kmer, val in kd.iteritems():
    if val > 7:
        print(kmer, val)

#print kd['TTTTTTTTGTTTAAGGAGGAAAGAACA']
print('passed')
print('done loading kmer', len(kd), 'kmers')
