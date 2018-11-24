import sys

k = int(sys.argv[1])
seqs = []
seq = ''
with open(sys.argv[2], 'r') as fh:
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                seqs.append(seq)
            seq = ''
        else:
            seq += line.strip()

for seq in seqs:
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        print kmer


