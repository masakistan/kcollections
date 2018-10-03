import kc


kmers = [
        "AAATCGTGGAAACGGTTCGGGTCCAGC",
        "CCATTATGATATCTGCCAGTTGCCACA",
        "TGGGCGGTGTGATGACGACGCTGATCG",
        "GCTATTCCGGCGTACCTTTGCCCTCCA",
        "CGCCGTAGAAAATGCCCATGGCAAGAA",
        "TTGGGGAATATATGCAGTATTGGGGAA",
        "GCTAAGCTATCGTCGTATCCGTAACAC",
        "GCTACCGTGAACGGTGCTACCTCCTTA",
        "GGTATACGGGAAGGCAGGCATTGGCTG",
        "CTGCTCGGTTTCCTCATCATCAAAATC",
        "CGTTACCGTGCAAAGCAGCCTCGATGC",
        "TGAATCTGTGTGGTGCACGCCGCACGG",
        "TCGGCTTTGCCACGTCCCGCCAGTTCA",
        "CCGGCGTTCCTTGATAACCCACGCATG",
        "ATGCGAATGGCAGCATTCATATTGGTC",
        "ATTTGCGCCATGGCAATGAAAAGCCAC",
        "ACCCGTTAAAATGAAATATAAGAGACG",
        "GAATGTATCAGCCGATGGTTCTACGAT",
        "ACGCAAACTTTTTGCGAAGGTGGCGTG",
        "GCTTTGATGAAAGCTTTTGGTGCGATG",
        "CGCTGACGTTGCCCCATGTGAGCGTGA",
        "AGTGCCGGACACATTGGATGTATGGTT",
        "TCCGTGGTTGGCGCAGCGGAGGCGCTT",
        "GTAACGGTGCGGGCTGACGCGTACAGG",
        "CCAACCGTCTGGCGGAGCTGGCCCAGC",
        "GCGCCGTTGTTCGACCACTTTATCGAG",
        "TACGGTCGCCATATACAAGTAGTGCTG",
        "AACCCGAAAAACGGTCGTCTGATTGTT",
        "ATCCGCAAACACCAGATCGCTTTAGGG",
        "GGTTCCCGCTGGCGCAATTGAAAACTT",
        "TAATCGACGCCGGAAGGTTTGATCACA",
        "ATATTTAACGACAGCGCGTGCAAATTG",
        "TCATACTTTTTCCATTTCAATTAACCG",
        "TCGCCGACCGGTTCGGTCAATGCCGCC",
        "GCGTGGTGCCCAGCGGTTTCAACACCA",
        "AACGCCTCAGAATACTTTACTGGGGCT",
        "TGATCGAATAACTAATACGGTTCTCTG"
        ]

other_kmers = [
        "GAGGCCAACGCCCATAATGCGGGCTGT",
        "GATTGCCGGTGATGCCGCCGGAATGTG",
        "GGCTCCCGCTTGCACGATCAACCGCCC",
        "CGGAAATCGTAAAAACGGATTTCATAA",
        "TGCATAATGCAGCCATCCTGAATATTG",
        "CGTTGGCGGTGCGCTGCTGGAGCAACT",
        "ACGCCGGTGAGCTGGCCGGAGCCGGGA",
        "TAGAAACCCACTTGAGCGGCGGACGAT",
        "GGGGGCTGCGACTGGTGACCGATGCCG",
        "AAGAGCTGCAAGAAAGCTTCGGTGGCC",
        "TCTCAGGCGGTACGTAACGAAGCAAAA",
        "TTCTGGAGATGCAATGAAGATTATTAC",
        "CTTTAAAGCACTCTTCAATTTGGGTAA",
        "GATGTGTGGAAACTGACGGTCAAAAAC",
        "GAACTCCGCTGAAAATTATGCCATAGG",
        "CATAACCAGCGCGATGCATTCGCGGAA",
        "ATTAAGATAAATCTTACCATTTCTACG",
        "GAAATTCTCAATAAATGCGGTAACTTA",
        "TGCGCGCGCGGTTCCAGCGTTTGGGTA",
        "TGGGTATTGAAGATCAGGCGGGCAGGA"
        ]

def get_kmers( v, prefix = "", depth = 0 ):
    for uc_bkmer in v.get_uc().get_bkmers():
        yield prefix + uc_bkmer.get_seq()
    for cc in v.get_ccs():
        for clust_idx, sfc in enumerate( cc.get_suf_clust_data() ):
            for suffix in get_kmers( sfc.get_child_vertex(), prefix + cc.prefix_from_clust( clust_idx ), depth + 1 ):
                yield suffix

k = 27
kd = kc.Kdict( k, kc.calc_bk( k ) )


print 'Insert kmers...'
for kmer in kmers:
    kd.insert( kmer )
assert len( kmers ) == kd.size(), 'ERROR: Kmer insertion failed'
print '\tPass!'


print 'Insert more kmers...'
for kmer in other_kmers:
    kd.insert( kmer )
assert len( kmers ) + len( other_kmers ) == kd.size(), 'ERROR: Kmer insertion failed'
print '\tPass!'


print 'Remove kmers...'
for kmer in other_kmers:
    kd.remove( kmer )
assert len( kmers ) == kd.size(), 'ERROR: Kmer removal failed!'
print '\tPass!'


print 'Check for kmers'
for kmer in kmers:
    #print '\tkmer', kmer, 'contained', kd.contains( kmer )
    assert kd.contains( kmer ), 'ERROR: kd does not contain kmer ' + kmer


for i, kmer in enumerate( get_kmers( kd.get_root() ) ):
    #print i, kmer
    assert kmer in kmers, 'ERROR: Kmer exists in kd but not in list of inserted kmers'

print 'Passed all tests!'


