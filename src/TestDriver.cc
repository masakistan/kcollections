#include <iostream>
#include "Kdict.h"

int n_test_insert_kmers = 37;
char* test_insert_kmers[] = {
    ( char* ) "AAATCGTGGAAACGGTTCGGGTCCAGC",
    ( char* ) "CCATTATGATATCTGCCAGTTGCCACA",
    ( char* ) "TGGGCGGTGTGATGACGACGCTGATCG",
    ( char* ) "TACCTTCCGGCGTACCTTTGCCCTCCA",
    ( char* ) "CGCCGTAGAAAATGCCCATGGCAAGAA",
    ( char* ) "TTGGGGAATATATGCAGTATTGGGGAA",
    ( char* ) "CGATAGCTATCGTCGTATCCGTAACAC",
    ( char* ) "GCTACCGTGAACGGTGCTACCTCCTTA",
    ( char* ) "GGTATACGGGAAGGCAGGCATTGGCTG",
    ( char* ) "CTGCTCGGTTTCCTCATCATCAAAATC",
    ( char* ) "CGTTACCGTGCAAAGCAGCCTCGATGC",
    ( char* ) "TGAATCTGTGTGGTGCACGCCGCACGG",
    ( char* ) "TCGGCTTTGCCACGTCCCGCCAGTTCA",
    ( char* ) "CCGGCGTTCCTTGATAACCCACGCATG",
    ( char* ) "ATGCGAATGGCAGCATTCATATTGGTC",
    ( char* ) "ATTTGCGCCATGGCAATGAAAAGCCAC",
    ( char* ) "ACCCGTTAAAATGAAATATAAGAGACG",
    ( char* ) "GAATGTATCAGCCGATGGTTCTACGAT",
    ( char* ) "ACGCAAACTTTTTGCGAAGGTGGCGTG",
    ( char* ) "GCTTTGATGAAAGCTTTTGGTGCGATG",
    ( char* ) "CGCTGACGTTGCCCCATGTGAGCGTGA",
    ( char* ) "AGTGCCGGACACATTGGATGTATGGTT",
    ( char* ) "TCCGTGGTTGGCGCAGCGGAGGCGCTT",
    ( char* ) "GTAACGGTGCGGGCTGACGCGTACAGG",
    ( char* ) "CCAACCGTCTGGCGGAGCTGGCCCAGC",
    ( char* ) "GCGCCGTTGTTCGACCACTTTATCGAG",
    ( char* ) "TACGGTCGCCATATACAAGTAGTGCTG",
    ( char* ) "AACCCGAAAAACGGTCGTCTGATTGTT",
    ( char* ) "ATCCGCAAACACCAGATCGCTTTAGGG",
    ( char* ) "GGTTCCCGCTGGCGCAATTGAAAACTT",
    ( char* ) "TAATCGACGCCGGAAGGTTTGATCACA",
    ( char* ) "ATATTTAACGACAGCGCGTGCAAATTG",
    ( char* ) "TCATACTTTTTCCATTTCAATTAACCG",
    ( char* ) "TCGCCGACCGGTTCGGTCAATGCCGCC",
    ( char* ) "GCGTGGTGCCCAGCGGTTTCAACACCA",
    ( char* ) "AACGCCTCAGAATACTTTACTGGGGCT",
    ( char* ) "TGATCGAATAACTAATACGGTTCTCTG"
};

int n_test_bad_kmers = 20;
char* test_bad_kmers[] = {
    ( char* ) "GAGGCCAACGCCCATAATGCGGGCTGT",
    ( char* ) "GATTGCCGGTGATGCCGCCGGAATGTG",
    ( char* ) "GGCTCCCGCTTGCACGATCAACCGCCC",
    ( char* ) "CGGAAATCGTAAAAACGGATTTCATAA",
    ( char* ) "TGCATAATGCAGCCATCCTGAATATTG",
    ( char* ) "CGTTGGCGGTGCGCTGCTGGAGCAACT",
    ( char* ) "ACGCCGGTGAGCTGGCCGGAGCCGGGA",
    ( char* ) "TAGAAACCCACTTGAGCGGCGGACGAT",
    ( char* ) "GGGGGCTGCGACTGGTGACCGATGCCG",
    ( char* ) "AAGAGCTGCAAGAAAGCTTCGGTGGCC",
    ( char* ) "TCTCAGGCGGTACGTAACGAAGCAAAA",
    ( char* ) "TTCTGGAGATGCAATGAAGATTATTAC",
    ( char* ) "CTTTAAAGCACTCTTCAATTTGGGTAA",
    ( char* ) "GATGTGTGGAAACTGACGGTCAAAAAC",
    ( char* ) "GAACTCCGCTGAAAATTATGCCATAGG",
    ( char* ) "CATAACCAGCGCGATGCATTCGCGGAA",
    ( char* ) "ATTAAGATAAATCTTACCATTTCTACG",
    ( char* ) "GAAATTCTCAATAAATGCGGTAACTTA",
    ( char* ) "TGCGCGCGCGGTTCCAGCGTTTGGGTA",
    ( char* ) "TGGGTATTGAAGATCAGGCGGGCAGGA"
};


int test_cc_contain( CContainer* cc, Bkmer* bkmer )
{
    char* seq = bkmer->get_seq();
    std::cout << "\t\tCompressed container stored " << seq;
    free( seq );
    if( cc->may_contain( bkmer ) )
    {
        if( cc->contains_prefix( bkmer ) )
        {
            std::cout << " True!" << std::endl;
        }
        else
        {
            std::cout << " False positive!" << std::endl;
        }
    }
    else
    {
        std::cout << " False!" << std::endl;
    }
    return 0;
}

int main()
{
    // Testing Parameters
    int k = 27;
    int bk = calc_bk( k );

    std::cout << "Testing kmer dictionary implementation!" << std::endl;
    std::cout << "***************************************" << std::endl;


    /****************************************************************
     * Test serialization and deserialization capabilities
     * *************************************************************/
    char* kmer = ( char* ) "AAATTCCAACCGTGTCTTCTCCATTAG";
    std::cout << "\tTesting kmer serialization/deserialization..." << std::endl << std::flush;
    std::cout << "\t\tkmer: " << kmer << std::endl;

    Bkmer* bkmer = new Bkmer( k, bk, kmer );
    for( int i = 0; i < bk; i++ )
    {
        std::cout << "\t\t\t" << unsigned( bkmer->get_bseq()[ i ] ) << std::endl;
    }

    char* seq = bkmer->get_seq();
    std::cout << "\t\tdeserialized binary kmer: " << seq << std::endl;
    free( seq );
    delete bkmer;
    /****************************************************************
     * End serialization tests
     * *************************************************************/


    /****************************************************************
     * Test Uncompressed Container
     * *************************************************************/
    std::cout << "\tTesting uncompressed container..." << std::endl << std::flush;
    UContainer* uc = new UContainer();
    bkmer = new Bkmer( k, bk, kmer );

    for( int i = 0; i < 2; i++ )
    {
        uc->insert( bkmer );
    }
    std::cout << "\t\tUC size: " << uc->size() << std::endl;

    if( uc->contains( bkmer ) )
    {
        std::cout << "\t\tUC kmer contains success!" << std::endl;
    }
    else
    {
        std::cout << "\t\tUC kmer contains failed!" << std::endl;
    }
    delete bkmer;
    delete uc;
    /****************************************************************
     * End uncompressed container tests
     * *************************************************************/


    /****************************************************************
     * Test Compressed Container
     * *************************************************************/
    std::cout << "\tTesting compressed container..." << std::endl << std::flush;
    CContainer* cc = new CContainer();
    bkmer = new Bkmer( k, bk, kmer );
    for( int i = 0; i < 2; i++ )
    {
        std::cout << "\t\t\tinserting: " << kmer << std::endl << std::flush;
        cc->insert( bkmer );
    }
    test_cc_contain( cc, bkmer );

    char* bad_kmer = ( char* ) "AAATTCCAACCGTGTCTTCTCCATTAA";
    Bkmer* bad_bkmer = new Bkmer( k, bk, bad_kmer );
    test_cc_contain( cc, bad_bkmer );
    delete bkmer;
    delete bad_bkmer;
    delete cc;
    std::cout << "\t\tDone!" << std::endl;
    
    //std::cout << "\t\tUC size: " << uc->size() << std::endl;

    /****************************************************************
     * End compressed container tests
     * *************************************************************/


    /****************************************************************
     * Test kdict
     * *************************************************************/

    std::cout << "\tTesting kdict..." << std::endl << std::flush;
    Kdict* kdict = new Kdict( k, bk );

    std::cout << "\t\tTesting insert and contains..." << std::endl << std::flush;
    for( int i = 0; i < n_test_insert_kmers; i++ )
    {
        kdict->insert( test_insert_kmers[ i ] );
    }

    for( int i = 0; i < n_test_insert_kmers; i++ )
    {
        assert( kdict->contains( test_insert_kmers[ i ] ) );
    }
    std::cout << "\t\t\tFound all good kmers!" << std::endl;

    for( int i = 0; i < n_test_bad_kmers; i++ )
    {
        assert( !kdict->contains( test_bad_kmers[ i ] ) );
    }
    std::cout << "\t\t\tFound no bad kmers!" << std::endl;

    std::cout << "\t\t\tFound " << kdict->size() << " kmers present";
    std::cout << " (correct amount " << ( kdict->size() == n_test_insert_kmers ) << ")" << std::endl;
    
    std::cout << "\t\tTesting remove and contains..." << std::endl << std::flush;
    int n_kmers_after_removal = n_test_insert_kmers;
    for( int i = 0; i < n_test_insert_kmers; i++ )
    {
        if( i % 4 == 0 )
        {
            n_kmers_after_removal--;
            kdict->remove( test_insert_kmers[ i ] );
        }
    }
    std::cout << "\t\t\tFound " << kdict->size() << " kmers present";
    std::cout << " (correct amount " << ( kdict->size() == n_kmers_after_removal ) << ")" << std::endl;
 
    for( int i = 0; i < n_test_insert_kmers; i++ )
    {
        bool contains = kdict->contains( test_insert_kmers[ i ] );
        if( i % 4 == 0 )
        {
            assert( !contains );
        }
        else
        {
            assert( contains );
        }
    }
    std::cout << "\t\t\tCorrect kmers found!" << std::endl;

    // Get kmers
    std::cout << "\t\tTesting kmer seq retrieval..." << std::endl << std::flush;
    int count = 0;
    for( auto bkmer : kdict->get_bkmers() )
    {
        char* seq = bkmer.get_seq();
        std::cout << "\t\t\tseq " << count++ << ": " << seq << std::endl;
        free( seq );
    }

    delete kdict;

    /****************************************************************
     * End kdict tests
     * *************************************************************/


    return 0;
}


