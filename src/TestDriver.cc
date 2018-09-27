#include <iostream>
#include "Kdict.h"

int main()
{
    // Testing Parameters
    int k = 27;
    int bk = calc_bk( k );

    std::cout << "Testing kmer dictionary implementation!" << std::endl;
    std::cout << "***************************************" << std::endl;
    Kdict* kd = new Kdict( k, bk );


    /****************************************************************
     * Test serialization and deserialization capabilities
     * *************************************************************/
    char* kmer = ( char* ) "AAATTCCAACCGTGTCTTCTCCATTAG";
    std::cout << "\tTesting kmer serialization/deserialization..." << std::endl;
    std::cout << "\t\tkmer: " << kmer << std::endl;

    Bkmer* bkmer = new Bkmer( k, bk, kmer );
    for( int i = 0; i < bk; i++ )
    {
        std::cout << "\t\t\t" << unsigned( bkmer->get_bseq()[ i ] ) << std::endl;
    }

    std::cout << "\t\tdeserialized binary kmer: " << bkmer->get_seq() << std::endl;
    /****************************************************************
     * End serialization tests
     * *************************************************************/


    /****************************************************************
     * Test Uncompressed Container
     * *************************************************************/
    std::cout << "\tTesting uncompressed container..." << std::endl;
    UContainer* uc = new UContainer();

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
    /****************************************************************
     * End uncompressed container tests
     * *************************************************************/


    /****************************************************************
     * Test Compressed Container
     * *************************************************************/
    CContainer* cc = new CContainer();
    std::cout << "\tTesting compressed container..." << std::endl;
    Bkmer* bad_bkmer = new Bkmer( k, bk, "AAATTCCAACCGTGTCTTCTCCATTAA" );

    for( int i = 0; i < 2; i++ )
    {
        cc->insert( bkmer );
    }
    
    std::cout << "\t\tkmer may contain " << bkmer->get_seq();
    if( cc->may_contain( bkmer ) )
    {
        std::cout << " True!" << std::endl;
    }
    else
    {
        std::cout << " False!" << std::endl;
    }

    std::cout << "\t\tkmer may contain " << bad_bkmer->get_seq();
    if( cc->may_contain( bad_bkmer ) )
    {
        std::cout << " True!" << std::endl;
    }
    else
    {
        std::cout << " False!" << std::endl;
    }

    
    //std::cout << "\t\tUC size: " << uc->size() << std::endl;

    /****************************************************************
     * End compressed container tests
     * *************************************************************/

    return 0;
}
