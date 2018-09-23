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
    Bkmer* bkmer = new Bkmer( k, bk, kmer );
    for( int i = 0; i < bk; i++ )
    {
        std::cout << "\t" << unsigned( bkmer->get_bseq()[ i ] ) << std::endl;
    }
    std::cout << "\tdeserialized binary kmer: " << bkmer->deserialize_seq() << std::endl;
    /****************************************************************
     * End serialization tests
     * *************************************************************/

    /****************************************************************
     * Test Uncompressed Container
     * *************************************************************/
    UContainer* uc = new UContainer();
    uc->insert( bkmer );
    if( uc->contains_kmer( bkmer ) )
    {
        std::cout << "\tUC kmer contains success!" << std::endl;
    }
    else
    {
        std::cout << "\tUC kmer contains failed!" << std::endl;
    }
    /****************************************************************
     * End uncompressed container tests
     * *************************************************************/



    return 0;
}
