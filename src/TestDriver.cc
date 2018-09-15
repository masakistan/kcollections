#include <iostream>
#include "Kdict.h"

int main()
{
    // Testing Parameters
    int k = 27;
    int bk = calc_bk( k );

    std::cout << "Testing kmer dictionary implementation!" << std::endl;
    Kdict* kd = new Kdict( k, bk );

    /****************************************************************
     * Test serialization and deserialization capabilities
     * *************************************************************/
    char* kmer = ( char* ) "AAATTCCAACCGTGTCTTCTCCATTAG";
    Bkmer* bkmer = serialize_kmer( kmer, k, bk );
    for( int i = 0; i < bk; i++ )
    {
        std::cout << unsigned( bkmer->bseq[ i ] ) << std::endl;
    }
    std::cout << "deserialized binary kmer: " << deserialize_bkmer( bkmer, k, bk ) << std::endl;
    /****************************************************************
     * End serialization tests
     * *************************************************************/


    return 0;
}
