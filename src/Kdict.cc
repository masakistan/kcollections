#include "Kdict.h"

Kdict::Kdict( int k, int bk ) : m_k( k ), m_bk( bk )
{
}

void Kdict::insert( char* kmer )
{
    uint8_t* bkmer = serialize_kmer( kmer, m_k, m_bk );
}


