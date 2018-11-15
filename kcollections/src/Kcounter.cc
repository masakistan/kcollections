#include "Kcounter.h"

Kcounter::Kcounter( int k ) : Kdict( k )
{
}

Kcounter::~Kcounter()
{
    free_kcontainer( kc );
}

void Kcounter::insert( char* kmer, int count )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcounter" );
    kcontainer_insert( kc, kmer, count );
}

int Kcounter::get( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcounter" );
    return kcontainer_get( kc, kmer );
}
