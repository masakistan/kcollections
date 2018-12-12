#include "Kcounter.h"


Kcounter::Kcounter( int k )
{
    kc = create_kcontainer( k );
    m_k = k;
}

Kcounter::~Kcounter()
{
    free_kcontainer( kc );
}

void Kcounter::clear()
{
    free_kcontainer( kc );
    kc = create_kcontainer( m_k );
}

void Kcounter::insert( char* kmer, int count )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcounter" );
    kcontainer_add( kc, kmer, count );
}

int Kcounter::get( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcounter" );
    return kcontainer_get( kc, kmer );
}

bool Kcounter::contains( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcounter" );
    return kcontainer_contains( kc, kmer );
}

uint64_t Kcounter::size()
{
    return kcontainer_size( kc );
}

void Kcounter::remove( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcounter" );
    kcontainer_remove( kc, kmer );
}

void Kcounter::add_seq( char* seq, uint32_t length )
{
    kcontainer_add_seq( kc, seq, length );
}
