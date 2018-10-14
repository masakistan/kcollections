#include "Kset.h"

Kset::Kset( int k )
{
    kc = create_kcontainer( k );
    m_k = k;
}

Kset::~Kset()
{
    free_kcontainer( kc );
}

void Kset::clear()
{
    free_kcontainer( kc );
    kc = create_kcontainer( m_k );
}

void Kset::insert( char* kmer )
{
    kcontainer_insert( kc, kmer );
}

bool Kset::contains( char* kmer )
{
    return kcontainer_contains( kc, kmer );
}

uint64_t Kset::size()
{
    return kcontainer_size( kc );
}

void Kset::remove( char* kmer )
{
    kcontainer_remove( kc, kmer );
}


