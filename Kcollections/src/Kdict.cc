#include "Kdict.h"


Kdict::Kdict( int k )
{
    kc = create_kcontainer( k );
    m_k = k;
}

Kdict::~Kdict()
{
    free_kcontainer( kc );
}

void Kdict::clear()
{
    free_kcontainer( kc );
    kc = create_kcontainer( m_k );
}

void Kdict::insert( char* kmer, py::handle* obj )
{
    kcontainer_insert( kc, kmer, obj );
}

py::handle* Kdict::get( char* kmer )
{
    return kcontainer_get( kc, kmer );
}

bool Kdict::contains( char* kmer )
{
    return kcontainer_contains( kc, kmer );
}

uint64_t Kdict::size()
{
    return kcontainer_size( kc );
}

void Kdict::remove( char* kmer )
{
    kcontainer_remove( kc, kmer );
}


