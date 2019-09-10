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

void Kdict::add( char* kmer, py::object obj )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kdict" );
    kcontainer_add( kc, kmer, obj );
}

py::object* Kdict::get( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kdict" );
    return kcontainer_get( kc, kmer );
}

bool Kdict::contains( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kdict" );
    return kcontainer_contains( kc, kmer );
}

uint64_t Kdict::size()
{
    return kcontainer_size( kc );
}

void Kdict::remove( char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kdict" );
    kcontainer_remove( kc, kmer );
}

