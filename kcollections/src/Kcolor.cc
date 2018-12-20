#include "Kcolor.h"


Kcolor::Kcolor( int k )
{
    kc = create_kcontainer( k );
    m_k = k;
}

Kcolor::~Kcolor()
{
    free_kcontainer( kc );
}

void Kcolor::clear()
{
    free_kcontainer( kc );
    kc = create_kcontainer( m_k );
}

void Kcolor::insert(const char* kmer, uint32_t color)
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcolor" );
    kcontainer_add( kc, kmer, color );
}

uint32_t* Kcolor::get( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcolor" );
    return kcontainer_get( kc, kmer );
}

bool Kcolor::contains( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcolor" );
    return kcontainer_contains( kc, kmer );
}

uint64_t Kcolor::size()
{
    return kcontainer_size( kc );
}

void Kcolor::remove( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kcolor" );
    kcontainer_remove( kc, kmer );
}

void Kcolor::add_seq(const char* seq, uint32_t length, uint32_t color)
{
    kcontainer_add_seq( kc, seq, length, color );
}

void Kcolor::parallel_add_seq(const char* seq, uint32_t length, uint32_t color) {
    parallel_kcontainer_add_seq(kc, seq, length, color);
}
