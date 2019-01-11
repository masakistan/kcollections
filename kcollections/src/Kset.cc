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

/*void Kset::add( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kset" );
    kcontainer_add( kc, kmer );
}*/

PgData* Kset::get( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kdict" );
    return kcontainer_get( kc, kmer );
}

bool Kset::contains( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kset" );
    return kcontainer_contains( kc, kmer );
}

uint64_t Kset::size()
{
    return kcontainer_size( kc );
}

void Kset::remove( const char* kmer )
{
    CHECK_KMER_LENGTH( kmer, m_k, "Kset" );
    kcontainer_remove( kc, kmer );
}

/*void Kset::add_seq(const char* seq, uint32_t length)
{
    kcontainer_add_seq(kc, seq, length);
}*/

void Kset::parallel_add_seq(char* seq, uint32_t length, uint16_t gidx) {
    if(gidx + 1 > chromoBoundaries.size()) {
        chromoBoundaries.push_back(std::vector<uint32_t>());
        chromoBoundaries[gidx].push_back(0);
    }

    parallel_kcontainer_add_seq(kc, seq, length, gidx, chromoBoundaries[gidx].back());
    chromoBoundaries[gidx].push_back(length);
}
