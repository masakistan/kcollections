#include "Kdict.h"

Kdict::Kdict( int k, int bk ) : m_k( k ), m_bk( bk )
{
    root = new Vertex();
    K = m_k;
    BK = bk;
}

Kdict::~Kdict()
{
    delete root;
}

void Kdict::insert( char* kmer )
{
    Bkmer* bkmer = serialize_kmer( kmer, m_k, m_bk );
    root->insert( bkmer );
}


