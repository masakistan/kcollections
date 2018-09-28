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
    std::cout << "inserting: " << kmer << std::endl;
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    root->insert( bkmer );
}

bool Kdict::contains( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    return root->contains( bkmer );
    delete bkmer;
}


