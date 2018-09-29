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
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    root->insert( bkmer );
    delete bkmer;
}

bool Kdict::contains( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    bool res = root->contains( bkmer );
    delete bkmer;
    return res;
}

size_t Kdict::size()
{
    return root->size();
}

void Kdict::remove( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    root->remove( bkmer );
    delete bkmer;
}

void Kdict::clear()
{
    delete root;
    root = new Vertex();
}


