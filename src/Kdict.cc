#incldue "Kdict.h"

Kdict::Kdict( int k )
{
    m_k = k;
    
    // set binary k size
    m_bk = k / 4;
    if( k % 4 > 0 )
    {
        m_k++;
    }
}

void Kdict::insert( char* kmer )
{
    uint8_t* bkmer = serializeKmer( kmer, m_k, m_bk );
}


