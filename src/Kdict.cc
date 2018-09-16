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
    insert( root, bkmer );
}

void Kdict::insert( Vertex* v, Bkmer* bkmer )
{
    // check if item is in the uncompressed container
    UContainer* uc = v->get_uc();
    if( uc->contains_kmer( bkmer ) )
    {
        // can skip, kmer already added
        return;
    }
    
    // check all compressed containers
    std::vector< CContainer* >* ccs = v->get_cc();
    for( CContainer* cc : *ccs )
    {
        // check if item is possibly in a compressed container
        if( cc->may_contain( ) )
        {
            // check if item is actually in compressed container
            /*if( cc.contains_prefix( ) )
            {
            }
            else
            {
            }*/
        }
    }
    
    // add to uncompressed container
        // burst container if necessary
}


