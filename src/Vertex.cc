#include "Vertex.h"

Vertex::Vertex()
{
    ccs = new std::vector< CContainer*>();
    uc = new UContainer();
}

Vertex::~Vertex()
{
    delete uc;
    delete ccs;
}

void Vertex::insert( Bkmer* bkmer )
{
    // get prefix position
    int pre_pos = SFPX_LEN;
    uint8_t* sfpx = ( uint8_t* ) calloc( SFPX_LEN, sizeof( uint8_t ) );
    for( int i = 0; i < SFPX_LEN; i++ )
    {
        sfpx[ i ] = bkmer->get_bseq()[ i ];
    }

    // check if item is in the uncompressed container
    if( uc->contains_kmer( bkmer ) )
    {
        // can skip, kmer already added
        return;
    }
    
    // check all compressed containers
    for( CContainer* cc : *ccs )
    {
        // check if item is possibly in a compressed container
        if( cc->may_contain( sfpx ) )
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
    
    if( uc->is_full() )
    {
        // burst container if necessary
        burst_uc( bkmer );
    }
    else
    {
        // add to uncompressed container
        uc->insert( *bkmer );
    }
}

void Vertex::burst_uc( Bkmer* bkmer )
{
    UContainer* nuc = new UContainer();
    CContainer* ncc = new CContainer();

    uint8_t* sfpx;

    for( Bkmer uc_bkmer : *( uc->get_bkmers() ) )
    {
        sfpx = uc_bkmer.get_prefix( Container::get_prefix_length() );

    }

    uc = nuc;
    ccs->push_back( ncc );
}


