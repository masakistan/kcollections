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
    insert( this, bkmer );
}

void Vertex::insert( Vertex* v, Bkmer* bkmer )
{
    // get prefix position
    int sfpx_len = Container::get_prefix_length();
    Bkmer* sfpx = bkmer->get_prefix( sfpx_len );

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
            if( cc->contains_prefix( sfpx ) )
            {
                sfpx = bkmer->emit_prefix( sfpx_len );
                Vertex* child_vertex = cc->get_child_of( sfpx );
                insert( child_vertex, bkmer );
                return;
            }
            // if it was a false positive
            else
            {
                cc->insert( sfpx );
                Vertex* child_vertex = cc->get_child_of( sfpx );
                bkmer->emit_prefix( sfpx_len );
                insert( child_vertex, bkmer );
            }
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
        uc->insert( bkmer );
    }
}

void Vertex::burst_uc( Bkmer* bkmer )
{
    UContainer* nuc = new UContainer();
    CContainer* ncc = new CContainer();
    uc->insert( bkmer );

    Bkmer* sfpx;

    for( Bkmer* uc_bkmer : *( uc->get_bkmers() ) )
    {
        sfpx = uc_bkmer->get_prefix( Container::get_prefix_length() );

        // check if compressed container is at capacity
        if( !ncc->is_full() )
        {
            uc_bkmer->emit_prefix( Container::get_prefix_length() );
            ncc->get_child_of( sfpx )->insert( uc_bkmer );
        }
        else if( !nuc->is_full() )
        {
            nuc->insert( uc_bkmer );
        }
        else
        {
            /// this is an error state
        }
    }

    uc = nuc;
    ccs->push_back( ncc );
}


