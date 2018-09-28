#include "Vertex.h"

Vertex::Vertex()
{
    m_ccs = new std::vector< CContainer*>();
    m_uc = new UContainer();
    //m_terminal_colors = new std::vector< bool >();
    //m_cardinality = 0;
}

Vertex::~Vertex()
{
    delete m_uc;
    delete m_ccs;
    //delete m_terminal_colors;
}

void Vertex::insert( Bkmer* bkmer )
{
    insert( this, bkmer );
}

void Vertex::insert( Vertex* v, Bkmer* bkmer )
{
    //std::cout << "\ttrying to insert: " << bkmer->get_seq() << std::endl;
    // get prefix position
    int sfpx_len = Container::get_prefix_length();
    Bkmer* sfpx = bkmer->get_prefix( sfpx_len );

    // check if item is in the uncompressed container
    if( m_uc->contains( bkmer ) )
    {
        // can skip, kmer already added
        return;
    }
    
    // check all compressed containers
    for( CContainer* cc : *m_ccs )
    {
        //std::cout << "\tchecking: " << sfpx->get_seq() << std::endl;
        // check if item is possibly in a compressed container
        if( cc->may_contain( sfpx ) )
        {
            //std::cout << "\t\tcc may contain " << std::endl;
            // check if item is actually in compressed container
            if( cc->contains_prefix( sfpx ) )
            {
                //std::cout << "\t\t\tcc contains" << std::endl;
                sfpx = bkmer->emit_prefix( sfpx_len );
                Vertex* child_vertex = cc->get_child_of( sfpx );
                insert( child_vertex, bkmer );
                return;
            }
            // if it was a false positive
            else
            {
                ////std::cout << "\t\t\tcc false positive" << std::endl;
                cc->insert( sfpx );
                Vertex* child_vertex = cc->get_child_of( sfpx );
                bkmer->emit_prefix( sfpx_len );
                insert( child_vertex, bkmer );
            }
        }
    }
    
    if( m_uc->is_full() )
    {
        // burst container if necessary
        burst_uc( bkmer );
    }
    else
    {
        // add to uncompressed container
        m_uc->insert( bkmer );
    }
}

bool Vertex::contains( Bkmer* bkmer )
{
    return contains( this, bkmer );
}

bool Vertex::contains( Vertex* v, Bkmer* bkmer )
{
    UContainer* uc = v->get_uc();

    if( uc->contains( bkmer ) )
    {
        return true;
    }

    std::vector< CContainer* >* ccs = v->get_ccs();
    int sfpx_length = Container::get_prefix_length();
    Bkmer* sfpx = bkmer->get_prefix( sfpx_length );

    for( CContainer* cc : *( ccs ) )
    {
        if( cc->may_contain( sfpx ) && cc->contains_prefix( sfpx ) )
        {
            sfpx = bkmer->emit_prefix( sfpx_length );
            Vertex* child_vertex = cc->get_child_of( sfpx );
            return contains( child_vertex, bkmer );
        }
    }

    return false;
}

void Vertex::burst_uc( Bkmer* bkmer )
{
    //std::cout << "bursting " << bkmer->get_seq() << std::endl;
    UContainer* nuc = new UContainer();
    CContainer* ncc = new CContainer();
    m_uc->insert( bkmer );

    Bkmer* sfpx;

    for( Bkmer uc_bkmer : *( m_uc->get_bkmers() ) )
    {
        //std::cout << "\tbursting processing bkmer: " << uc_bkmer.get_seq() << std::endl;
        sfpx = uc_bkmer.get_prefix( Container::get_prefix_length() );

        // check if compressed container is at capacity
        if( !ncc->is_full() )
        {
            ncc->insert( sfpx );
            uc_bkmer.emit_prefix( Container::get_prefix_length() );
            //std::cout << "\tnew uc bkmer: " << uc_bkmer.get_seq() << std::endl;
            ncc->get_child_of( sfpx )->insert( &uc_bkmer );
        }
        else if( !nuc->is_full() )
        {
            nuc->insert( &uc_bkmer );
        }
        else
        {
            /// this is an error state
        }
    }

    m_uc = nuc;
    m_ccs->push_back( ncc );
    //std::cout << "Done bursting!" << std::endl;
}


