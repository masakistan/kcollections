#include "Vertex.h"

Vertex::Vertex()
{
    m_ccs = new std::vector< std::unique_ptr< CContainer > >();
    m_uc = std::make_unique< UContainer >();
    //m_terminal_colors = new std::vector< bool >();
    //m_cardinality = 0;
}

Vertex::~Vertex()
{
    delete m_ccs;
    //delete m_terminal_colors;
}

void Vertex::insert( Bkmer* bkmer )
{
    insert( this, bkmer );
}

void Vertex::insert( Vertex* v, Bkmer* bkmer )
{
    // get prefix position
    int sfpx_len = Container::get_prefix_length();
    std::unique_ptr< Bkmer > sfpx = bkmer->get_prefix( sfpx_len );

    // check if item is in the uncompressed container
    if( m_uc->contains( bkmer ) )
    {
        // can skip, kmer already added
        return;
    }
    
    // check all compressed containers
    for( std::unique_ptr< CContainer >& cc : *m_ccs )
    {
        // check if item is possibly in a compressed container
        if( cc->may_contain( sfpx.get() ) )
        {
            // check if item is actually in compressed container
            if( cc->contains_prefix( sfpx.get() ) )
            {
                sfpx = bkmer->emit_prefix( sfpx_len );
                Vertex* child_vertex = cc->get_child_of( sfpx.get() );
                insert( child_vertex, bkmer );
                return;
            }
            // if it was a false positive
            else
            {
                cc->insert( sfpx.get() );
                Vertex* child_vertex = cc->get_child_of( sfpx.get() );
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

bool Vertex::contains( Vertex* v, Bkmer* bkmer ) const
{
    UContainer* uc = v->get_uc();

    if( uc->contains( bkmer ) )
    {
        return true;
    }

    std::vector< std::unique_ptr< CContainer > >* ccs = v->get_ccs();
    int sfpx_length = Container::get_prefix_length();
    std::unique_ptr< Bkmer > sfpx = bkmer->get_prefix( sfpx_length );

    for( std::unique_ptr< CContainer >& cc : *( ccs ) )
    {
        if( cc->may_contain( sfpx.get() ) && cc->contains_prefix( sfpx.get() ) )
        {
            sfpx = bkmer->emit_prefix( sfpx_length );
            Vertex* child_vertex = cc->get_child_of( sfpx.get() );
            return contains( child_vertex, bkmer );
        }
    }

    return false;
}

void Vertex::burst_uc( Bkmer* bkmer )
{
    std::unique_ptr< UContainer > nuc = std::make_unique< UContainer >(); //new UContainer();
    std::unique_ptr< CContainer > ncc = std::make_unique< CContainer >(); //new CContainer();
    m_uc->insert( bkmer );

    std::unique_ptr< Bkmer > sfpx;

    for( Bkmer uc_bkmer : *( m_uc->get_bkmers() ) )
    {
        sfpx = uc_bkmer.get_prefix( Container::get_prefix_length() );

        // check if compressed container is at capacity
        if( !ncc->is_full() )
        {
            ncc->insert( sfpx.get() );
            uc_bkmer.emit_prefix( Container::get_prefix_length() );
            ncc->get_child_of( sfpx.get() )->insert( &uc_bkmer );
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

    m_uc = std::move( nuc );
    m_ccs->push_back( std::move( ncc ) );
}

size_t Vertex::size()
{
    // check if the kmer exists
    size_t t_size = m_uc->size();

    for( std::unique_ptr< CContainer >& cc : *( m_ccs ) )
    {
        t_size += cc->size();
    }

    return t_size;
}

void Vertex::remove( Bkmer* bkmer )
{
    remove( this, bkmer );
}

void Vertex::remove( Vertex* v, Bkmer* bkmer )
{
    // check if the kmer exists
    UContainer* uc = v->get_uc();

    if( uc->contains( bkmer ) )
    {
        uc->remove( bkmer );
        return;
    }

    std::vector< std::unique_ptr< CContainer > >* ccs = v->get_ccs();
    int sfpx_length = Container::get_prefix_length();
    std::unique_ptr< Bkmer > sfpx = bkmer->get_prefix( sfpx_length );

    for( std::unique_ptr< CContainer >& cc : *( ccs ) )
    {
        if( cc->may_contain( sfpx.get() ) && cc->contains_prefix( sfpx.get() ) )
        {
            sfpx = bkmer->emit_prefix( sfpx_length );
            Vertex* child_vertex = cc->get_child_of( sfpx.get() );
            return remove( child_vertex, bkmer );
        }
    }

    return;
}


