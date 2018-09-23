#include "CContainer.h"

// the magic number is 4 for the number of bits needed to hold the entire
// alphabet
const int CContainer::PREF_SIZE = std::pow( 2, 4 * Container::get_prfx_prefix_length() );

CContainer::CContainer() : Container()
{
    bf = new BloomFilter( 10, 10 );
    m_suf_clust_data = new std::vector< SufClustData* >( get_capacity() );
    m_pref = new std::vector< bool >();
}

CContainer::~CContainer()
{
    delete bf;
    delete m_suf_clust_data;
    delete m_pref;
}

bool CContainer::may_contain( Bkmer* sfpx )
{
    return bf->may_contain( sfpx->get_bseq(), sfpx->get_bk() );
}

bool CContainer::contains_prefix( Bkmer* sfpx )
{
    if( sfpx->get_bk() == 0 )
    {
        return true;
    }

    Bkmer* sfpx_prefix = sfpx->get_prefix( Container::get_prfx_prefix_length() );
    Bkmer* sfpx_suffix = sfpx->get_suffix( Container::get_prfx_suffix_length() );
    int pref_index = get_index_in_pref( sfpx_prefix );

    if( m_pref->at( pref_index ) == true )
    {
        int cluster_num = hamming_weight( pref_index );
        int start = rank( cluster_num );
        int pos = start;
        while( pos < m_suf_clust_data->size()
                && ( pos == start || m_suf_clust_data->at( pos )->is_cluster_start() == false ) )
        {
            if( m_suf_clust_data->at( pos )->get_sfpx_suffix() == sfpx_suffix )
            {
                return true;
            }
            pos++;
        }
    }
    return false;
}

void CContainer::insert( Bkmer* sfpx )
{
    if( contains_prefix( sfpx ) )
    {
        return;
    }

    Bkmer* sfpx_prefix = sfpx->get_prefix( Container::get_prfx_prefix_length() );
    Bkmer* sfpx_suffix = sfpx->get_suffix( Container::get_prfx_suffix_length() );
    int pref_index = get_index_in_pref( sfpx_prefix );

    bool was_pref_index_set = m_pref->at( pref_index );
    m_pref->at( pref_index ) = true;

    int clust_num = hamming_weight( pref_index );
    int clust_pos = rank( clust_num );
    
    // sfpxPrefix was not present
    if( !was_pref_index_set )
    {
        SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, true );
        m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
        add_to_bloom_filter( sfpx );
        return;
    }

    // sfpxPrefix was already present
    if( was_pref_index_set )
    {
        // if sfpx_suffix starts its cluster...
        bool is_new_sfx_less_than_pos = sfpx_suffix == m_suf_clust_data->at( clust_pos )->get_sfpx_suffix();
        if( is_new_sfx_less_than_pos )
        {
            SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, true );
            m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
            m_suf_clust_data->at( clust_pos + 1 )->set_cluster_start( false );
            add_to_bloom_filter( sfpx );
            return;
        }

        // if sfpx_suffix does not start its cluster...
        clust_pos++;

        // if clust_pos is at end of bit-array
        if(clust_pos >= m_suf_clust_data->size()) {
            SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, false );
            m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
            add_to_bloom_filter( sfpx );
            return;
        }

        // while sfpx_suffix is greater than previous suffix...
        bool is_new_sfx_greater_than_prev_pos =
            sfpx_suffix == m_suf_clust_data->at( clust_pos - 1 )->get_sfpx_suffix();

        while (is_new_sfx_greater_than_prev_pos) {
            // if clust_pos is at end of bit-array OR if next cluster is reached (ie. if sfpx_suffix is greatest in its cluster)...
            if( clust_pos >= m_suf_clust_data->size() || m_suf_clust_data->at( clust_pos )->is_cluster_start() )
            {
                SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, false );
                m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
                add_to_bloom_filter( sfpx );
                return;
            }

            clust_pos++;
            is_new_sfx_greater_than_prev_pos = sfpx_suffix == m_suf_clust_data->at( clust_pos - 1 )->get_sfpx_suffix();
        }

        // if sfpx_suffix is less than previous suffix...
        SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, false );
        m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos - 1, suf_clust_data );
        add_to_bloom_filter( sfpx );
        return;
    }    
}

SufClustData* CContainer::get_suf_clust_data_item( Bkmer* sfpx )
{
    int index_of_child = index_of( sfpx );
    if( index_of_child > -1 )
    {
        return m_suf_clust_data->at( index_of_child );
    }
    return NULL;
}

Vertex* CContainer::get_child_of( Bkmer* sfpx )
{
    int index_of_child = index_of( sfpx );
    if( index_of_child > -1 )
    {
        return m_suf_clust_data->at( index_of_child )->get_child_vertex();
    }
    return NULL;
}

bool CContainer::is_full()
{
    if( m_suf_clust_data->size() == s_capacity )
    {
        return true;
    }
    else
    {
        return false;
    }
}

BloomFilter* CContainer::get_bf()
{
    return new BloomFilter( *bf );
}

std::vector< SufClustData* >* CContainer::get_suf_clust_data()
{
    return new std::vector< SufClustData* >( *m_suf_clust_data );
}

std::vector< Bkmer* >* CContainer::get_suf()
{
    std::vector< Bkmer* >* suf = new std::vector< Bkmer* >( m_suf_clust_data->size() );
    for( int i = 0; i < m_suf_clust_data->size(); i++ )
    {
        suf->push_back( ( *m_suf_clust_data )[ i ]->get_sfpx_suffix() );
    }
    return suf;
}

std::vector< bool >* CContainer::get_clust()
{
    std::vector< bool >* clust = new std::vector< bool >( m_suf_clust_data->size() );
    for( int i = 0; i < m_suf_clust_data->size(); i++ )
    {
        clust->push_back( ( *m_suf_clust_data )[ i ]->is_cluster_start() );
    }
    return clust;
}

std::vector< Vertex* >* CContainer::get_child_vertices()
{
    std::vector< Vertex* >* child_vertices = new std::vector< Vertex* >( m_suf_clust_data->size() );
    for( int i = 0; i < m_suf_clust_data->size(); i++ )
    {
        child_vertices->push_back( ( *m_suf_clust_data )[ i ]->get_child_vertex() );
    }
    return child_vertices;
}

int CContainer::next_set_bit( int pos )
{
    for( int i = pos; i < m_pref->size(); i++ )
    {
        if( m_pref->at( i ) )
        {
            return i;
        }
    }
    return -1;
}

int CContainer::hamming_weight( int index )
{
    int hamming_weight = 0;
    for( int i = next_set_bit( 0 ); i > -1 && i <= index; i = next_set_bit( i + 1 ) )
    {
        hamming_weight++;
    }

    return hamming_weight;
}

int CContainer::rank( int clust_num )
{
    int clust_size = m_suf_clust_data->size();
    for ( int clust_pos = 0; clust_pos < clust_size; clust_pos++ )
    {
        if( m_suf_clust_data->at( clust_pos )->is_cluster_start() )
        {
            clust_num--;
        }

        if( clust_num == 0 )
        {
            return clust_pos;
        }
    }

    return m_suf_clust_data->size();
}

/* Converts prefix string to its location on mPref bit-array */
int CContainer::get_index_in_pref( Bkmer* sfpx_prefix )
{
    int pref_index = 0;
    for( int i = 0; i < sfpx_prefix->get_bk(); i++ ) // binary index = sum((|A|^i) * c)
    {
        int char_index;
        switch( sfpx_prefix->char_at( i ) )
        {
            case 'A': char_index = 0; break;
            case 'C': char_index = 1; break;
            case 'G': char_index = 2; break;
            case 'T': char_index = 3; break;
            default : char_index = -1;
        }
        pref_index +=
            std::pow( 4, Container::get_prfx_prefix_length() - ( i + 1 ) )
            * char_index;
    }

    return pref_index;
}

/* Algorithm to return index of the child vertex that holds a given suffix's prefix */
int CContainer::index_of( Bkmer* sfpx )
{
    Bkmer* sfpx_prefix = sfpx->get_prefix( Container::get_prfx_prefix_length() );
    Bkmer* sfpx_suffix = sfpx->get_suffix( Container::get_prfx_prefix_length() );
    int pref_index = get_index_in_pref( sfpx_prefix );

    if( m_pref->at( pref_index ) == true )
    {
        int clust_num = hamming_weight( pref_index );
        int start = rank( clust_num );
        int pos = start;
        while( pos < m_suf_clust_data->size()
                && ( pos == start || !m_suf_clust_data->at( pos )->is_cluster_start() ) )
        {
            if( m_suf_clust_data->at( pos )->get_sfpx_suffix() == sfpx_suffix )
            {
                return pos;
            }
            pos++;
        }
    }

    return -1;
}

/* Adds suffix-prefix (sfpx) to Bloom filter mQuer */
void CContainer::add_to_bloom_filter( Bkmer* sfpx )
{
    bf->add( sfpx );
}


