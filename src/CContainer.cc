#include "CContainer.h"
#include "Vertex.h"

// length of prefix = 1 for 1 uint
// each prefix holds 4 bases
// there are 4 different possible bases
const int CContainer::PREF_SIZE = std::pow( Container::get_prfx_prefix_length() * 4, 4 );

CContainer::CContainer() : Container()
{
    bf = new BloomFilter( 256, 5 );
    m_suf_clust_data = new std::vector< std::unique_ptr< SufClustData > >();
    m_pref = new std::vector< bool >( PREF_SIZE );
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

    std::unique_ptr< Bkmer > sfpx_prefix = sfpx->get_prefix( Container::get_prfx_prefix_length() );
    std::unique_ptr< Bkmer > sfpx_suffix = sfpx->get_suffix( Container::get_prfx_suffix_length() );
    int pref_index = get_index_in_pref( sfpx_prefix.get() );

    if( m_pref->at( pref_index ) )
    {
        int cluster_num = hamming_weight( pref_index );
        int start = rank( cluster_num );
        int pos = start;
        while( pos < m_suf_clust_data->size()
                && ( pos == start || m_suf_clust_data->at( pos )->is_cluster_start() == false ) )
        {
            if( *( m_suf_clust_data->at( pos )->get_sfpx_suffix() ) == *sfpx_suffix )
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

    std::unique_ptr< Bkmer >sfpx_prefix = sfpx->get_prefix( Container::get_prfx_prefix_length() );
    std::unique_ptr< Bkmer >sfpx_suffix = sfpx->get_suffix( Container::get_prfx_suffix_length() );
    int pref_index = get_index_in_pref( sfpx_prefix.get() );

    bool was_pref_index_set = m_pref->at( pref_index );
    m_pref->at( pref_index ) = true;

    int clust_num = hamming_weight( pref_index );
    int clust_pos = rank( clust_num );
    
    // sfpxPrefix was not present
    if( !was_pref_index_set )
    {
        //SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, true );
        //m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
        m_suf_clust_data->emplace( m_suf_clust_data->begin() + clust_pos, new SufClustData( sfpx_suffix.get(), true ) );
        add_to_bloom_filter( sfpx );
        return;
    }

    // sfpxPrefix was already present
    if( was_pref_index_set )
    {
        // if sfpx_suffix starts its cluster...
        bool is_new_sfx_less_than_pos = *sfpx_suffix < *m_suf_clust_data->at( clust_pos )->get_sfpx_suffix();
        if( is_new_sfx_less_than_pos )
        {
            //SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, true );
            //m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
            m_suf_clust_data->emplace( m_suf_clust_data->begin() + clust_pos, new SufClustData( sfpx_suffix.get(), true ) );
            m_suf_clust_data->at( clust_pos + 1 )->set_cluster_start( false );
            add_to_bloom_filter( sfpx );
            return;
        }

        // if sfpx_suffix does not start its cluster...
        clust_pos++;

        // if clust_pos is at end of bit-array
        if(clust_pos >= m_suf_clust_data->size()) {
            //SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, false );
            //m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
            m_suf_clust_data->emplace( m_suf_clust_data->begin() + clust_pos, new SufClustData( sfpx_suffix.get(), false ) );
            add_to_bloom_filter( sfpx );
            return;
        }

        // while sfpx_suffix is greater than previous suffix...
        bool is_new_sfx_greater_than_prev_pos =
            *sfpx_suffix > *m_suf_clust_data->at( clust_pos - 1 )->get_sfpx_suffix();

        while (is_new_sfx_greater_than_prev_pos) {
            // if clust_pos is at end of bit-array OR if next cluster is reached (ie. if sfpx_suffix is greatest in its cluster)...
            if( clust_pos >= m_suf_clust_data->size() || m_suf_clust_data->at( clust_pos )->is_cluster_start() )
            {
                //SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, false );
                //m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos, suf_clust_data );
                m_suf_clust_data->emplace( m_suf_clust_data->begin() + clust_pos, new SufClustData( sfpx_suffix.get(), false ) );
                add_to_bloom_filter( sfpx );
                return;
            }

            clust_pos++;
            is_new_sfx_greater_than_prev_pos = *sfpx_suffix > *m_suf_clust_data->at( clust_pos - 1 )->get_sfpx_suffix();
        }

        // if sfpx_suffix is less than previous suffix...
        //SufClustData* suf_clust_data = new SufClustData( sfpx_suffix, false );
        //m_suf_clust_data->insert( m_suf_clust_data->begin() + clust_pos - 1, suf_clust_data );
        m_suf_clust_data->emplace( m_suf_clust_data->begin() + clust_pos, new SufClustData( sfpx_suffix.get(), false ) );
        add_to_bloom_filter( sfpx );
        return;
    }    
}

/*std::unique_ptr< SufClustData > CContainer::get_suf_clust_data_item( Bkmer* sfpx )
{
    int index_of_child = index_of( sfpx );
    if( index_of_child > -1 )
    {
        return m_suf_clust_data->at( index_of_child );
    }
    return NULL;
}*/

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
    std::vector< Vertex* >* child_vertices = new std::vector< Vertex* >();
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

int CContainer::pref_index_from_hamming_weight( int clust_num )
{
    for( int i = 0; i < m_pref->size(); i++ )
    {
        if( m_pref->at( i ) )
        {
            clust_num--;
            if( clust_num == 0 )
            {
                return i;
            }
        }
    }
    return -1;
}

int CContainer::clust_num_from_rank( int clust_pos )
{
    int clust_num = 0;
    for( ; clust_pos >= 0; clust_pos-- )
    {
        if( m_suf_clust_data->at( clust_pos )->is_cluster_start() )
        {
            clust_num++;
        }
    }
    return clust_num;
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

char* CContainer::prefix_from_clust( int clust_pos )
{
    int clust_num = clust_num_from_rank( clust_pos );
    int pref_index = pref_index_from_hamming_weight( clust_num );
    return index_to_pref( pref_index );
}

/* Converts prefix string to its location on mPref bit-array
 * this should always be a single 8*/
unsigned int CContainer::get_index_in_pref( Bkmer* sfpx_prefix )
{
    assert ( sfpx_prefix->get_bk() == 1 );
    return unsigned( sfpx_prefix->get_bseq()[ 0 ] );
}

char* CContainer::index_to_pref( uint8_t index )
{
    return Bkmer::deserialize_seq( 4, 1, &index );
}

/* Algorithm to return index of the child vertex that holds a given suffix's prefix */
int CContainer::index_of( Bkmer* sfpx )
{
    std::unique_ptr< Bkmer >sfpx_prefix = sfpx->get_prefix( Container::get_prfx_prefix_length() );
    std::unique_ptr< Bkmer >sfpx_suffix = sfpx->get_suffix( Container::get_prfx_prefix_length() );
    int pref_index = get_index_in_pref( sfpx_prefix.get() );

    if( m_pref->at( pref_index ) == true )
    {
        int clust_num = hamming_weight( pref_index );
        int start = rank( clust_num );
        int pos = start;
        while( pos < m_suf_clust_data->size()
                && ( pos == start || !m_suf_clust_data->at( pos )->is_cluster_start() ) )
        {
            if( *( m_suf_clust_data->at( pos )->get_sfpx_suffix() ) == *sfpx_suffix )
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

size_t CContainer::size()
{
    size_t t_size = 0;
    for( std::unique_ptr< SufClustData >& sfc : *m_suf_clust_data )
    {
        t_size += sfc->get_child_vertex()->size();
    }
    return t_size;
}


