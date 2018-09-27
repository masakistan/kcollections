#include "Bkmer.h"

Bkmer::Bkmer( int k, int bk, char* kmer )
{
    m_k = k;
    m_bk = bk;
    m_bseq = ( uint8_t* ) calloc( bk, sizeof( uint8_t ) );
    serialize_kmer( kmer );
}

void Bkmer::serialize_kmer( char* kmer )
{
    for( int pos = 0; pos < m_k; pos++ )
    {
        switch( kmer[ pos ] )
        {
            case 'a': break;
            case 'A': break;
            case 'c': m_bseq[ pos / 4 ] |= MASK_INSERT[ 0 ][ pos % 4 ]; break;
            case 'C': m_bseq[ pos / 4 ] |= MASK_INSERT[ 0 ][ pos % 4 ]; break;
            case 'g': m_bseq[ pos / 4 ] |= MASK_INSERT[ 1 ][ pos % 4 ]; break;
            case 'G': m_bseq[ pos / 4 ] |= MASK_INSERT[ 1 ][ pos % 4 ]; break;
            case 't': m_bseq[ pos / 4 ] |= MASK_INSERT[ 2 ][ pos % 4 ]; break;
            case 'T': m_bseq[ pos / 4 ] |= MASK_INSERT[ 2 ][ pos % 4 ]; break;
            default: std::ostringstream errMsg;
                     errMsg << "Runtime Error: bad symbol \"" << kmer[ pos ] << "\" when serializing kmer";
                     throw std::runtime_error( errMsg.str() );
        }
    }
}

Bkmer::Bkmer( const Bkmer& other )
{
    int o_bk = other.get_bk();
    uint8_t* o_bseq = other.get_bseq();
    m_bseq = ( uint8_t* ) calloc( o_bk, sizeof( uint8_t ) );
    for( int i = 0; i < o_bk; i++ )
    {
        m_bseq[ i ] = o_bseq[ i ];
    }
}

Bkmer* Bkmer::emit_prefix( int len )
{
    Bkmer* sfpx = get_prefix( len );

    if( sfpx == NULL )
    {
        return sfpx;
    }

    for( int i = 0; i < m_bk; i++ )
    {
        m_bseq[ i ] = m_bseq[ i + len ];
    }
    m_bk -= len;
    m_k = m_k - ( len * 4 );
    resize();

    return sfpx;
}

Bkmer* Bkmer::get_prefix( int len )
{
    if( m_k == 0 )
    {
        return NULL;
    }

    Bkmer* sfpx = new Bkmer( *this );
    sfpx->set_bk( len );
    sfpx->set_k( m_k - len * 4 );
    for( int i = 0; i < len; i++ )
    {
        sfpx->m_bseq[ i ] = m_bseq[ i ];
    }
    sfpx->resize();
    return sfpx;
}

Bkmer* Bkmer::get_suffix( int pos )
{
    if( m_k == 0 )
    {
        return NULL;
    }

    Bkmer* sf = new Bkmer( *this );
    for( int i = pos; i < m_bk; i++ )
    {
        sf->m_bseq[ i - pos ] = m_bseq[ i ];
    }
    sf->set_bk( m_bk - pos );
    sf->set_k( m_k - sf->get_bk() * 4 );
    return sf;
}

void Bkmer::resize()
{
    uint8_t* nbseq = ( uint8_t* ) calloc( m_bk, sizeof( uint8_t ) );
    for( int i = 0; i < m_bk; i ++ )
    {
        nbseq[ i ] = m_bseq[ i ];
    }
    free( m_bseq );
    m_bseq = nbseq;
}

Bkmer::~Bkmer()
{
    free( m_bseq );
}

bool Bkmer::operator<( const Bkmer& other ) const
{
    for( int i = 0; i < BK; i++ )
    {
        if(  unsigned( m_bseq[ i ] ) < unsigned( other.get_bseq()[ i ] ) )
        {
            return true;
        }
        else if( unsigned( m_bseq[ i ] ) > unsigned( other.get_bseq()[ i ] ) )
        {
            return false;
        }
        else
        {
            // continue search
        }
    }
    return false;
}

size_t Bkmer::get_bk() const
{
    return m_bk;
}

size_t Bkmer::get_k() const
{
    return m_k;
}

uint8_t* Bkmer::get_bseq() const
{
    return m_bseq;
}

void Bkmer::set_bseq( uint8_t* bseq )
{
    m_bseq = bseq;
}

bool Bkmer::operator==( const Bkmer& other ) 
{
    if( m_k != other.get_k() || m_bk != other.get_bk() )
    {
        return false;
    }
    
    for( int i = 0; i < m_bk; i ++ )
    {
        if( m_bseq[ i ] != other.get_bseq()[ i ] )
        {
            return false;
        }
    }

    return true;
}

char Bkmer::char_at( int pos )
{
    int bases_to_process, cpos;
    uint8_t pseq;

    for( int i = 0; i < m_bk; i++ )
    {
        if( i == m_bk - 1 )
        {
            bases_to_process = m_k % 4;
        }
        else
        {
            bases_to_process = 4;
        }

        pseq = m_bseq[ i ];
        for( int j = 0; j < bases_to_process; j++ )
        {
            cpos = ( i * 4 ) + j;
            if( cpos == pos )
            {
                return COMP_TO_ASCII[ pseq & 0x3 ];
            }
        }
    }

    return 'N';
}

char* Bkmer::deserialize_seq()
{
    int bases_processed = 0, j, pos, bases_to_process;
    char* kmer = ( char* ) malloc( sizeof( char ) * ( m_k + 1 ) );
    uint8_t tbkmer;
    
    for( int i = 0; i < m_bk; i ++ )
    {
        if( i == m_bk - 1 )
        {
            bases_to_process = m_k % 4;
        }
        else
        {
            bases_to_process = 4;
        }

        tbkmer = m_bseq[ i ];
        for( j = 0; j < bases_to_process; j ++ )
        {
            pos = ( i * 4 ) + j;
            kmer[ pos ] = COMP_TO_ASCII[ tbkmer & 0x3 ];
            tbkmer >>= 2;
        }
    }
    kmer[ m_k ] = '\0';

    return kmer;
}

char* Bkmer::get_seq()
{
    return deserialize_seq();
}


