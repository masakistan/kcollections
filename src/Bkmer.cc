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

    m_bk = other.get_bk();
    m_k = other.get_k();
    m_bseq = ( uint8_t* ) calloc( m_bk, sizeof( uint8_t ) );
    
    uint8_t* o_bseq = other.get_bseq();
    for( int i = 0; i < m_bk; i++ )
    {
        m_bseq[ i ] = o_bseq[ i ];
    }
}

std::unique_ptr< Bkmer > Bkmer::emit_prefix( int len )
{
    std::unique_ptr< Bkmer > sfpx = get_prefix( len );

    if( sfpx == NULL )
    {
        return sfpx;
    }

    for( int i = 0; i < m_bk - len; i++ )
    {
        m_bseq[ i ] = m_bseq[ i + len ];
    }
    m_bk -= len;
    m_k = m_k - ( len * 4 );
    resize();

    return sfpx;
}

std::unique_ptr< Bkmer > Bkmer::get_prefix( int len )
{
    if( m_k == 0 )
    {
        return NULL;
    }

    std::unique_ptr< Bkmer > sfpx = std::make_unique< Bkmer >( *this );
    sfpx->set_bk( len );
    sfpx->set_k( len * 4 );
    for( int i = 0; i < len; i++ )
    {
        sfpx->m_bseq[ i ] = m_bseq[ i ];
    }
    sfpx->resize();
    return sfpx;
}

std::unique_ptr< Bkmer > Bkmer::get_suffix( int pos )
{
    if( m_k == 0 )
    {
        return NULL;
    }

    std::unique_ptr< Bkmer > sf = std::make_unique< Bkmer >( *this );
    int bytes = m_bk - pos;
    int removed = m_bk - bytes;
    for( int i = 0; i < bytes; i++ )
    {
        sf->m_bseq[ 0 ] = m_bseq[ i + pos ];
    }
    sf->set_bk( bytes );
    sf->set_k( m_k - ( removed * 4 ) );
    return sf;
}

void Bkmer::resize()
{
    uint8_t* nbseq = ( uint8_t* ) calloc( m_bk, sizeof( uint8_t ) );
    for( int i = 0; i < m_bk; i++ )
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
    for( int i = 0; i < m_bk; i++ )
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
    free( m_bseq );
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

char* Bkmer::deserialize_seq( int k, int bk, uint8_t* bseq )
{
    int bases_processed = 0, j, pos, bases_in_byte, bases_to_process = k;
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );
    uint8_t tbkmer;
    
    for( int i = 0; i < bk; i ++ )
    {
        if( bases_to_process > 4 )
        {
            bases_in_byte = 4;
        }
        else
        {
            bases_in_byte = bases_to_process;
        }

        tbkmer = bseq[ i ];
        for( j = 0; j < bases_in_byte; j++ )
        {
            pos = ( i * 4 ) + j;
            kmer[ pos ] = COMP_TO_ASCII[ tbkmer & 0x3 ];
            tbkmer >>= 2;
        }
        bases_to_process -= bases_in_byte;
    }
    kmer[ k ] = '\0';

    return kmer;
}

char* Bkmer::get_seq() const
{
    return deserialize_seq( m_k, m_bk, m_bseq );
}


