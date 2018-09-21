#include "Bkmer.h"

Bkmer::Bkmer( int k, int bk )
{
    m_bseq = ( uint8_t* ) calloc( BK, sizeof( uint8_t ) );
    m_k = K;
    m_bk = BK;
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

bool Bkmer::operator<( Bkmer& other ) const
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


