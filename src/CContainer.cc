#include "CContainer.h"

CContainer::CContainer() : Container()
{
    bf = new BloomFilter( 10, 10 );
    m_suf_clust_data = new std::vector< SufClustData* >( get_capacity() );
}

CContainer::~CContainer()
{
    delete bf;
    delete m_suf_clust_data;
}

bool CContainer::may_contain( uint8_t* sfpx )
{
    return bf->may_contain( sfpx, SFPX_LEN );
}

bool CContainer::contains_prefix( uint8_t* sfpx )
{
    return false;
}

void CContainer::insert( uint8_t* sfpx )
{
    bf->add( sfpx, SFPX_LEN );
}

Vertex* CContainer::get_child_of( uint8_t* spfx )
{
    return NULL;
}



