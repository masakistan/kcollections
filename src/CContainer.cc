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

bool CContainer::may_contain( Bkmer* sfpx )
{
    return bf->may_contain( sfpx->get_bseq(), sfpx->get_bk() );
}

bool CContainer::contains_prefix( Bkmer* sfpx )
{
    return false;
}

void CContainer::insert( Bkmer* sfpx )
{
    bf->add( sfpx );
}

Vertex* CContainer::get_child_of( Bkmer* spfx )
{
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


