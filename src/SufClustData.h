#include "SufClustData.h"

SufClustData::SufClustData( uint8_t* sfpx_suffix, bool starts_cluster )
{
    m_sfpx_suffix = sfpx_suffix;
    m_starts_cluster = starts_cluster;
    m_child_vertex = new Vertex();
}

SufClustData::~SufClustData()
{
    free( m_sfpx_suffix );
    delete m_child_vertex;
}




