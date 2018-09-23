#include "SufClustData.h"
#include "Vertex.h"

SufClustData::SufClustData( Bkmer* sfpx_suffix, bool starts_cluster )
{
    m_sfpx_suffix = sfpx_suffix;
    m_starts_cluster = starts_cluster;
    //m_child_vertex = NULL;
    m_child_vertex = new Vertex();
}

SufClustData::~SufClustData()
{
    free( m_sfpx_suffix );
    delete m_child_vertex;
}

void SufClustData::set_child_vertex( Vertex* vertex )
{
    m_child_vertex = vertex;
}

Bkmer* SufClustData::get_sfpx_suffix() const
{
    return m_sfpx_suffix;
}

bool SufClustData::is_cluster_start() const
{
    return m_starts_cluster;
}

Vertex* SufClustData::get_child_vertex() const
{
    return m_child_vertex;
}

void SufClustData::set_cluster_start( bool starts_cluster )
{
    m_starts_cluster = starts_cluster;
}


