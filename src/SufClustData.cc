#include "SufClustData.h"
#include "Vertex.h"

SufClustData::SufClustData( Bkmer* sfpx_suffix, bool starts_cluster )
{
    m_sfpx_suffix = std::make_unique< Bkmer >( *sfpx_suffix );
    m_starts_cluster = starts_cluster;
    m_child_vertex = new Vertex();
}

SufClustData::~SufClustData()
{
    delete m_child_vertex;
}

void SufClustData::set_child_vertex( Vertex* vertex )
{
    m_child_vertex = vertex;
}

Bkmer* SufClustData::get_sfpx_suffix() const
{
    return m_sfpx_suffix.get();
}

bool SufClustData::is_cluster_start() const
{
    return m_starts_cluster;
}

Vertex* SufClustData::get_child_vertex()
{
    return m_child_vertex;
}

void SufClustData::set_cluster_start( bool starts_cluster )
{
    m_starts_cluster = starts_cluster;
}


