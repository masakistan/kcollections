#pragma once

#include <stdlib.h>

class Vertex;

class SufClustData
{
    private:
        uint8_t* m_sfpx_suffix;
        bool m_starts_cluster;
        Vertex* m_child_vertex;

    public:
        SufClustData( uint8_t* sfpx_suffix, bool starts_cluster );
        ~SufClustData();
        uint8_t* get_sfpx_suffix() const;
        bool is_cluster_start() const;
        Vertex* get_child_vertex() const;
        void set_cluster_start( bool starts_cluster );
        
};


