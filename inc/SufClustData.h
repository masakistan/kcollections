#pragma once

#include <stdlib.h>
#include "Bkmer.h"

class Vertex;

class SufClustData
{
    private:
        Bkmer* m_sfpx_suffix;
        bool m_starts_cluster;
        Vertex* m_child_vertex;

    public:
        SufClustData( Bkmer* sfpx_suffix, bool starts_cluster );
        ~SufClustData();
        Bkmer* get_sfpx_suffix() const;
        bool is_cluster_start() const;
        Vertex* get_child_vertex() const;
        void set_cluster_start( bool starts_cluster );
        void set_child_vertex( Vertex* vertex );
        
};


