#pragma once

#include <vector>
#include "Container.h"
#include "UContainer.h"
#include "CContainer.h"
#include "SufClustData.h"
#include "Bkmer.h"

class Vertex
{
    private:
        int depth;

        // Containers
        UContainer* m_uc;
        std::vector< CContainer* >* m_ccs;
        std::vector< bool >* m_terminal_colors;
        int m_cardinality;

    public:
        Vertex();
        ~Vertex();

        std::vector< bool >* get_terminal_colors() { return m_terminal_colors; }
        std::vector< CContainer* >* get_ccs() { return m_ccs; }
        UContainer* get_uc() { return m_uc; }
        int num_terminal_colors() { return m_cardinality; }


        void insert( Vertex* v, Bkmer* bkmer );
        void insert( Bkmer* bkmer );
        bool contains( Bkmer* bkmer );
        bool contains_sequence( Bkmer* bkmer );
        bool contains( Vertex* v, Bkmer* bkmer );

        void burst_uc( Bkmer* bkmer );

};


