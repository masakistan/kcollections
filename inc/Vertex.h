#pragma once

#include <vector>
#include "UContainer.h"
#include "CContainer.h"
#include "SufClustData.h"
#include "Bkmer.h"

class Vertex
{
    private:
        int depth;

        // Containers
        UContainer* uc;
        std::vector< CContainer* >* ccs;

        Vertex* get_child_of( CContainer* cc, Bkmer* sfpx );

    public:
        Vertex();
        ~Vertex();
        void insert( Vertex* v, Bkmer* bkmer );
        void insert( Bkmer* bkmer );
        void burst_uc( Bkmer* bkmer );

};


