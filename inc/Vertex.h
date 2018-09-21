#pragma once

#include <vector>
#include "Container.h"
#include "UContainer.h"
#include "CContainer.h"
#include "Bkmer.h"

class Vertex
{
    private:
        int depth;

        // Containers
        UContainer* uc;
        std::vector< CContainer* >* ccs;

    public:
        Vertex();
        ~Vertex();
        void insert( Vertex* v, Bkmer* bkmer );
        void insert( Bkmer* bkmer );
        void burst_uc( Bkmer* bkmer );

};


