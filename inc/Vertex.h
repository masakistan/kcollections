#pragma once

#include <vector>
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
        void insert( Bkmer* bkmer );
        void burst_uc( Bkmer* bkmer );

};


