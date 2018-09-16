#pragma once

#include <vector>
#include "UContainer.h"
#include "CContainer.h"

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
        UContainer* get_uc();
        std::vector< CContainer* >* get_cc();

};


