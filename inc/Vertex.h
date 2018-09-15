#pragma once

#include "UContainer.h"
#include "CContainer.h"

class Vertex
{
    private:
        int depth;

        // Containers

        UContainer* uc;
        //std::vector< CompressedContainer* >* ccs;

    public:
        Vertex();
        ~Vertex();
        UContainer* get_uc();

};


