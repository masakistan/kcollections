#pragma once

#include "UncompressedContainer.h"
#include "CompressedContainer.h"

class Vertex
{
    private:
        int depth;

        // Containers

        UncompressedContainer* uc;
        //std::vector< CompressedContainer* >* ccs;

    public:
        Vertex();
        ~Vertex();

};


