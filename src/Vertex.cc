#include "Vertex.h"

Vertex::Vertex()
{
    //ccs = new std::vector< CompressedContainer*>();
    uc = new UncompressedContainer();
}

Vertex::~Vertex()
{
    delete uc;
}


