#include "Vertex.h"

Vertex::Vertex()
{
    //ccs = new std::vector< CompressedContainer*>();
    uc = new UContainer();
}

Vertex::~Vertex()
{
    delete uc;
}

UContainer* Vertex::get_uc()
{
    return uc;
}


