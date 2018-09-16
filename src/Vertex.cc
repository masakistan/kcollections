#include "Vertex.h"

Vertex::Vertex()
{
    ccs = new std::vector< CContainer*>();
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

std::vector< CContainer* >* Vertex::get_cc()
{
    return ccs;
}


