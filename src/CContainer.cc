#include "CContainer.h"

CContainer::CContainer() : Container()
{
    bf = new BloomFilter( 10, 10 );
}

CContainer::~CContainer()
{
    delete bf;
}

bool CContainer::may_contain()
{
    return bf->may_contain();
}


