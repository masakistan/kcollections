#include "CContainer.h"

CContainer::CContainer() : Container()
{
    bf = new BloomFilter( 10, 10 );
}

CContainer::~CContainer()
{
    delete bf;
}

bool CContainer::may_contain( Bkmer* bkmer )
{
    return bf->may_contain( bkmer->bseq, BK );
}


