#pragma once

#include "Container.h"
#include "BloomFilter.h"
#include "Bkmer.h"
#include "Globals.h"

class CContainer : public Container
{
    private:
        BloomFilter* bf;

    public:
        CContainer();
        ~CContainer();
        bool may_contain( Bkmer* bkmer );

};


