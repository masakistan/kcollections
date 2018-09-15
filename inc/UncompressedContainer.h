#pragma once

#include "Container.h"

class UncompressedContainer : public Container
{
    public:
        UncompressedContainer();
        void insert();
        //bool contains( uint8_t* bkmer );
};


