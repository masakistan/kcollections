#pragma once

#include <vector>
#include "Container.h"
#include "BloomFilter.h"
#include "Bkmer.h"
#include "Globals.h"
#include "SufClustData.h"


class CContainer : public Container
{
    private:
        BloomFilter* bf;
        bool contains_prefix( uint8_t* sfpx );
        //std::vector

    public:
        CContainer();
        ~CContainer();
        bool may_contain( uint8_t* sfpx );
        void insert( uint8_t* sfpx );
        Vertex* get_child_of( uint8_t* spfx );

};


