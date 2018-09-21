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
        std::vector< SufClustData* >* m_suf_clust_data; 

    public:
        CContainer();
        ~CContainer();
        bool may_contain( Bkmer* sfpx );
        bool contains_prefix( Bkmer* sfpx );
        void insert( Bkmer* sfpx );
        Vertex* get_child_of( Bkmer* spfx );
        bool is_full();

};


