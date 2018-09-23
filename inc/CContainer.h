#pragma once

#include <vector>
#include <cmath>
#include "Container.h"
#include "BloomFilter.h"
#include "Bkmer.h"
#include "Globals.h"
#include "SufClustData.h"


class CContainer : public Container
{
    private:
        BloomFilter* bf;
        std::vector< bool >* m_pref;
        std::vector< SufClustData* >* m_suf_clust_data; 
        const static int PREF_SIZE = (int) Math.pow(2, (double) Container.numBitsNeededForAlphabet() * Container.getPrfxPrefixLength());

        int get_index_in_pref( Bkmer* bkmer );
        int hamming_weight( int index );
        int rank( int clust_num );
        int get_index_in_pref( Bkmer* bkmer );
        int index_of( Bkmer* bkmer );
        void add_to_bloom_filter( Bkmer* bkmer );

    public:
        CContainer();
        ~CContainer();
        bool may_contain( Bkmer* sfpx );
        Vertex* get_child_of( Bkmer* spfx );
        bool contains_prefix( Bkmer* sfpx );
        void insert( Bkmer* sfpx );
        bool is_full();

        int size() { return m_suf_clust_data->size(); }
        BloomFilter* get_bf();
        std::vector< bool >* get_suf_clust_data();
        std::vector< Bkmer* >* get_suf();
        std::vector< bool >* get_clust();
        std::vector< Vector* >* get_child_vertices();
        static int get_pref_size() { return PREF_SIZE; }

};


