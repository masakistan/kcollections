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
        std::vector< std::unique_ptr< SufClustData > >* m_suf_clust_data; 
        const static int PREF_SIZE;

        unsigned int get_index_in_pref( Bkmer* bkmer );
        int hamming_weight( int index );
        int rank( int clust_num );
        int index_of( Bkmer* bkmer );
        void add_to_bloom_filter( Bkmer* bkmer );
        int next_set_bit( int pos );

    public:
        CContainer();
        ~CContainer();
        bool may_contain( Bkmer* sfpx );
        Vertex* get_child_of( Bkmer* spfx );
        bool contains_prefix( Bkmer* sfpx );
        void insert( Bkmer* sfpx );
        bool is_full();
        char* index_to_pref( uint8_t index );

        //int size() { return m_suf_clust_data->size(); }
        BloomFilter* get_bf();
        //std::vector< std::unique_ptr< SufClustData > >* get_suf_clust_data();
        //std::vector< Bkmer* >* get_suf();
        std::vector< bool >* get_clust();
        std::vector< Vertex* >* get_child_vertices();
        static int get_pref_size() { return PREF_SIZE; }
        //std::unique_ptr< SufClustData > get_suf_clust_data_item( Bkmer* sfpx );
        size_t size();

};


