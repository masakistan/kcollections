#pragma once

#include <memory>
#include <vector>
#include <bitset>
#include <cmath>
#include "Container.h"
#include "BloomFilter.h"
#include "Bkmer.h"
#include "SufClustData.h"


class __attribute__ ((__packed__)) CContainer : public Container
{
    private:
        BloomFilter* bf;
        std::vector< bool >* m_pref;
        std::vector< SufClustData* >* m_suf_clust_data;
        std::bitset< 256 > cluster_starts;
        const static int PREF_SIZE;

        unsigned int get_index_in_pref( Bkmer* bkmer );
        int hamming_weight( int index );
        int rank( int clust_num );
        int clust_num_from_rank( int clust_pos );
        int pref_index_from_hamming_weight( int clust_num );
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
        char* prefix_from_clust( int clust_pos );

        std::vector< SufClustData* >* get_suf_clust_data() { return m_suf_clust_data; }
        BloomFilter* get_bf();
        std::vector< Vertex* >* get_child_vertices();
        static int get_pref_size() { return PREF_SIZE; }
        size_t size();

};


