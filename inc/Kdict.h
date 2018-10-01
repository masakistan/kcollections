#pragma once

#include <stdlib.h>
#include <vector>
//#include <co2/generator.hpp>
//#include <co2/recursive_generator.hpp>
#include <boost/coroutine2/all.hpp>
#include "Globals.h"
#include "Bkmer.h"
#include "Vertex.h"
#include "Helper.h"

typedef boost::coroutines2::coroutine<char*> coro_t;

class Kdict
{
    private:
        Vertex* root;

    public:
        const int m_k, m_bk;

        Kdict( int k, int bk );
        ~Kdict();
        void insert( char* kmer );
        bool contains( char* kmer );
        size_t size();
        void remove( char* kmer );
        void clear();
        Vertex* get_root() { return root; }
        //void get_kmers( coro_t::push_type& yield ) { get_kmers( yield, root ); }
        
        static void get_kmers(coro_t::push_type& yield, Vertex* v)
        {
            for( Bkmer uc_bkmer : *v->get_uc()->get_bkmers() )
            {
                yield( uc_bkmer.get_seq() );
            }

            for( std::unique_ptr< CContainer>& cc : *v->get_ccs() )
            {
                for( std::unique_ptr< SufClustData >& sfc : *cc->get_suf_clust_data() )
                {
                    get_kmers( yield, sfc->get_child_vertex() );
                }
            }
        }
};


