#pragma once

#include <stdlib.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <boost/coroutine2/all.hpp>
#include "Globals.h"
#include "Bkmer.h"
#include "Vertex.h"
#include "Helper.h"

namespace py = pybind11;
//typedef boost::coroutines2::coroutine<char*> coro_t;

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
        
        /*static void get_kmers( coro_t::push_type& yield, Kdict* kdict )
        {
            char* seq = ( char* ) malloc( sizeof( char ) * ( kdict->m_k + 1 ) );
            seq[ kdict->m_k ] = '\0';
            yield_kmers( yield, kdict->get_root(), seq, 0 );
            free( seq );
        }

        static void yield_kmers( coro_t::push_type& yield, Vertex* v, char* seq, int pos )
        {
            for( Bkmer uc_bkmer : *v->get_uc()->get_bkmers() )
            {
                //std::cout << "uc" << std::endl << std::flush;
                char* tseq = uc_bkmer.get_seq();
                strcpy( seq + pos, tseq );
                free( tseq );
                yield( seq );
            }

            for( CContainer* cc : *v->get_ccs() )
            {
                //std::unique_ptr< SufClustData >& sfc : *cc->get_suf_clust_data() )
                //std::cout << "cc" << std::endl << std::flush;
                std::vector< SufClustData* >* sfcs = cc->get_suf_clust_data();
                for( int i = 0; i < sfcs->size(); i++ )
                {
                    //std::cout << "\tclust " << i << std::endl << std::flush;
                    SufClustData* sfc = sfcs->at( i );
                    char* prefix = cc->prefix_from_clust( i );
                    strcpy( &seq[ pos ], prefix );
                    free( prefix );
                    yield_kmers( yield, sfc->get_child_vertex(), seq, pos + 4);
                    //std::cout << "\tend clust " << i << std::endl;
                }
                //std::cout << "end cc" << std::endl;
            }
        }*/
};


