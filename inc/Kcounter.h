#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kcounter
{
    private:
        Kcontainer* kc;
        int m_k;
    public:
        Kcounter( const int k );
        ~Kcounter();
        void insert( const char* kmer, int count );
        bool contains( const char* kmer );
        void clear();
        uint64_t size();
        void remove( const char* kmer );
        void add_seq( const char* seq, uint32_t length );
        int get( const char* kmer );
        int get_k() { return m_k; }
        Kcontainer* get_kc() { return kc; }
        void parallel_add_init(int threads) {
            parallel_kcontainer_add_init(kc, threads);
        };
        void parallel_add(const char* kmer) {
            parallel_kcontainer_add(kc, kmer);
        }
        void parallel_add_join() {
            parallel_kcontainer_add_join(kc);
        }

        void parallel_add_seq(const char* seq, uint32_t length);

        char* get_uc_kmer( Vertex* v, int k, int idx )
        {
            UC* uc = &v->uc;
            int bk = calc_bk( k );
            int suffix_idx = bk * idx;
            //std::cout << bk << "\t" << suffix_idx << "\t" << idx << std::endl;
            //for(int i = 0; i < bk; i++) {
            //    std::cout << (unsigned) uc->suffixes[0] << std::endl;
            //}
            return deserialize_kmer( k, bk, &uc->suffixes[ suffix_idx ] );
        }

        int get_uc_size( Vertex* v )
        {
            return v->uc.size;
        }

        Vertex* get_root() { return &kc->v; }
        int get_vs_size( Vertex* v ){ return v->vs_size; }
        Vertex* get_child_vertex( Vertex* v, int idx )
        {
            return &v->vs[idx];
        }
        char* get_child_suffix( Vertex* v, int idx )
        {
            return kcontainer_get_child_suffix(v, idx);
        }
};
