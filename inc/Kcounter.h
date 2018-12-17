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
        void insert( char* kmer, int count );
        bool contains( char* kmer );
        void clear();
        uint64_t size();
        void remove( char* kmer );
        void add_seq( char* seq, uint32_t length );
        int get( char* kmer );
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

        void parallel_add_seq(char* seq, uint32_t length);

        char* get_uc_kmer( Vertex* v, int k, int idx )
        {
            UC* uc = &v->uc;
            int bk = calc_bk( k );
            int suffix_idx = bk * idx;
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
            uint256_t verts = v->pref_pres;
            uint8_t i = 0;
            while(i < idx) {
                if(verts & 0x1) {
                    idx++;
                }
                verts >>= 1;
            }
            return deserialize_kmer(4, 1, &i);
        }
};
