#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kdict
{
    private:
        Kcontainer* kc;
        int m_k;
    public:
        Kdict( const int k );
        ~Kdict();
        void add( char* kmer, py::handle* obj );
        bool contains( char* kmer );
        void clear();
        uint64_t size();
        void remove( char* kmer );
        py::handle* get( char* kmer );
        int get_k() { return m_k; }
        Kcontainer* get_kc() { return kc; }

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
        int get_vs_size( Vertex* v ){ return v->vs_size; }
        Vertex* get_child_vertex( Vertex* v, int idx )
        {
            return &v->vs[idx];
        }
        char* get_child_suffix( Vertex* v, int idx )
        {
            return kcontainer_get_child_suffix(v, idx);
        }
        Vertex* get_root() { return &kc->v; }
};
