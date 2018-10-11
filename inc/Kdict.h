#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kdict
{
    private:
        Kcontainer* kc;
        int k;
    public:
        Kdict( const int k );
        ~Kdict();
        void insert( char* kmer, py::object* obj );
        bool contains( char* kmer );
        void clear();
        uint64_t size();
        void remove( char* kmer );
        py::object* get( char* kmer );
};


