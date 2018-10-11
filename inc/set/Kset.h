#pragma once

#include <pybind11/pybind11.h>
#include "set/Kcontainer.h"

namespace py = pybind11;

class Kset
{
    private:
        Kcontainer* kc;
        int k;
    public:
        Kset( const int k );
        ~Kset();
        void insert( char* kmer );
        bool contains( char* kmer );
        void clear();
        uint64_t size();
        void remove( char* kmer );
        //void remove( char* kmer );
};


