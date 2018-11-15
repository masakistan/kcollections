#pragma once

#include <pybind11/pybind11.h>
#include "Kdict.h"

namespace py = pybind11;

class Kcounter: public Kdict
{
    public:
        Kcounter( const int k );
        ~Kcounter();
        void insert( char* kmer, int count );
        int get( char* kmer );
};
