#pragma once

#include <stdlib.h>
#include <vector>
#include <co2/generator.hpp>
#include <co2/recursive_generator.hpp>
#include "Globals.h"
#include "Bkmer.h"
#include "Vertex.h"
#include "Helper.h"


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
        //auto get_bkmers(Vertex*v) -> co2::generator< Bkmer >;
        auto get_bkmers() { return root->get_bkmers(); }
        
};


