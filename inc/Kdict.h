#pragma once

#include "SerializeKmer.h"
#include "Helper.h"

class Kdict
{
    private:

    public:
        const int m_k, m_bk;

        Kdict( int k, int bk );
        void insert( char* kmer );
};


