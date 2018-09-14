#pragma once

#include "SerializeKmer.h"

class Kdict
{
    private:
        int m_k, m_bk;

    public:
        Kdict( int k );
        insert( char* kmer );

};


