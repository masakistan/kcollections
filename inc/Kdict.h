#pragma once

#include <stdlib.h>
#include <vector>
#include "Globals.h"
#include "Bkmer.h"
#include "SerializeKmer.h"
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

};


