#pragma once

#include <set>

#include "Container.h"
#include "Bkmer.h"

class UContainer : public Container
{
    private:
        std::set< Bkmer >* bkmers;

    public:
        UContainer();
        void insert( Bkmer bkmer );
        bool contains_kmer( Bkmer* bkmer );
};


