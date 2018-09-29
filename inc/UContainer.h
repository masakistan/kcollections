#pragma once

#include <set>

#include "Container.h"
#include "Bkmer.h"

class UContainer : public Container
{
    private:
        std::set< Bkmer >* m_bkmers;

    public:
        UContainer();
        ~UContainer();
        int size() { return m_bkmers->size(); }
        bool contains( Bkmer* bkmer );
        void insert( Bkmer* bkmer );
        bool is_full();
        std::set< Bkmer >* get_bkmers();
        void remove( Bkmer* bkmer );
};


