#pragma once

#include <vector>
#include <iostream>
#include "Container.h"
#include "UContainer.h"
#include "CContainer.h"
#include "SufClustData.h"
#include "Bkmer.h"

class Vertex
{
    private:
        int depth;

        // Containers
        std::unique_ptr< UContainer > m_uc;
        std::vector< std::unique_ptr< CContainer > >* m_ccs;

        UContainer* get_uc() const { return m_uc.get(); }
        std::vector< std::unique_ptr< CContainer > >* get_ccs() const { return m_ccs; }

    public:
        Vertex();
        ~Vertex();

        void insert( Vertex* v, Bkmer* bkmer );
        void insert( Bkmer* bkmer );
        bool contains( Bkmer* bkmer );
        bool contains( Vertex* v, Bkmer* bkmer ) const;
        void burst_uc( Bkmer* bkmer );
        size_t size();
        void remove( Vertex* v, Bkmer* bkmer );
        void remove( Bkmer* bkmer );

};


