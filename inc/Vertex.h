#pragma once

#include <vector>
#include <iostream>
//#include <co2/recursive_generator.hpp>
//#include "aco.h"
#include "Container.h"
#include "UContainer.h"
#include "CContainer.h"
#include "SufClustData.h"
#include "Bkmer.h"


class Vertex
{
    private:
        // Containers
        UContainer* m_uc;
        std::vector< CContainer* >* m_ccs;

    public:
        Vertex();
        ~Vertex();
        UContainer* get_uc() const { return m_uc; }
        std::vector< CContainer* >* get_ccs() const { return m_ccs; }

        void burst_uc( Bkmer* bkmer );
        size_t size();
        
        void insert( Bkmer* bkmer );
        static void insert( Vertex* v, Bkmer* bkmer );
        
        bool contains( Bkmer* bkmer );
        static bool contains( const Vertex* v, Bkmer* bkmer );
        
        void remove( Bkmer* bkmer );
        static void remove( Vertex* v, Bkmer* bkmer );
};


