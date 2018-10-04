#pragma once

#include <vector>
#include <iostream>
#include <memory>
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
        //std::vector< CContainer* >* m_ccs;
        std::array< CContainer*, 26 >* m_ccs;
        uint8_t m_ccs_pos = 0;

    public:
        Vertex();
        ~Vertex();
        UContainer* get_uc() const { return m_uc; }
        std::array< CContainer*, 26 >* get_ccs() { return m_ccs; }
        //std::vector< CContainer* >* get_ccs() const { return m_ccs; }

        void burst_uc( Bkmer* bkmer );
        size_t size();
        
        void insert( Bkmer* bkmer );
        static void insert( Vertex* v, Bkmer* bkmer );
        
        bool contains( Bkmer* bkmer );
        static bool contains( Vertex* v, Bkmer* bkmer );
        
        void remove( Bkmer* bkmer );
        static void remove( Vertex* v, Bkmer* bkmer );
};


