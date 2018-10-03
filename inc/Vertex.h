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


        void insert( Vertex* v, Bkmer* bkmer );
        void insert( Bkmer* bkmer );
        bool contains( Bkmer* bkmer );
        bool contains( Vertex* v, Bkmer* bkmer ) const;
        void burst_uc( Bkmer* bkmer );
        size_t size();
        void remove( Vertex* v, Bkmer* bkmer );
        void remove( Bkmer* bkmer );
        /*auto get_bkme=s() { return get_bkmers( this ); }
        static auto get_bkmers( Vertex* v ) -> co2::recursive_generator< char* >;*/
};


