#pragma once

#include <assert.h>
#include "Bkmer.h"

class Container
{
    protected:
        static uint16_t s_capacity;
        static uint8_t s_prefix_length;
        static uint8_t s_prfx_prefix_length;
        static uint8_t s_prfx_suffix_length;

    private:

    public:
        Container(){}
        static void set_capacity( int capacity ){ s_capacity = capacity; }
        static int get_capacity(){ return s_capacity; }
        static int get_prefix_length(){ return s_prefix_length; }
        static int get_prfx_prefix_length(){ return s_prfx_prefix_length; }
        static int get_prfx_suffix_length(){ return s_prfx_suffix_length; }
        virtual void insert( Bkmer* sfpx ) = 0;

};


