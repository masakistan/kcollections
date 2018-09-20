#pragma once

#include "Bkmer.h"

class Container
{
    protected:
        static int s_capacity;

    private:

    public:
        Container(){}
        void set_capacity( int capacity ){ s_capacity = capacity; }
        int get_capacity(){ return s_capacity; }
        //virtual void insert( uint8_t* sfpx ) = 0;

};


