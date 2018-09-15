#pragma once

#include "Bkmer.h"

class Container
{
    private:
        static int s_capacity;

    public:
        Container(){}
        void set_capacity( int capacity ){ s_capacity = capacity; }
        int get_capacity(){ return s_capacity; }
        virtual void insert( Bkmer bkmer ) = 0;

};


