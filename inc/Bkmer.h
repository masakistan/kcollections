#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "Globals.h"

class Bkmer
{
    private:
        uint8_t* m_bseq;
        unsigned int m_k, m_bk;
    
    public:
        Bkmer( int k, int bk );
        Bkmer( const Bkmer& other );
        uint8_t* emit_prefix( int len );
        uint8_t* get_prefix( int len );
        ~Bkmer();
        bool operator<( Bkmer& other ) const;
        int get_bk() const;
        int get_k() const;
        uint8_t* get_bseq() const;
        void set_bseq( uint8_t* bseq );
};


