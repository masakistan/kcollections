#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "Globals.h"

class Bkmer
{
    private:
        uint8_t* m_bseq;
        size_t m_k, m_bk;
        void resize();
    
    public:
        Bkmer( int k, int bk );
        Bkmer( const Bkmer& other );
        Bkmer* emit_prefix( int len );
        Bkmer* get_prefix( int len );
        Bkmer* get_suffix( int pos );
        ~Bkmer();
        bool operator<( const Bkmer& other );
        bool operator==( const Bkmer& other );
        size_t get_bk() const;
        size_t get_k() const;
        void set_bk( size_t bk ){ m_bk = bk; }
        void set_k( size_t k ) { m_k = k; }
        uint8_t* get_bseq() const;
        void set_bseq( uint8_t* bseq );
};


