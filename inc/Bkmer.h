#pragma once

#include <stdlib.h>
#include <sstream>
#include <stdint.h>
#include "Globals.h"

static const uint8_t MASK_INSERT[ 3 ][ 4 ] = {
        { 1, 4, 16, 64 },
        { 2, 8, 32, 128 },
        { 3, 12, 48, 192 }
    };

static const char COMP_TO_ASCII[4] = {'A', 'C', 'G', 'T'};

class Bkmer
{
    private:
        uint8_t* m_bseq;
        size_t m_k, m_bk;
        void resize();
        void serialize_kmer( char* kmer );
    
    public:
        Bkmer( int k, int bk, char* kmer );
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
        char char_at( int pos );
        char* deserialize_seq();
};


