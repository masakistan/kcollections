#include "Vertex.h"


void init_vertex( Vertex* v )
{
    v->start = false;
    v->cc = NULL;
    v->cc_size = 0;
    init_uc( &( v->uc ) );
}

bool vertex_contains( Vertex* v, uint8_t* bseq, int k, int depth )
{
    int uc_idx = uc_contains( &( v->uc ), k, depth, bseq );
    if( uc_idx < 0 )
    {
        return true;
    }


    if( v->cc != NULL )
    {
        int size = sizeof( *v->cc ) / sizeof( CC );
        bool res = false;
        CC* cc;
        for( int i = 0; i < v->cc_size; i++ )
        {
            cc = &v->cc[ i ];
            if( cc_may_contain( cc, bseq ) )
            {
                int idx = cc_contains_prefix( cc, bseq );
                if ( idx > -1 )
                {
                    Vertex* child = get_child_of( cc, bseq, idx );
                    return vertex_contains( child, &bseq[ 1 ], k - 4, depth + 1 );
                }
            }
        }
    }

    return false;
}

void burst_uc( Vertex* v, int k, int depth )
{
    CC* cc = ( CC* ) malloc( sizeof( CC ) );
    int suffix_size = calc_bk( k );
    init_cc( cc, suffix_size );

    uint8_t* suffixes = v->uc.suffixes;
    int idx;
    for( int i = 0; i < CAPACITY; i++ )
    {
        idx = i * suffix_size;

        uint8_t* bseq = &suffixes[ idx ];
        uint8_t prefix = bseq[ 0 ];
        uint8_t* suffix = &bseq[ 1 ];

        // check if full TODO finish
        if( !false )
        {
            Vertex* child = cc_insert( cc, k, depth, bseq );
            //Vertex* child = get_child_of( cc, bseq, index_of( cc, bseq ) );
            vertex_insert( child, suffix, k - 4, depth + 1 );
        }
    }

    if( v->cc == NULL )
    {
        v->cc = cc;
        v->cc_size++;
    }
    else
    {
        v->cc = ( CC* ) realloc( v->cc, ( v->cc_size + 1 ) * sizeof( CC ) );
        std::memcpy( &v->cc[ v->cc_size ], cc, sizeof( CC ) );
        v->cc_size++;
        free( cc );
    }

    free_uc( &( v->uc ) );
    init_uc( &( v->uc ) );
    //v->uc.size = 0;
}

void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
{
    int uc_idx = uc_contains( &( v->uc ), k, depth, bseq );
    if( uc_idx < 0 )
    {
        return;
    }
    
    for( int i = 0; i < v->cc_size; i++ )
    {
        CC* cc = &v->cc[ i ];
        if( cc_may_contain( cc, bseq ) )
        {
            uint8_t* suffix = &bseq[ 1 ];
            int idx = cc_contains_prefix( cc, bseq );
            Vertex* child;
            if( idx > -1 )
            {
                child = get_child_of( cc, bseq, idx );
            }
            else
            {
                child = cc_insert( cc, k, depth, bseq );
            }
            vertex_insert( child, suffix, k - 4, depth + 1 );
            return;
        }
    }

    if( v->uc.size < CAPACITY - 1 )
    {
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
    }
    else
    {
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
        burst_uc( v, k, depth );
    }
}

void free_vertex( Vertex* v )
{
    free_uc( &( v->uc ) );
    if( v->cc != NULL )
    {
        for( int i = 0; i < v->cc_size; i++ )
        {
            free_cc( &v->cc[ i ] );
        }
        free( v->cc );
    }
}

