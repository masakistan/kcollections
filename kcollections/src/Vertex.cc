#include "Vertex.h"


void init_vertex( Vertex* v )
{
    v->start = false;
    v->cc = NULL;
    v->cc_size = 0;
    init_uc( &( v->uc ) );
}

CC* get_cc( Vertex* v, int idx )
{
    return &v->cc[ idx ];
}

#if defined KDICT || defined KCOUNTER
#if KDICT
py::handle* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOUNTER
int vertex_get_counter( Vertex* v, uint8_t* bseq, int k, int depth )
#endif
{
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
#if KDICT
        return &v->uc.objs[ uc_idx ];
#elif KCOUNTER
        return v->uc.count;
#endif
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
#if KDICT
                    return vertex_get( child, &bseq[ 1 ], k - 4, depth + 1 );
#elif KCOUNTER
                    return vertex_get_counter( child, &bseq[ 1 ], k - 4, depth + 1 );
#endif
                }
            }
        }
    }

#if KCOUNTER
    // add a vertex to be 0 if key is not found
    int default_count = 0;
    vertex_insert(v, bseq, k, depth, default_count);
    return default_count;
#endif

    throw pybind11::key_error( "Key not in dictionary!" );
}
#endif

void vertex_remove( Vertex* v, uint8_t* bseq, int k, int depth )
{
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
        uc_remove( &( v->uc ), calc_bk( k ), uc_idx );
        return;
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
                    vertex_remove( child, &bseq[ 1 ], k - 4, depth + 1 );
                }
            }
        }
    }

    throw pybind11::key_error( "Key not found!" );
}

bool vertex_contains( Vertex* v, uint8_t* bseq, int k, int depth )
{
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    if( sres.first )
    {
        return true;
    }


    if( v->cc != NULL )
    {
        int size = sizeof( *v->cc ) / sizeof( CC );
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
#if KDICT
    py::handle* objs = v->uc.objs;
#elif KCOUNTER
    int count = v->uc.count;
#endif
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
#if KDICT
            vertex_insert( child, suffix, k - 4, depth + 1, &objs[ i ] );
#elif KSET
            vertex_insert( child, suffix, k - 4, depth + 1 );
#elif KCOUNTER
            vertex_insert( child, suffix, k - 4, depth + 1, count );
#endif
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

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::handle* obj )
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOUNTER
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, int count )
#endif
{
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
        // set the object here
#if KDICT
        v->uc.objs[ uc_idx ].dec_ref();
        std::memcpy(
                &v->uc.objs[ uc_idx ],
                obj,
                sizeof( py::handle )
                );
        v->uc.objs[ uc_idx ].inc_ref();
#elif KCOUNTER
        v->uc.count = count;
#endif
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
#if KDICT
            vertex_insert( child, suffix, k - 4, depth + 1, obj );
#elif KSET
            vertex_insert( child, suffix, k - 4, depth + 1 );
#elif KCOUNTER
            vertex_insert( child, suffix, k - 4, depth + 1, count );
#endif
            return;
        }
    }

    if( v->uc.size < CAPACITY - 1 )
    {
#if KDICT
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#elif KCOUNTER
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, 0 );
#endif
    }
    else
    {
#if KDICT
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#elif KCOUNTER
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, 0 );
#endif
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

uint64_t vertex_size( Vertex* v )
{
    uint64_t c = v->uc.size;
    for( int i = 0; i < v->cc_size; i++ )
    {
        for( int j = 0; j < v->cc[ i ].size; j++ )
        {
            c += vertex_size( &v->cc[ i ].child_suffixes[ j ].v );
        }
    }
    return c;
}

