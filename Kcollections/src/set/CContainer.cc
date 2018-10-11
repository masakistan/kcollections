#include "set/CContainer.h"
#include "set/Vertex.h"

void init_cc( CC* cc, int suffix_size )
{
    cc->child_suffixes = NULL;
    cc->suffix_size = suffix_size;
    
    cc->size = 0;
    memset( &cc->pref, 0, 32 );
    init_bf( &cc->bf );
}

int rank( CC* cc, int clust_num )
{
    int clust_size = cc->size;
    for( int clust_pos = 0; clust_pos < clust_size; clust_pos++ )
    {
        if( cc->child_suffixes[ clust_pos ].v.start )
        {
            clust_num--;
        }

        if( clust_num == 0 )
        {
            return clust_pos;
        }
    }

    return cc->size;
}

int hamming_weight( CC* cc, int index )
{
    int hamming_weight = 0;
    for (
            int i = next_set_bit( cc->pref, 0, 256 );
            i > -1 && i <= index;
            i = next_set_bit( cc->pref, i + 1, 256 )
            )
    {
        hamming_weight++;
    }

    return hamming_weight;
}

Vertex* insert_cluster( CC* cc, int idx, uint8_t suffix, bool start )
{
    int to_move = cc->size - idx;
    
    if( cc->child_suffixes == NULL )
    {
        cc->child_suffixes = ( CS* ) calloc( 1, sizeof( CS ) );
    }
    else
    {
        cc->child_suffixes = ( CS* ) realloc(
                cc->child_suffixes,
                ( cc->size + 1 ) * sizeof( CS )
                );
    }


    Vertex* v = ( Vertex* ) malloc( sizeof( Vertex ) );
    init_vertex( v );
    v->start = start;
    if( to_move > 0 )
    {
        std::memmove(
                &cc->child_suffixes[ idx + 1 ],
                &cc->child_suffixes[ idx ],
                to_move * sizeof( CS )
                );
    }

    std::memcpy(
            &cc->child_suffixes[ idx ].v,
            v,
            sizeof( Vertex )
            );
    free_vertex( v );
    free( v );

    cc->child_suffixes[ idx ].suffix = suffix;
    cc->size++;
    return &cc->child_suffixes[ idx ].v;
}

Vertex* cc_insert( CC* cc, int k, int depth, uint8_t* sfpx )
{
    int idx = cc_contains_prefix( cc, sfpx );
    if( idx > -1 )
    {
        return get_child_of( cc, sfpx, idx );
    }
    uint8_t prefix = sfpx[ 0 ];
    uint8_t suffix = sfpx[ 0 ];

    //int suffix_size = calc_bk( k ) - depth;
    int pref_index = unsigned( prefix );
    bool was_pref_index_set = testbit( cc->pref, pref_index );
    setbit( cc->pref, pref_index );


    int clust_num = hamming_weight( cc, pref_index );
    int clust_pos = rank( cc, clust_num );
    Vertex* child;

    // sfpxPrefix was not present
    if( !was_pref_index_set )
    {
        child = insert_cluster( cc, clust_pos, suffix, true );
        add_to_bloom_filter( &cc->bf, sfpx, 1 );
        return child;
    }

    // sfpxPrefix was already present
    if( was_pref_index_set )
    {
        // if sfpx_suffix starts its cluster...
        if( unsigned( suffix ) < unsigned( cc->child_suffixes[ clust_pos ].suffix ) )
        {
            child = insert_cluster( cc, clust_pos, suffix, false );
            cc->child_suffixes[ clust_pos + 1 ].v.start = false;
            add_to_bloom_filter( &cc->bf, sfpx, 1 );
            return child;
        }

        // if sfpx_suffix does not start its cluster...
        clust_pos++;

        // if clust_pos is at end of bit-array
        //if( clust_pos >= cc->size )
        if( unsigned( suffix ) >= unsigned( cc->child_suffixes[ clust_pos ].suffix ) )
        {
            child = insert_cluster( cc, clust_pos, suffix, false );
            add_to_bloom_filter( &cc->bf, sfpx, 1 );
            return child;
        }

        // while sfpx_suffix is greater than previous suffix...
        while ( unsigned( suffix ) > unsigned( cc->child_suffixes[ clust_pos - 1 ].suffix ) )
        {
            // if clust_pos is at end of bit-array OR if next cluster is reached (ie. if sfpx_suffix is greatest in its cluster)...
            if(
                    clust_pos >= cc->size
                    || cc->child_suffixes[ clust_pos ].v.start
                    )
            {
                child = insert_cluster( cc, clust_pos, suffix, false );
                add_to_bloom_filter( &cc->bf, sfpx, 1 );
                return child;
            }

            clust_pos++;
        }

        // if sfpx_suffix is less than previous suffix...
        std::cout << "less than previous suffix!" << std::endl;
        child = insert_cluster( cc, clust_pos - 1, suffix, false );
        add_to_bloom_filter( &cc->bf, sfpx, 1 );
        return child;
    }
}

int cc_contains_prefix( CC* cc, uint8_t* sfpx )
{
    uint8_t prefix = sfpx[ 0 ];
    uint8_t suffix = sfpx[ 0 ];
    int pref_index = unsigned( prefix );

    if( testbit( cc->pref, pref_index ) )
    {
        int cluster_num = hamming_weight( cc, pref_index );
        int start = rank( cc, cluster_num );
        int pos = start;
        while( pos < cc->size
                && (
                    pos == start
                    || cc->child_suffixes[ pos ].v.start == false
                    )
                )
        {
            if( unsigned( suffix ) == unsigned( cc->child_suffixes[ pos ].suffix ) )
            {
                return pos;
            }
            pos++;
        }
    }

    return -1;
}

bool cc_may_contain( CC* cc, uint8_t* bseq )
{
    bool res = bf_may_contain( &cc->bf, bseq, 1);
    //std::cout << "\t\t\tmay contain: " << res << std::endl;
    return res;
}

void cc_shrink( CC* cc )
{
}

void free_cc( CC* cc )
{
    //free( cc->suffixes );
    for( int i = 0; i < cc->size; i++ )
    {
        free_vertex( &cc->child_suffixes[ i ].v );
    }
    free( cc->child_suffixes );
    //free_bf( &( cc->bf ) );
}

Vertex* get_child_of( CC* cc, uint8_t* sfpx, int idx )
{
    if( idx > -1 )
    {
        return &cc->child_suffixes[ idx ].v;
    }
    return NULL;
}

int index_of( CC* cc, uint8_t* sfpx )
{
    uint8_t prefix = sfpx[ 0 ];
    uint8_t suffix = sfpx[ 0 ];
    int pref_index = unsigned( prefix );

    if( testbit( cc->pref, pref_index ) )
    {
        int cluster_num = hamming_weight( cc, pref_index );
        int start = rank( cc, cluster_num );
        int pos = start;
        while( pos < cc->size
                && (
                    pos == start
                    || !cc->child_suffixes[ pos ].v.start
                   )
             )
        {
            if( unsigned( cc->child_suffixes[ pos ].suffix ) == unsigned( suffix ) )
            {
                return pos;
            }
            pos++;
        }
    }

    return -1;
}


