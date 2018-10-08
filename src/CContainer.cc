#include "CContainer.h"
#include "Vertex.h"

void init_cc( CC* cc, int suffix_size )
{
    //cc->prefixes = ( uint8_t* ) calloc( CAPACITY, sizeof( uint8_t ) );
    cc->children = ( Vertex* ) calloc( CAPACITY, sizeof( Vertex ) );
    cc->suffixes = ( uint8_t* ) calloc( CAPACITY, sizeof( uint8_t ) );
    cc->suffix_size = suffix_size;
    
    /*for( int i; i < CAPACITY; i++ )
    {
        init_vertex( &cc->children[ i ] );
    }*/

    cc->size = 0;
    memset( &cc->pref, 0, 32 );
    init_bf( &cc->bf );
    /*std::cout << "init pref" << std::endl;
    for( int i = 0; i < 256; i++ )
    {
        std::cout << i << "\t" << testbit( cc->pref, i ) << std::endl;
    }*/
}

int rank( CC* cc, int clust_num )
{
    int clust_size = cc->size;
    for( int clust_pos = 0; clust_pos < clust_size; clust_pos++ )
    {
        if( cc->children[ clust_pos ].start )
        {
            clust_num--;
        }

        //std::cout << "\trank: " << clust_pos << std::endl;

        if( clust_num == 0 )
        {
            return clust_pos;
        }
    }

    return cc->size;
}

int hamming_weight( CC* cc, int index )
{
    //std::cout << "max hamming index: " << index << std::endl;

    /*for( int i = 0; i < 256; i++ )
    {
        std::cout << i << "\t" << testbit( cc->pref, i ) << std::endl;
    }*/

    int hamming_weight = 0;
    for (
            int i = next_set_bit( cc->pref, 0 );
            i > -1 && i <= index;
            //i = next_set_bit( cc->pref, i + 1 )
            )
    {
        i = next_set_bit( cc->pref, i + 1 );
        //std::cout << "\thamming: " << i << "\t" << i << std::endl;
        hamming_weight++;
    }

    return hamming_weight;
}

void insert_cluster( CC* cc, int idx, uint8_t suffix, bool start )
{
    /*std::cout << "size: " << cc->size << "\tindex: " << idx << std::endl << std::flush;
    std::memmove(
            &cc->prefixes[ idx + 1 ],
            &cc->prefixes[ idx ],
            ( cc->size - idx ) * sizeof( uint8_t )
            );
    cc->prefixes[ idx ] = prefix;
    std::cout << "\tprefix: " << unsigned( cc->prefixes[ idx ] ) << std::endl;*/

    int to_move = cc->size - idx;
    //std::cout << "idx: " << idx << "\tshifting: " << to_move << std::endl;

    Vertex* v = ( Vertex* ) malloc( sizeof( Vertex ) );
    init_vertex( v );
    v->start = start;
    if( to_move > 0 )
    {
        std::memmove(
                &cc->children[ idx + 1 ],
                &cc->children[ idx ],
                ( cc->size - idx ) * sizeof( Vertex )
                );
        std::memmove(
                &cc->suffixes[ idx + 1 ],
                &cc->suffixes[ idx ],
                ( cc->size - idx ) * sizeof( uint8_t )
                );
    }

    std::memmove(
            &cc->children[ idx ],
            v,
            sizeof( Vertex )
            );
    free_vertex( v );
    free( v );

    cc->suffixes[ idx ] = suffix;
    //std::cout << "setting cluster idx: " << unsigned( suffix ) << "\t" << idx << std::endl;

    /*int suffix_idx = cc->suffix_size * idx;
    int suffix_bytes_to_move = ( cc->size - idx ) * cc->suffix_size;
    std::memmove(
            &cc->suffixes[ suffix_idx + cc->suffix_size ],
            &cc->suffixes[ suffix_idx ],
            suffix_bytes_to_move
            );
    std::memmove(
            &cc->suffixes[ suffix_idx ],
            suffix,
            cc->suffix_size
            );*/
    cc->size++;
}

uint8_t* get_suffix( CC* cc, int idx )
{
    //int pos = idx * cc->suffix_size;
    return &cc->suffixes[ idx ];
}

void cc_insert( CC* cc, int k, int depth, uint8_t* sfpx )
{
    //std::cout << "cc insert: " << deserialize_kmer( k, calc_bk( k ), sfpx ) << std::endl;
    if( cc_contains_prefix( cc, sfpx ) )
    {
        return;
    }
    uint8_t prefix = sfpx[ 0 ];
    uint8_t suffix = sfpx[ 0 ];

    //int suffix_size = calc_bk( k ) - depth;
    int pref_index = unsigned( prefix );
    bool was_pref_index_set = testbit( cc->pref, pref_index );
    setbit( cc->pref, pref_index );

    //std::cout << "setting: " << pref_index << "\t" << testbit( cc->pref, pref_index ) << std::endl;

    int clust_num = hamming_weight( cc, pref_index );
    int clust_pos = rank( cc, clust_num );
    //Vertex* nv = ( Vertex* ) malloc( sizeof( Vertex ) );
    //std::cout << "\tclust num: " << clust_num << "\tclust_pos: " << clust_pos << "\tcc size: " << cc->size << std::endl;

    //std::cout << "\t" << nv << "\t" << &nv->uc << "\twas set: " << was_pref_index_set << "\tpref index: " << pref_index << "\tclust_num: " << clust_num << "\tclust_pos: " << clust_pos << std::endl;

    // sfpxPrefix was not present
    if( !was_pref_index_set )
    {
        insert_cluster( cc, clust_pos, suffix, true );
        add_to_bloom_filter( &cc->bf, sfpx, 1 );
        return;
    }

    // sfpxPrefix was already present
    if( was_pref_index_set )
    {
        // if sfpx_suffix starts its cluster...
        //int sfx_cmp = std::memcmp( &suffix, get_suffix( cc, clust_pos ), 1 );
        //if( sfx_cmp < 0 )
        if( unsigned( suffix ) < unsigned( cc->suffixes[ clust_pos ] ) )
        {
            insert_cluster( cc, clust_pos, suffix, false );
            cc->children[ clust_pos + 1 ].start = false;
            add_to_bloom_filter( &cc->bf, sfpx, 1 );
            return;
        }

        // if sfpx_suffix does not start its cluster...
        clust_pos++;

        // if clust_pos is at end of bit-array
        //if( clust_pos >= cc->size )
        if( unsigned( suffix ) >= unsigned( cc->suffixes[ clust_pos ] ) )
        {
            insert_cluster( cc, clust_pos, suffix, false );
            add_to_bloom_filter( &cc->bf, sfpx, 1 );
            return;
        }

        // while sfpx_suffix is greater than previous suffix...
        //sfx_cmp = std::memcmp( &suffix, get_suffix( cc, clust_pos - 1 ), 1 );
        while ( unsigned( suffix ) > unsigned( cc->suffixes[ clust_pos - 1 ] ) )
        {
            // if clust_pos is at end of bit-array OR if next cluster is reached (ie. if sfpx_suffix is greatest in its cluster)...
            if(
                    clust_pos >= cc->size
                    || cc->children[ clust_pos ].start
                    )
            {
                insert_cluster( cc, clust_pos, suffix, false );
                add_to_bloom_filter( &cc->bf, sfpx, 1 );
                return;
            }

            clust_pos++;
            //sfx_cmp = std::memcmp( &suffix, get_suffix( cc, clust_pos - 1 ), 1 );
        }

        // if sfpx_suffix is less than previous suffix...
        std::cout << "less than previous suffix!" << std::endl;
        insert_cluster( cc, clust_pos - 1, suffix, false );
        add_to_bloom_filter( &cc->bf, sfpx, 1 );
        return;
    }
}

bool cc_contains_prefix( CC* cc, uint8_t* sfpx )
{
    uint8_t prefix = sfpx[ 0 ];
    uint8_t suffix = sfpx[ 0 ];

    int pref_index = unsigned( prefix );
    //std::cout << "\t\t\tpref: " << pref_index << "\t" << testbit( cc->pref, pref_index ) << "\t" << cc->size << std::endl;
    /*for( int i = 0 ; i < cc->size; i++ )
    {
        std::cout << i << "\t" << unsigned( cc->suffixes[ i ] ) << std::endl;
    }*/

    if( testbit( cc->pref, pref_index ) )
    {
        int cluster_num = hamming_weight( cc, pref_index );
        int start = rank( cc, cluster_num );
        int pos = start;
        while( pos < cc->size
                && (
                    pos == start
                    || cc->children[ pos ].start == false
                    )
                )
        {
            //std::cout << "checking pos: " << pos << "\t" << unsigned( suffix ) << "\t" << unsigned( cc->suffixes[ pos ] )<< std::endl;
            //int cmp = std::memcmp( &suffix, get_suffix( cc, pos ), 1 );
            if( unsigned( suffix ) == unsigned( cc->suffixes[ pos ] ) )
            {
                //std::cout << "\t\t\tcontains prefix: true" << std::endl;
                return true;
            }
            pos++;
        }
    }

    //std::cout << "\t\t\tcontains prefix: false" << std::endl;
    return false;
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
    free( cc->suffixes );
    for( int i = 0; i < cc->size; i++ )
    {
        free_vertex( &cc->children[ i ] );
    }
    free( cc->children );
    free_bf( &( cc->bf ) );
}

Vertex* get_child_of( CC* cc, uint8_t* sfpx )
{
    int idx = index_of( cc, sfpx );
    //std::cout << "found child index: " << index_of << std::endl;
    if( idx > -1 )
    {
        return &cc->children[ idx ];
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
                    || !cc->children[ pos ].start
                   )
             )
        {
            if( cc->suffixes[ pos ] == suffix )
            {
                return pos;
            }
            pos++;
        }
    }

    return -1;
}


