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
    //char* seq = deserialize_kmer( k, calc_bk( k ), bseq );
    //std::cout << "\nfor: " << seq << "\n\tchecking in uc: " << v->uc.size << "\t" << depth << std::endl;
    //int ak = k - ( depth * 4 );
    //std::cout << "\n" << depth << " searching for: " << deserialize_kmer( ak, calc_bk( ak ), bseq ) << std:: endl;
    //std::cout << "ucs present" << std::endl;
    //print( &( v->uc ), k, depth );
    if( uc_contains( &( v->uc ), k, depth, bseq ) )
    {
        return true;
    }


    if( v->cc != NULL )
    {
        int size = sizeof( *v->cc ) / sizeof( CC );
        bool res = false;
        //std::cout << "\tneed check " << size << " cc nodes" << std::endl;
        CC* cc;
        for( int i = 0; i < v->cc_size; i++ )
        {
            //std::cout << "\t\tchecking: " << i << std::endl;
            cc = &v->cc[ i ];
            if( cc_may_contain( cc, bseq ) && cc_contains_prefix( cc, bseq ) )
            {
                return vertex_contains( get_child_of( cc, bseq ), &bseq[ 1 ], k - 4, depth + 1 );
            }
        }
    }

    //std::cout << depth << " not found in vertex! " << depth << std::endl;
    return false;
}

void burst_uc( Vertex* v, int k, int depth )
{
    //std::cout << "\nbursting at depth: " << depth << std::endl;
    CC* cc = ( CC* ) malloc( sizeof( CC ) );
    int suffix_size = calc_bk( k );
    init_cc( cc, suffix_size );
    //std::cout << "bf size: " << cc->bf.m_nHashes << std::endl;

    uint8_t* suffixes = v->uc.suffixes;
    int idx;
    for( int i = 0; i < CAPACITY; i++ )
    {
        idx = i * suffix_size;
        char* seq = deserialize_kmer( k, calc_bk( k ), &suffixes[ idx ] );
        //std::cout << "\tburst kmer: " << seq << std::endl;
        free( seq );

        uint8_t* bseq = &suffixes[ idx ];
        uint8_t prefix = bseq[ 0 ];
        uint8_t* suffix = &bseq[ 1 ];

        // check if full TODO finish
        if( !false )
        {
            cc_insert( cc, k, depth, bseq );
            Vertex* child = get_child_of( cc, bseq );
            //char* prefix = deserialize_kmer( 4, 1, bseq );
            //std::cout << prefix << "\t" << child << std::endl;
            vertex_insert( child, suffix, k - 4, depth + 1 );
        }
    }

    if( v->cc == NULL )
    {
        //std::cout << "\tstarting cc list" << std::endl;
        v->cc = cc;
        v->cc_size++;
    }
    else
    {
        //std::cout << "\tappending to cc list" << std::endl;
        v->cc = ( CC* ) realloc( v->cc, ( v->cc_size + 1 ) * sizeof( CC ) );
        std::memmove( &v->cc[ v->cc_size ], cc, sizeof( CC ) );
        //std::cout << "inserting into " << v->cc_size << std::endl;
        v->cc_size++;
        //free_cc( cc );
        free( cc );
    }

    //std::cout << "\tnum of ccs: " << v->cc_size << "\t" << depth << std::endl;
    
    free_uc( &( v->uc ) );
    init_uc( &( v->uc ) );
    //std::cout << "done bursting!" << std::endl;
}

void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
{
    //std::cout << "uc size before: " << v->uc.size << std::endl;
    /*char* seq = deserialize_kmer( k, calc_bk( k ), bseq );
    for( int i = 0; i < depth; i++ )
    {
        std::cout << "  ";
    }
    std::cout << "vertex inserting kmer: " << seq << std::endl;
    free( seq );*/

    if( uc_contains( &( v->uc ), k, depth, bseq ) )
    {
        return;
    }
    
    for( int i = 0; i < v->cc_size; i++ )
    {
        //std::cout << "\tchecking cc: " << i << std::endl;
        CC* cc = &v->cc[ i ];
        if( cc_may_contain( cc, bseq ) )
        {
            //std::cout << "\t\tputting into cc: " << i << std::endl;
            //std::cout << "\t\tmay contain" << std::endl;
            uint8_t* suffix = &bseq[ 1 ];
            if( cc_contains_prefix( cc, bseq ) )
            {
                //std::cout << "\t\t\tcontains!" << std::endl;
            }
            else
            {
                //std::cout << "\t\t\tfalse positive!" << std::endl;
                cc_insert( cc, k, depth, bseq );
            }
            Vertex* child = get_child_of( cc, bseq );
            vertex_insert( child, suffix, k - 4, depth + 1 );
            return;
        }
    }

    if( v->uc.size < CAPACITY - 1 )
    {
        //std::cout << "\t\tinserting into uc!" << std::endl;
        uc_insert( &( v->uc ), bseq, k, depth );
    }
    else
    {
        //std::cout << "about to burst!" << std::endl;
        uc_insert( &( v->uc ), bseq, k, depth );
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

