#include "Vertex.h"


void init_vertex( Vertex* v )
{
    v->start = false;
    v->pref_pres = 0;
    v->vs_size = 0;
    v->vs = NULL;
    v->cc = NULL;
    v->cc_size = 0;
    init_uc( &( v->uc ) );
}

CC* get_cc( Vertex* v, int idx )
{
    return &v->cc[ idx ];
}

#if KDICT
py::handle* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth )
{
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
        return &v->uc.objs[ uc_idx ];
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
                    return vertex_get( child, &bseq[ 1 ], k - 4, depth + 1 );
                }
            }
        }
    }

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
    
    //throw pybind11::key_error( "Key not found!" );
}

bool vertex_contains( Vertex* v, uint8_t* bseq, int k, int depth )
{
    //std::cout << "checking vs " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    if( sres.first )
    {
        return true;
    }

    if(v->vs_size > 0) {
        uint8_t prefix = bseq[0];
        uint256_t exists = v->pref_pres;
        exists <<= (256 - (unsigned) prefix);
        exists >>= (256 - (unsigned) prefix);
        //std::cout << "\tvs exists!\t" << v->vs_size << "\t" << (unsigned) prefix << std::endl;
        if((v->pref_pres >> (unsigned) prefix) & 0x1) {
            // get child
            //std::cout << "found it!" << std::endl;
            int vidx = __builtin_popcount(exists);
            Vertex* child = &v->vs[vidx];
            return vertex_contains(child, &bseq[1], k - 4, depth + 1);
        } else {
            //std::cout << "couldn't find a vs" << std::endl;
            return false;
        }
    }

    return false;
}

void burst_uc( Vertex* v, int k, int depth )
{
    //std::cout << "bursting!" << std::endl;
    int suffix_size = calc_bk( k );

    //Vertex* nv = (Vertex*) malloc(sizeof(Vertex));
    //init_vertex(nv);

    uint8_t* suffixes = v->uc.suffixes;
#if KDICT
    py::handle* objs = v->uc.objs;
#endif
    int idx;
    for( int i = 0; i < v->uc.size; i++ )
    {
        //std::cout << "*******************************************" << std::endl;
        idx = i * suffix_size;
        char* dseq = deserialize_kmer(k, calc_bk(k), &v->uc.suffixes[idx]);
        //std::cout << i << "\t" << dseq << std::endl;

        uint8_t* bseq = &suffixes[ idx ];
        uint8_t prefix = bseq[ 0 ];
        uint8_t* suffix = &bseq[ 1 ];
        uint8_t bits_to_shift = (unsigned) prefix;

        uint256_t pre_verts = v->pref_pres;
        pre_verts <<= (256 - (unsigned) bits_to_shift);
        pre_verts >>= (256 - (unsigned) bits_to_shift);

        //std::cout << "prefix: " << (unsigned) prefix << "\t" << (unsigned) bits_to_shift << std::endl;
        int vidx = __builtin_popcount(pre_verts);
        //std::cout << "vidx: " << vidx << std::endl;

        // check if there is already a vertex that represents this prefix
        if(!((v->pref_pres >> (unsigned) bits_to_shift) & 0x1)) {
            //std::cout << "\tvertex does not exist" << std::endl;
            //std::cout << "\tpref pres: " << v->pref_pres << std::endl;

            if(v->vs == NULL) {
                v->vs = (Vertex*) calloc(1, sizeof(Vertex));
            } else {
                v->vs = (Vertex*) realloc(v->vs, (v->vs_size + 1) + sizeof(Vertex));
            }
            // move any previous vertices if necessary
            if(vidx < v->vs_size) {
                std::memmove(&v->vs[vidx + 1],
                             &v->vs[vidx],
                             (v->vs_size - vidx) * sizeof(Vertex)
                             );
            }

            //std::cout << "old pref pres: " << v->pref_pres << "\t" << (unsigned) bits_to_shift << std::endl;
            v->pref_pres |= ((uint256_t) 0x1 << (unsigned) bits_to_shift);
            //std::cout << "new pref pres: " << v->pref_pres << std::endl;

            // insert a vertex at vidx
            init_vertex(&v->vs[vidx]);

            // increment size
            v->vs_size++;
        } else {
            //std::cout << "\tusing existing vertex" << std::endl;
        }

        Vertex* child = &v->vs[vidx];

#if KDICT
        vertex_insert( child, suffix, k - 4, depth + 1, &objs[ i ] );
#elif KSET
        vertex_insert( child, suffix, k - 4, depth + 1 );
#endif
    }

    //std::cout << "about to free uc" << std::endl;
    free_uc( &( v->uc ) );
    init_uc( &( v->uc ) );
    //v->uc.size = 0;
    //std::cout << "done bursting" << std::endl;
}

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::handle* obj )
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
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
#endif
            return;
        }
    }

    if( v->uc.size < CAPACITY - 1 )
    {
        //std::cout << "\tbefore uc size: " << v->uc.size << std::endl;
        //std::cout << "\tuc insert: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
#if KDICT
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#endif
        //std::cout << "\tafter uc size: " << v->uc.size << std::endl;
    }
    else
    {
#if KDICT
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#endif
        burst_uc( v, k, depth );
    }
}

void free_vertex( Vertex* v )
{
    free_uc( &( v->uc ) );
    for(int i = 0; i < v->vs_size; i++) {
        //std::cout << "\tremove vertex " << i << std::endl;
        free_vertex(&v->vs[i]);
    }
}

uint64_t vertex_size( Vertex* v )
{
    uint64_t c = v->uc.size;
    //std::cout << "vs size: " << v->vs_size << std::endl;
    for(int i = 0; i < v->vs_size; i++) {
        c += vertex_size(&v->vs[i]);
    }

    return c;
}

