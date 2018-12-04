#include "Vertex.h"


void init_vertex( Vertex* v )
{
    v->start = false;
    v->pref_pres = (uint256_t) 0;
    v->vs_size = 0;
    v->vs = NULL;
    init_uc( &( v->uc ) );
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

    if(v->vs != NULL) {
        uint8_t prefix = bseq[0];
        uint256_t pre_verts = v->pref_pres;
        pre_verts <<= (256 - (unsigned) prefix);
        pre_verts >>= (256 - (unsigned) prefix);

        if((v->pref_pres >> (unsigned) prefix) & 0x1) {
            //int vidx = __builtin_popcount(pre_verts);
            uint32_t* t = (uint32_t*) &pre_verts;
            int vidx = 0;
            for(int i = 0; i < 8; i++) {
                vidx += __builtin_popcount(t[i]);
            }
            Vertex* child = &v->vs[vidx];
            return vertex_get(child, &bseq[1], k - 4, depth + 1);
        }
    }

    throw pybind11::key_error( "Key not in dictionary!" );
}
#endif

void vertex_remove( Vertex* v, uint8_t* bseq, int k, int depth )
{
    //std::cout << "removing: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
    //std::cout << "\tnum uc items: " << v->uc.size << std::endl;
    for(int i = 0; i < v->uc.size; i++) {
        int idx = i * calc_bk(k);
        //std::cout << "\t\t" << deserialize_kmer(k, calc_bk(k), &v->uc.suffixes[idx]) << std::endl;
    }

    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
        //std::cout << "\tfound item in uc!" << std::endl;
        uc_remove( &( v->uc ), calc_bk( k ), uc_idx );
        return;
    }

    if(v->vs != NULL) {
        uint8_t prefix = bseq[0];
        //std::cout << "\tchecking prefix: " << (unsigned) prefix << std::endl;

        if((v->pref_pres >> (unsigned) prefix) & 0x1) {
            uint256_t pre_verts = v->pref_pres;
            pre_verts <<= (256 - (unsigned) prefix);
            pre_verts >>= (256 - (unsigned) prefix);
            //int vidx = __builtin_popcount(pre_verts);
            uint32_t* t = (uint32_t*) &pre_verts;
            int vidx = 0;
            for(int i = 0; i < 8; i++) {
                vidx += __builtin_popcount(t[i]);
            }
            //std::cout << "\t\tfound vertex to traverse " << vidx << std::endl;
            Vertex* child = &v->vs[vidx];
            vertex_remove(child, &bseq[1], k - 4, depth + 1);
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
            int vidx = 0;
            uint32_t* t = (uint32_t*) &exists;

            for(int i = 0; i < 8; i++) {
                vidx += __builtin_popcount(t[i]);
            }
            //int vidx = __builtin_popcount(exists);
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
        //char* dseq = deserialize_kmer(k, calc_bk(k), &v->uc.suffixes[idx]);
        //std::cout << i << "\t" << dseq << std::endl;

        uint8_t* bseq = &suffixes[ idx ];
        uint8_t prefix = bseq[ 0 ];
        uint8_t* suffix = &bseq[ 1 ];
        uint8_t bits_to_shift = (unsigned) prefix;

        uint256_t pre_verts = v->pref_pres;
        pre_verts <<= (256 - (unsigned) bits_to_shift);
        pre_verts >>= (256 - (unsigned) bits_to_shift);

        //std::cout << "prefix: " << (unsigned) prefix << "\t" << (unsigned) bits_to_shift << std::endl;
        uint32_t* t = (uint32_t*) &pre_verts;
        int vidx = 0;
        for(int i = 0; i < 8; i++) {
            vidx += __builtin_popcount(t[i]);
        }

        // check if there is already a vertex that represents this prefix
        if(!((v->pref_pres >> (unsigned) bits_to_shift) & 0x1)) {
            //std::cout << "vidx: " << vidx << "\t" << (unsigned) prefix << "\t" << deserialize_kmer(4,1,&prefix)<< std::endl;
            //std::cout << "\tvertex does not exist, place at " << vidx << std::endl;
            //std::cout << "\tpref pres: " << v->pref_pres << std::endl;

            if(v->vs == NULL) {
                //std::cout << "\t\tcreating new vertex at position 0" << std::endl;
                v->vs = (Vertex*) calloc(1, sizeof(Vertex));
            } else {
                //std::cout << "\t\tcreating space for a new vertex, new size " << v->vs_size + 1 << std::endl;
                v->vs = (Vertex*) realloc(v->vs, (v->vs_size + 1) * sizeof(Vertex));
            }
            // move any previous vertices if necessary
            if(vidx < v->vs_size) {
                //std::cout << (unsigned) prefix << std::endl;
                //std::cout << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
                //std::cout << "moving vertices " << v->vs_size << "\t" << vidx << "\t" << "\t" << (v->vs_size - vidx)<< std::endl;
                std::memmove(&v->vs[vidx + 1],
                             &v->vs[vidx],
                             (v->vs_size - vidx) * sizeof(Vertex)
                             );
            }

            //std::cout << "old pref pres: " << v->pref_pres << "\t" << (unsigned) bits_to_shift << std::endl;
            v->pref_pres |= ((uint256_t) 0x1 << (unsigned) bits_to_shift);
            //std::cout << "new pref pres: " << v->pref_pres << std::endl;
            //std::cout << " *******************************" << std::endl;

            // insert a vertex at vidx
            //std::cout << "\t\t init vertex at " << vidx << std::endl;
            init_vertex(&v->vs[vidx]);
            v->vs_size++;
    /*for(int i = 0; i < v->vs_size; i++) {
        std::cout << "vertex " << i << std::endl;
        for(int j = 0; j < v->vs[i].uc.size; j++) {
            int idx = j * calc_bk(k - 4);
            std::cout << "\t" << deserialize_kmer(k - 4, calc_bk(k - 4), &v->vs[i].uc.suffixes[idx]) << std::endl;
        }
    }
    std::cout << "done!" << std::endl;*/

            // increment size
        } else {
            //std::cout << "\t" << vidx << "\tusing existing vertex" << std::endl;
        }

        Vertex* child = &v->vs[vidx];
        //std::cout << "vidx: " << vidx << std::endl;

#if KDICT
        vertex_insert( child, suffix, k - 4, depth + 1, &objs[ i ] );
#elif KSET
        vertex_insert( child, suffix, k - 4, depth + 1 );
#endif
    }

    //std::cout << "about to free uc" << std::endl;
    //v->uc.size = 0;
    //std::cout << "done bursting\n\n\n" << std::endl;
    /*for(int i = 0; i < v->vs_size; i++) {
        std::cout << "vertex " << i << std::endl;
        for(int j = 0; j < v->vs[i].uc.size; j++) {
            int idx = j * calc_bk(k - 4);
            std::cout << "\t" << deserialize_kmer(k - 4, calc_bk(k - 4), &v->vs[i].uc.suffixes[idx]) << std::endl;
        }
    }*/
    //std::cout << "done bursting" << std::endl;
    free_uc( &( v->uc ) );
    init_uc( &( v->uc ) );
}

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::handle* obj )
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
#endif
{
    //std::cout << "vertex insertion: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
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

    uint8_t prefix = bseq[ 0 ];
    uint256_t pre_verts = v->pref_pres;
    pre_verts <<= (256 - (unsigned) prefix);
    pre_verts >>= (256 - (unsigned) prefix);
    if((v->pref_pres >> (unsigned) prefix) & 0x1) {
        //int vidx = __builtin_popcount(pre_verts);
        uint32_t* t = (uint32_t*) &pre_verts;
        int vidx = 0;
        for(int i = 0; i < 8; i++) {
            vidx += __builtin_popcount(t[i]);
        }
        Vertex* child = &v->vs[vidx];
#if KDICT
        vertex_insert(child, &bseq[1], k - 4, depth + 1, obj);
#elif KSET
        vertex_insert(child, &bseq[1], k - 4, depth + 1);
#endif
        return;
    }

    if( v->uc.size < CAPACITY - 1 )
    {
        //std::cout << "*********************************" << std::endl;
        //std::cout << &v->uc << std::endl;
        //std::cout << "\tuc insert: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
        //std::cout << "\tbefore uc size: " << v->uc.size << "\t" << uc_idx << std::endl;
#if KDICT
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
        uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#endif
        //std::cout << "\tafter uc size: " << v->uc.size << std::endl;
        //std::cout << "\t\tlast kmer added: ";
        //std::cout << deserialize_kmer(k, calc_bk(k), &v->uc.suffixes[0]) << std::endl;
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
        free_vertex(&v->vs[i]);
    }
    free(v->vs);
}

uint64_t vertex_size( Vertex* v )
{
    uint64_t c = v->uc.size;
    //std::cout << &v->uc << ":" << v->uc.size << std::endl;
    for(int i = 0; i < v->vs_size; i++) {
        c += vertex_size(&v->vs[i]);
    }

    return c;
}

