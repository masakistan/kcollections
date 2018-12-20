#include "Vertex.h"


void init_vertex( Vertex* v )
{
    //v->start = false;
    v->pref_pres = (uint256_t) 0;
    v->vs_size = 0;
    v->vs = NULL;
    init_uc( &( v->uc ) );
}

#if defined KDICT || defined KCOUNTER || KCOLOR
#if KDICT
py::handle* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOUNTER
int vertex_get_counter( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOLOR
uint32_t* vertex_get_colors(Vertex* v, uint8_t* bseq, int k, int depth)
#endif
{
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
#if KDICT
        return &v->uc.objs[ uc_idx ];
#elif KCOUNTER
        return v->uc.counts[ uc_idx ];
#elif KCOLOR
        uint64_t card = roaring_bitmap_get_cardinality(v->uc.colors[uc_idx]);
        uint32_t* colors = (uint32_t*) malloc(card * sizeof(uint32_t));
        roaring_bitmap_to_uint32_array(v->uc.colors[uc_idx], colors);
        return colors;
#endif
    }

    if(v->vs != NULL) {
        uint8_t prefix = bseq[0];

        if((v->pref_pres >> (unsigned) prefix) & 0x1) {
            int vidx = calc_vidx(v->pref_pres, prefix);
            Vertex* child = &v->vs[vidx];
#if KDICT
            return vertex_get(child, &bseq[1], k - 4, depth + 1);
#elif KCOUNTER
            return vertex_get_counter( child, &bseq[ 1 ], k - 4, depth + 1 );
#elif KCOLOR
            return vertex_get_colors(child, &bseq[1], k - 4, depth + 1);
#endif
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
            int vidx = calc_vidx(v->pref_pres, prefix);
            //std::cout << "\t\tfound vertex to traverse " << vidx << std::endl;
            Vertex* child = &v->vs[vidx];
            vertex_remove(child, &bseq[1], k - 4, depth + 1);
        }
    }

    throw pybind11::key_error( "Key not found!" );
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
        //std::cout << "\tvs exists!\t" << v->vs_size << "\t" << (unsigned) prefix << std::endl;
        if((v->pref_pres >> (unsigned) prefix) & 0x1) {
            // get child
            //std::cout << "found it!" << std::endl;
            int vidx = calc_vidx(v->pref_pres, prefix);
            Vertex* child = &v->vs[vidx];
            return vertex_contains(child, &bseq[1], k - 4, depth + 1);
        } else {
            //std::cout << "couldn't find a vs" << std::endl;
            return false;
        }
    }

    return false;
}

int calc_vidx(uint256_t vertices, uint8_t bts) {
    vertices <<= (256 - (unsigned) bts);
    uint32_t* t = (uint32_t*) &vertices;
    int vidx = __builtin_popcount(t[0]);
    vidx += __builtin_popcount(t[1]);
    vidx += __builtin_popcount(t[2]);
    vidx += __builtin_popcount(t[3]);
    vidx += __builtin_popcount(t[4]);
    vidx += __builtin_popcount(t[5]);
    vidx += __builtin_popcount(t[6]);
    vidx += __builtin_popcount(t[7]);
    return vidx;
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
#elif KCOUNTER
    int* counts = v->uc.counts;
#elif KCOLOR
    roaring_bitmap_t** colors = v->uc.colors;
#endif
    int idx;
    for( int i = 0; i < v->uc.size; i++ )
    {
        idx = i * suffix_size;

        uint8_t* bseq = &suffixes[ idx ];
        uint8_t prefix = bseq[ 0 ];
        uint8_t* suffix = &bseq[ 1 ];
        uint8_t bits_to_shift = (unsigned) prefix;
        int vidx = calc_vidx(v->pref_pres, prefix);
        //std::cout << "burst: " << deserialize_kmer(k, calc_bk(k), bseq) << "\t" << (unsigned) bseq[0] << std::endl;

        // check if there is already a vertex that represents this prefix
        if(!((v->pref_pres >> (unsigned) bits_to_shift) & 0x1)) {

            if(v->vs == NULL) {
                v->vs = (Vertex*) calloc(1, sizeof(Vertex));
            } else {
                v->vs = (Vertex*) realloc(v->vs, (v->vs_size + 1) * sizeof(Vertex));
            }
            // move any previous vertices if necessary
            if(vidx < v->vs_size) {
                std::memmove(&v->vs[vidx + 1],
                             &v->vs[vidx],
                             (v->vs_size - vidx) * sizeof(Vertex)
                             );
            }

            v->pref_pres |= ((uint256_t) 0x1 << (unsigned) bits_to_shift);
            // insert a vertex at vidx
            init_vertex(&v->vs[vidx]);
            // increment size
            v->vs_size++;
        } else {
            //std::cout << "\t" << vidx << "\tusing existing vertex" << std::endl;
        }

        Vertex* child = &v->vs[vidx];
#if KDICT
        vertex_insert( child, suffix, k - 4, depth + 1, &objs[ i ] );
#elif KSET
        vertex_insert( child, suffix, k - 4, depth + 1 );
#elif KCOUNTER
        vertex_insert( child, suffix, k - 4, depth + 1, counts[ i ] );
#elif KCOLOR
        vertex_insert(child, suffix, k - 4, depth + 1, colors[i]);
#endif
    }

#if KCOLOR
    free_uc(&v->uc, true);
#else
    free_uc(&v->uc);
#endif
    init_uc( &( v->uc ) );
}

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::handle* obj )
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOUNTER
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, int count )
#elif KCOLOR
void vertex_insert(Vertex* v, uint8_t* bseq, int k, int depth, roaring_bitmap_t* colors)
#endif
{
    uint8_t prefix = bseq[ 0 ];
    if((v->pref_pres >> (unsigned) prefix) & 0x1) {
        int vidx = calc_vidx(v->pref_pres, prefix);
        Vertex* child = &v->vs[vidx];
#if KDICT
        vertex_insert( child, &bseq[1], k - 4, depth + 1, obj );
#elif KSET
        vertex_insert( child, &bseq[1], k - 4, depth + 1 );
#elif KCOUNTER
        vertex_insert( child, &bseq[1], k - 4, depth + 1, count );
#elif KCOLOR
        vertex_insert(child, &bseq[1], k - 4, depth + 1, colors);
#endif
        return;
    }

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
#elif KCOUNTER
        std::memcpy(
                &v->uc.counts[ uc_idx ],
                &count,
                sizeof( int )
                );
#elif KCOLOR
        roaring_bitmap_or_inplace(v->uc.colors[uc_idx], colors);
        roaring_bitmap_free(colors);
#endif
        return;
    }

#if KDICT
    uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
    uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#elif KCOUNTER
    uc_insert( &( v->uc ), bseq, k, depth, uc_idx, count );
#elif KCOLOR
    uc_insert(&(v->uc), bseq, k, depth, uc_idx, colors);
#endif

    if(v->uc.size == CAPACITY)
    {
        //std::cout << "bursting" << std::endl;
        burst_uc( v, k, depth );
    }
}

void free_vertex( Vertex* v )
{
#if KCOLOR
    free_uc(&v->uc, false);
#else
    free_uc(&v->uc);
#endif
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

