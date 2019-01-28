#include "Vertex.h"


void init_vertex( Vertex* v )
{
    v->pref_pres = (uint256_t) 0;
    v->vs_size = 0;
    v->vs = NULL;
    init_uc( &( v->uc ) );
}

PgData* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth )
{
    uint8_t prefix = bseq[0];
    //std::cout << "prefix: " << (unsigned) prefix << std::endl;

    if((v->pref_pres >> (unsigned) prefix) & 0x1) {
        int vidx = calc_vidx(v->pref_pres, prefix);
        Vertex* child = &v->vs[vidx];
        return vertex_get(child, &bseq[1], k - 4, depth + 1);
    }

    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    //std::cout << "found: " << sres.first << std::endl;
    if( sres.first )
    {
        int uc_idx = sres.second;
        return &v->uc.data[ uc_idx ];
    }

    /*std::cout << "trying to find: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
    for(int i = 0; i < v->uc.size; i++) {
        int idx = calc_bk(k) * i;
        std::cout << i << "\t" << deserialize_kmer(calc_bk(k) * 4, calc_bk(k), &v->uc.suffixes[idx]) << std::endl;
    }
    sres = binary_search_debug(v->uc.suffixes, v->uc.size, k, bseq);

    std::cout << "Error: could not find " << deserialize_kmer(calc_bk(k) * 4, calc_bk(k), bseq) << std::endl;*/
    throw pybind11::key_error( "Key not in dictionary!" );
}

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
    uint8_t prefix = bseq[0];
    //std::cout << "\tvs exists!\t" << v->vs_size << "\t" << (unsigned) prefix << std::endl;
    if((v->pref_pres >> (unsigned) prefix) & 0x1) {
	// get child
	//std::cout << "found it!" << std::endl;
	int vidx = calc_vidx(v->pref_pres, prefix);
	Vertex* child = &v->vs[vidx];
	return vertex_contains(child, &bseq[1], k - 4, depth + 1);
    }

    //std::cout << "checking vs " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    if( sres.first )
    {
        return true;
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
    PgData* data = v->uc.data;
#if KDICT
    py::handle* objs = v->uc.objs;
#elif KCOUNTER
    count_dtype* counts = v->uc.counts;
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
        vertex_insert( child, suffix, k - 4, depth + 1, (void*) &data[i], true);
#elif KCOUNTER
        vertex_insert( child, suffix, k - 4, depth + 1, counts[ i ] );
#endif
    }

    free_uc( &( v->uc ) );
    init_uc( &( v->uc ) );
}

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::handle* obj )
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, void* data, bool bursting)
#elif KCOUNTER
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, count_dtype count )
#endif
{
    uint8_t prefix = bseq[ 0 ];
    if((v->pref_pres >> (unsigned) prefix) & 0x1) {
        int vidx = calc_vidx(v->pref_pres, prefix);
        Vertex* child = &v->vs[vidx];
#if KDICT
        vertex_insert( child, &bseq[1], k - 4, depth + 1, obj );
#elif KSET
        vertex_insert( child, &bseq[1], k - 4, depth + 1, data, bursting);
#elif KCOUNTER
        vertex_insert( child, &bseq[1], k - 4, depth + 1, count );
#endif
        return;
    }

    // NOTE: the kmer already exists!
    std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
    int uc_idx = sres.second;
    if( sres.first )
    {
        // set the object here

        // NOTE check if genome already seen
        PgData* kvals = &v->uc.data[uc_idx];
        std::tuple<uint16_t, uint32_t, bool, uint8_t*>* cdata = (std::tuple<uint16_t, uint32_t, bool, uint8_t*>*) data;
        uint16_t gidx = std::get<0>(*cdata);
        if(kvals->genomes & 0x1 << gidx) {
            if((unsigned) kvals->counts->back() == 1) {
                kvals->counts->back() = 2;
            }
        } else {
            kvals->genomes |= 0x1 << gidx;
            kvals->counts->push_back(1);
            uint32_t pos = std::get<1>(*cdata);
            bool reverse = std::get<2>(*cdata);
            kvals->coords->push_back(pos);
            kvals->orientation |= reverse << kvals->size;
            kvals->size += 1;
        }
        return;
    }
    uc_insert( &( v->uc ), bseq, k, depth, uc_idx, data, bursting);

    if(v->uc.size == CAPACITY)
    {
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
    for(int i = 0; i < v->vs_size; i++) {
        c += vertex_size(&v->vs[i]);
    }

    return c;
}

