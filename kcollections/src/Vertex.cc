#include "Vertex.h"


void init_vertex( Vertex* v )
{
  v->start = false;
  v->pref_pres = (uint256_t) 0;
  v->vs_size = 0;
  v->vs = NULL;
  v->uc = new UC<int>();
}

#if defined(KDICT) || defined(KCOUNTER)
int vertex_get( Vertex* v, uint8_t* bseq, int k, int depth )
{
  uint8_t prefix = bseq[0];

  if((v->pref_pres >> (unsigned) prefix) & 0x1) {
    int vidx = calc_vidx(v->pref_pres, prefix);
    Vertex* child = &v->vs[vidx];
    return vertex_get(child, &bseq[1], k - 4, depth + 1);
  }

  std::pair< bool, int > sres = v->uc->uc_find(k, bseq);
  int uc_idx = sres.second;
  if(sres.first) {
    return v->uc->get_obj(uc_idx);
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
  uint8_t prefix = bseq[0];
  //std::cout << "\tchecking prefix: " << (unsigned) prefix << std::endl;

  if((v->pref_pres >> (unsigned) prefix) & 0x1) {
    int vidx = calc_vidx(v->pref_pres, prefix);
    //std::cout << "\t\tfound vertex to traverse " << vidx << std::endl;
    Vertex* child = &v->vs[vidx];
    vertex_remove(child, &bseq[1], k - 4, depth + 1);
  }

  //std::cout << "removing: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
  //std::cout << "\tnum uc items: " << v->uc->size << std::endl;
  //for(int i = 0; i < v->uc->get_size(); i++) {
  //int idx = i * calc_bk(k);
    //std::cout << "\t\t" << deserialize_kmer(k, calc_bk(k), &v->uc->suffixes[idx]) << std::endl;
  //}

  std::pair< bool, int > sres = v->uc->uc_find(k, bseq);
  int uc_idx = sres.second;
  if( sres.first )
    {
      //std::cout << "\tfound item in uc!" << std::endl;
      v->uc->uc_remove(calc_bk(k), uc_idx);
      return;
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
  std::pair< bool, int > sres = v->uc->uc_find(k, bseq);
  if( sres.first )
    {
      return true;
    }

  return false;
}

int calc_vidx(uint256_t vertices, uint8_t bts) {
  vertices <<= (256 - (unsigned) bts);
  uint64_t* t = (uint64_t*) &vertices;
  int vidx = __builtin_popcountll(t[0]);
  vidx += __builtin_popcountll(t[1]);
  vidx += __builtin_popcountll(t[2]);
  vidx += __builtin_popcountll(t[3]);
  //vidx += __builtin_popcount(t[4]);
  //vidx += __builtin_popcount(t[5]);
  //vidx += __builtin_popcount(t[6]);
  //vidx += __builtin_popcount(t[7]);
  return vidx;
}

#if defined(KDICT) || defined(KCOUNTER)
void burst_uc( Vertex* v, int k, int depth, std::function<int(int, int)>* merge_func )
#else
void burst_uc( Vertex* v, int k, int depth )
#endif
{
  //std::cout << "bursting!" << std::endl;
  int suffix_size = calc_bk( k );

  //Vertex* nv = (Vertex*) malloc(sizeof(Vertex));
  //init_vertex(nv);

  uint8_t* suffixes = v->uc->get_suffixes();
#if defined(KDICT) || defined(KCOUNTER)
  std::vector<int> objs = v->uc->get_objs();
#endif
  int idx;
  for( int i = 0; i < v->uc->get_size(); i++ )
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
#if defined(KDICT) || defined(KCOUNTER)
      vertex_insert( child, suffix, k - 4, depth + 1, objs[ i ], merge_func );
#elif KSET
      vertex_insert( child, suffix, k - 4, depth + 1 );
#endif
    }

  delete v->uc;
  v->uc = new UC<int>();
}

#if defined(KDICT) || defined(KCOUNTER)
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, int obj, std::function<int(int, int)>* merge_func )
#elif KSET
  void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
#endif
{
  //std::cout << "vertex insertion: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
  uint8_t prefix = bseq[ 0 ];
  if((v->pref_pres >> (unsigned) prefix) & 0x1) {
    int vidx = calc_vidx(v->pref_pres, prefix);
    Vertex* child = &v->vs[vidx];
#if defined(KDICT) || defined(KCOUNTER)
    vertex_insert( child, &bseq[1], k - 4, depth + 1, obj, merge_func );
#elif KSET
    vertex_insert( child, &bseq[1], k - 4, depth + 1 );
#endif
    return;
  }

  // NOTE: if the key already exisrts, we need to update it somehow
  std::pair< bool, int > sres = v->uc->uc_find(k, bseq);
  int uc_idx = sres.second;
  if( sres.first ) {
    // replace object here
#if defined(KDICT) || defined(KCOUNTER)
    if(merge_func != NULL) {
      int merged_obj = (*merge_func)(v->uc->get_obj(uc_idx), obj);
      v->uc->set_obj(uc_idx, merged_obj);
    } else {
      v->uc->set_obj(uc_idx, obj);
    }
#endif
    return;
  }

#if defined(KDICT) || defined(KCOUNTER)
  v->uc->uc_insert(bseq, k, uc_idx, obj);
#elif KSET
  v->uc->uc_insert(bseq, k, uc_idx);
#endif

  if(v->uc->get_size() == CAPACITY) {
#if defined(KDICT) || defined(KCOUNTER)
    burst_uc(v, k, depth, merge_func);
#else
    burst_uc(v, k, depth);
#endif
  }
}

void free_vertex( Vertex* v )
{
  delete v->uc;
  for(int i = 0; i < v->vs_size; i++) {
    free_vertex(&v->vs[i]);
  }
  free(v->vs);
}

uint64_t vertex_size( Vertex* v )
{
  uint64_t c = v->uc->get_size();
  //std::cout << &v->uc << ":" << v->uc->size << std::endl;
  for(int i = 0; i < v->vs_size; i++) {
    c += vertex_size(&v->vs[i]);
  }

  return c;
}

