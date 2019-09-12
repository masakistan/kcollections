#include "Vertex.h"


void init_vertex( Vertex* v )
{
  v->start = false;
  v->pref_pres = (uint256_t) 0;
  v->vs_size = 0;
  v->vs = NULL;
  init_uc( &( v->uc ) );
}

#if defined KDICT || defined KCOUNTER
#if KDICT
py::object* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOUNTER
  int vertex_get_counter( Vertex* v, uint8_t* bseq, int k, int depth )
#endif
{
  uint8_t prefix = bseq[0];

  if((v->pref_pres >> (unsigned) prefix) & 0x1) {
    int vidx = calc_vidx(v->pref_pres, prefix);
    Vertex* child = &v->vs[vidx];
#if KDICT
    return vertex_get(child, &bseq[1], k - 4, depth + 1);
#elif KCOUNTER
    return vertex_get_counter( child, &bseq[ 1 ], k - 4, depth + 1 );
#endif
  }

  std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
  int uc_idx = sres.second;
  if( sres.first )
    {
#if KDICT
      return &v->uc.objs[ uc_idx ];
#elif KCOUNTER
      return v->uc.counts[ uc_idx ];
#endif
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

#if defined(KDICT)
void burst_uc( Vertex* v, int k, int depth, std::function<py::object(py::object, py::object)>* merge_func )
#else
void burst_uc( Vertex* v, int k, int depth )
#endif
{
  //std::cout << "bursting!" << std::endl;
  int suffix_size = calc_bk( k );

  //Vertex* nv = (Vertex*) malloc(sizeof(Vertex));
  //init_vertex(nv);

  uint8_t* suffixes = v->uc.suffixes;
#if KDICT
  py::object* objs = v->uc.objs;
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
      vertex_insert( child, suffix, k - 4, depth + 1, objs[ i ], merge_func );
#elif KSET
      vertex_insert( child, suffix, k - 4, depth + 1 );
#elif KCOUNTER
      vertex_insert( child, suffix, k - 4, depth + 1, counts[ i ] );
#endif
    }

  free_uc( &( v->uc ) );
  init_uc( &( v->uc ) );
}

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::object obj, std::function<py::object(py::object, py::object)>* merge_func )
#elif KSET
  void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth )
#elif KCOUNTER
  void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, count_dtype count )
#endif
{
  //std::cout << "vertex insertion: " << deserialize_kmer(k, calc_bk(k), bseq) << std::endl;
  uint8_t prefix = bseq[ 0 ];
  if((v->pref_pres >> (unsigned) prefix) & 0x1) {
    int vidx = calc_vidx(v->pref_pres, prefix);
    Vertex* child = &v->vs[vidx];
#if KDICT
    vertex_insert( child, &bseq[1], k - 4, depth + 1, obj, merge_func );
#elif KSET
    vertex_insert( child, &bseq[1], k - 4, depth + 1 );
#elif KCOUNTER
    vertex_insert( child, &bseq[1], k - 4, depth + 1, count );
#endif
    return;
  }

  // NOTE: if the key already exisrts, we need to update it somehow
  std::pair< bool, int > sres = uc_find( &( v->uc ), k, depth, bseq );
  int uc_idx = sres.second;
  if( sres.first ) {
    //std::cout << "found previously!" << std::endl;
      // replace object here
#if KDICT
      //py::object merged_obj = obj;
      if(merge_func != NULL) {
	//std::cout << "merging!" << std::endl;
	//return;
	py::gil_scoped_acquire acquire;
	  py::object merged_obj = (*merge_func)(v->uc.objs[uc_idx], obj);
	  //std::cout << "merged: " << std::string(py::str(merged_obj)) << std::endl;
	  //obj.dec_ref();
	  v->uc.objs[uc_idx].dec_ref();
	  merged_obj.inc_ref();
	  py::gil_scoped_release release;

	std::memcpy(
		  &v->uc.objs[ uc_idx ],
		  &merged_obj,
		  sizeof( py::object )
		  );
	//std::cout << "output" << std::string(py::str(*(v->uc.objs[uc_idx]))) << std::endl;
	//py::gil_scoped_acquire acquire1;
	//v->uc.objs[ uc_idx ].inc_ref();
	//py::gil_scoped_release release1;
	//std::cout << "done merging" << std::endl;
	
	return;
      } else {
	std::cout << "not merging\t" << merge_func << std::endl;
      }

      {
	py::gil_scoped_acquire acquire;
	v->uc.objs[ uc_idx ].dec_ref();
	obj.inc_ref();
	py::gil_scoped_release release;
      }
      std::memcpy(
		  &v->uc.objs[ uc_idx ],
		  &obj,
		  sizeof( py::object )
		  );
      {
	//py::gil_scoped_acquire acquire;
	//v->uc.objs[ uc_idx ].inc_ref();
	//py::gil_scoped_release release;
      }
#elif KCOUNTER
      if(v->uc.counts[uc_idx] < MAXCOUNT) {
	std::memcpy(
		    &v->uc.counts[ uc_idx ],
		    &count,
		    sizeof(count_dtype)
		    );
      }
#endif
      return;
    }

#if KDICT
  uc_insert( &( v->uc ), bseq, k, depth, uc_idx, obj );
#elif KSET
  uc_insert( &( v->uc ), bseq, k, depth, uc_idx );
#elif KCOUNTER
  uc_insert( &( v->uc ), bseq, k, depth, uc_idx, count );
#endif

  if(v->uc.size == CAPACITY)
    {
      //std::cout << "bursting " << v->uc.size << std::endl;
#if defined(KDICT)      
      burst_uc( v, k, depth, merge_func );
#else
      burst_uc( v, k, depth );
#endif
      //std::cout << "\tafter" << v->uc.size << std::endl;
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

