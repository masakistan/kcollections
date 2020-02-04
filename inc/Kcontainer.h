#pragma once

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include "Vertex.h"
#include "helper.h"
//#include <jemalloc/jemalloc.h>
#include <math.h>
#include <functional>

#if defined(PYTHON)
namespace py = pybind11;
#endif

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif

struct ThreadGlobals {
#if defined(KSET)
  std::vector<std::vector<std::vector<uint8_t*>>>* kmers;
#elif defined(KDICT) || defined(KCOUNTER)
  std::vector<std::vector<std::vector<std::pair<uint8_t*, T>>>>* kmers;
#endif

#if defined(KDICT) || defined(KCOUNTER)
  std::function<T(T&, T&)>* merge_func;
  Vertex<T>** v;
#else
  Vertex** v;
#endif

  sem_t* signal_b;
  pthread_mutex_t** blocks;
  sem_t** rsignal;
  int bk, k, nthreads;
  pthread_t* p_threads;
  int* wbin;
  int* rbin;
  int work_queues;
  int bits_to_shift;
  size_t MAX_BIN_SIZE;
};

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif
struct ThreadInfo{
  int thread_id;
#if defined(KDICT) || defined(KCOUNTER)
  ThreadGlobals<T>* tg;
#else
  ThreadGlobals* tg;
#endif
};


#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif
class Kcontainer {

public:
  Kcontainer(const int k) : k(k) {
    tg = NULL;
    ti = NULL;
  }

  ~Kcontainer() {
    if(ti != NULL) {
      free(ti);
    }
    if(tg != NULL) {
      free(tg);
    }
  }

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & k;
    ar & v;
  }

  class KcontainerIterator {
  public:
    typedef std::input_iterator_tag iterator_category;    // iterator category
    typedef std::string value_type;           // value type
    typedef std::ptrdiff_t difference_type;           // difference type

#if defined(KDICT) || defined(KCOUNTER)
    typedef std::pair<std::string, T*>* pointer;
    typedef std::pair<std::string, T*>& reference;
#elif defined(KSET)
    typedef std::string* pointer;    // pointer
    typedef std::string& reference;           // reference
#endif

    KcontainerIterator() {}

#if defined(KDICT) || defined(KCOUNTER)
    KcontainerIterator(Vertex<T>& v, int k) : k(k)
#else
    KcontainerIterator(Vertex& v, int k) : k(k)
#endif
    {
      // check if there is an uncompressed container to start
      depth = 0;
      v_stack.push_back(&v);
      uc_stack.push_back(0);
      cc_stack.push_back(0);
#if defined(KDICT) || defined(KCOUNTER)
      kmer = std::pair<std::string, T*>(std::string(k, 'X'), NULL);
#else
      kmer = std::string(k, 'X');
#endif
      get_next();
      //std::cout << "creating iterator!\t" << kmer << std::endl;
    }
    
    reference operator*() {
      return kmer;
    }
    
    pointer operator->() {
      return &kmer;
    }
    
    KcontainerIterator& operator++() {
      get_next();
      return *this;
    }
    
    KcontainerIterator operator++(int) {
      KcontainerIterator oldValue = *this;
      get_next();
      return oldValue;
    }
    
    bool operator== (KcontainerIterator& right) {
      if(v_stack.size() != right.v_stack.size()) {
	return false;
      }

      if(v_stack.size() == 0) {
	return true;
      }
      
      if(v_stack.back() != right.v_stack.back() || uc_stack.back() != right.uc_stack.back() || cc_stack.back() != right.cc_stack.back()) {
	return false;
      }
      
      return true;
    }
    
    bool operator!= (KcontainerIterator& right) {
      return !(*this == right);
    }
    
  private:
    int depth, k;
#if defined(KDICT) || defined(KCOUNTER)
    std::pair<std::string, T*> kmer;
    std::vector<Vertex<T>*> v_stack;
#else
    std::string kmer;
    std::vector<Vertex*> v_stack;
#endif
    std::vector<int> uc_stack;
    std::vector<int> cc_stack;
    //std::vector<std::string> prefixes;
    
    void get_next() {
      int uc_idx = uc_stack.back();
      int cc_idx = cc_stack.back();
#if defined(KDICT) || defined(KCOUNTER)
      Vertex<T>* v_ptr = v_stack.back();
#else
      Vertex* v_ptr = v_stack.back();
#endif

      int n = k - (depth * 4);
      
      if((size_t) uc_idx < v_ptr->get_uc()->get_size()) {
	// NOTE: remove last n characters
	// NOTE: replace last n characters
#if defined(KDICT) || defined(KCOUNTER)
	kmer.first.replace(k - n, n, deserialize_kmer_to_string(n, v_ptr->get_uc()->get_suffix(n, uc_idx)));
	kmer.second = &v_ptr->get_uc()->get_obj(uc_idx);
#else
	kmer.replace(k - n, n, deserialize_kmer_to_string(n, v_ptr->get_uc()->get_suffix(n, uc_idx)));
#endif
	uc_stack.back()++;
	return;
      } else {
	if(cc_idx < v_ptr->get_vs_size()) {
#if defined(KDICT) || defined(KCOUNTER)
	  kmer.first.replace(k - n, k - n + 4, kcontainer_get_child_suffix(v_ptr, cc_idx));
#else
	  kmer.replace(k - n, k - n + 4, kcontainer_get_child_suffix(v_ptr, cc_idx));
#endif
	  
	  depth++;
	  v_stack.push_back(&v_ptr->get_vs()[cc_idx]);
	  cc_stack.back()++;
	  cc_stack.push_back(0);
	  uc_stack.push_back(0);
	  return get_next();
	}
      }

      depth--;
      v_stack.pop_back();
      cc_stack.pop_back();
      uc_stack.pop_back();

      if(v_stack.size() == 0) {
#if defined(KDICT) || defined(KCOUNTER)
	kmer.first = std::string("");
	kmer.second = NULL;
#else
	kmer = std::string("");
#endif
	return;
      }
      
      return get_next();
    }
  };
  
  typedef KcontainerIterator iterator;
  
#if defined(KDICT) || defined(KCOUNTER)
  Kcontainer<T>::iterator begin() {
    return Kcontainer<T>::iterator(v, k);
  }
  Kcontainer<T>::iterator end() {
    return Kcontainer<T>::iterator();
  }
  Vertex<T>* get_v() { return &v; }
  
#else
  Kcontainer::iterator begin() {
    return Kcontainer::iterator(v, k);
  }
  Kcontainer::iterator end() {
    return Kcontainer::iterator();
  }
  Vertex* get_v() { return &v; }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  static std::string kcontainer_get_child_suffix(Vertex<T>* v, int idx) {
#else
  static std::string kcontainer_get_child_suffix(Vertex* v, int idx) {
#endif
    uint256_t verts = v->get_pref_pres();
    uint8_t j = 0, i = 0;
    while(true) {
      if(verts & 0x1) {
	j++;
      }

      if(j > idx)
	break;

      if((unsigned) i == 255) {
	break;
      }

      verts >>= 1;
      i++;
    }
    char* kmer = deserialize_kmer(4, &i);
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  bool kcontainer_contains(const char* kmer)
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    int ret = serialize_kmer(kmer, k, bseq);
    bool res = false;
    if(ret != -1) {
      free(bseq);
      throw std::invalid_argument("Contains op: Could not serialize kmer, ambiguity bases present.");
    }
    res = v.vertex_contains(bseq, k);
    free(bseq);
    return res;
  }

  void kcontainer_remove(const char* kmer)
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    int ret = serialize_kmer(kmer, k, bseq);
    if(ret != -1) {
      free(bseq);
      throw std::invalid_argument("Remove op: Could not serialize kmer, ambiguity bases present.");
    }
    v.vertex_remove(bseq, k);
    free(bseq);
  }

#if defined(KDICT) || defined(KCOUNTER)
  void kcontainer_add(const char* kmer, T obj, std::function<T(T&, T&)>& merge_func)
#elif KSET
  void kcontainer_add(const char* kmer)
#endif
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    int ret = serialize_kmer(kmer, k, bseq);
    if(ret != -1) {
      free(bseq);
      throw std::invalid_argument("Add op: Could not serialize kmer, ambiguity bases present.");
    }
    
#if defined(KDICT) || defined(KCOUNTER)
    v.vertex_insert(bseq, k, obj, merge_func);
#elif KSET
    v.vertex_insert(bseq, k);
#endif
    free(bseq);
  }

#if defined KDICT || defined KCOUNTER
#if defined KDICT
  T& kcontainer_get(const char* kmer)
#elif defined KCOUNTER
  T kcontainer_get(const char* kmer)
#endif
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    int ret = serialize_kmer(kmer, k, bseq);
    if(ret != -1) {
      free(bseq);
      throw std::invalid_argument("Get op: Could not serialize kmer, ambiguity bases present.");
    } else {
#if defined KDICT
    T& res = v.vertex_get(bseq, k);
#elif defined KCOUNTER
    T res = v.vertex_get(bseq, k);
#endif
    free(bseq);
    return res;
    }
  }
#endif

  uint64_t kcontainer_size()
  {
    return v.get_vertex_size();
  }

#if defined(KDICT)
#if defined(PYTHON)
  void kcontainer_add_seq(const char* seq, uint32_t length, py::iterable& values, std::function<T(T&, T&)> &f)
#else
  template <typename Iterable>
  void kcontainer_add_seq(const char* seq, uint32_t length, Iterable& values, std::function<T(T, T)> &f)
#endif
#elif defined(KCOUNTER)
  void kcontainer_add_seq(const char* seq, uint32_t length, std::function<T(T&, T&)> &f)
#else
  void kcontainer_add_seq(const char* seq, uint32_t length)
#endif
  {
    int size64 = k / 32;
    if(k % 32 > 0) {
      size64++;
    }

    uint64_t* bseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
    uint8_t* bseq8 = (uint8_t*) bseq64;

    uint bk = calc_bk(k);
    uint8_t last_index = (k - 1) % 4;

    // serialize the first kmer
    int ret = serialize_kmer(seq, k, bseq8);
    int start = 0;
    while(ret != -1 && start + ret + k < length) {
      start += ret + 1;
      ret = serialize_kmer(&seq[start], k, bseq8);
    }

#if KSET
    v.vertex_insert(bseq8, k);
#elif defined(KCOUNTER)
    v.vertex_insert(bseq8, k, 1, f);
#elif defined(KDICT)
#if defined(PYTHON)
    auto iter = py::iter(values);
    v.vertex_insert(bseq8, k, (*iter).cast<T>(), f);
#else
    auto iter = values.begin();
    v.vertex_insert(bseq8, k, *iter, f);
#endif
#endif

    for(uint32_t j = start + k; j < length; j++) {
      bseq64[0] >>= 2;
      for(int i = 1; i < size64; i++) {
	bseq64[i - 1] |= (bseq64[i] << 62);
	bseq64[i] >>= 2;
      }

      ret = serialize_position(j, bk - 1, last_index, bseq8, seq);

      if(ret != -1) {
	start = j - k + 1;
	ret = k - 1;

	while(ret != -1 && start + ret + k < length) {
	  // NOTE: can't big shift part of it since a lot of could be wrong
	  start += ret + 1;
	  memset(bseq64, 0, sizeof(uint64_t) * size64);
	  ret = serialize_kmer(&seq[start], k, bseq8);
	}

	if(ret != -1) {
	  break;
	}

	j = start + k - 1;
      }
      
#if KSET
      v.vertex_insert(bseq8, k);
#elif defined(KCOUNTER)
      v.vertex_insert(bseq8, k, 1, f);
#elif defined(KDICT)
#if defined(PYTHON)
      std::advance(iter, 1);
      v.vertex_insert(bseq8, k, (*iter).cast<T>(), f);
#else
      v.vertex_insert(bseq8, k, *iter, f);
#endif
#endif
    }

    free(bseq64);
  }


  void tg_init(int threads) {
    tg->MAX_BIN_SIZE = 500;
    tg->work_queues = 10;

    // NOTE: initialize variables
    tg->nthreads = threads;
    tg->bits_to_shift = 8 - log2((double) threads);
    tg->bk = calc_bk(k);
    tg->k = k;

#if defined(KDICT) || defined(KCOUNTER)
    tg->v = (Vertex<T>**) calloc(threads, sizeof(Vertex<T>*));
#else
    tg->v = (Vertex**) calloc(threads, sizeof(Vertex*));
#endif

    tg->signal_b = (sem_t*) calloc(threads, sizeof(sem_t));
    tg->rsignal = (sem_t**) calloc(threads, sizeof(sem_t*));
    tg->p_threads = (pthread_t*) calloc(threads, sizeof(pthread_t));

    tg->wbin = (int*) calloc(threads, sizeof(int));
    tg->rbin = (int*) calloc(threads, sizeof(int));

    tg->blocks = (pthread_mutex_t**) calloc(threads, sizeof(pthread_mutex_t*));
    std::string appName = "kcollections";

#if defined(KDICT) || defined(KCOUNTER)
    tg->kmers = new std::vector<std::vector<std::vector<std::pair<uint8_t*, T>>>>();
#elif defined(KSET)
    tg->kmers = new std::vector<std::vector<std::vector<uint8_t*>>>();
#endif
    // NOTE: initialize per thread variables
    // NOTE: initialize mutexes
    for(int i = 0; i < threads; i++) {
#if defined(KSET)
      tg->kmers->push_back(std::vector<std::vector<uint8_t*>>());
#elif defined(KDICT) || defined(KCOUNTER)
      tg->kmers->push_back(std::vector<std::vector<std::pair<uint8_t*, T>>>());
#endif

      tg->blocks[i] = (pthread_mutex_t*) calloc(tg->work_queues, sizeof(pthread_mutex_t));
      for(int j = 0; j < tg->work_queues; j++) {
        pthread_mutex_init(&tg->blocks[i][j], NULL);

#if defined(KSET)
        (*tg->kmers)[i].push_back(std::vector<uint8_t*>());
#elif defined(KDICT) || defined(KCOUNTER)
        (*tg->kmers)[i].push_back(std::vector<std::pair<uint8_t*, T>>());
#endif

      }
#if defined(KDICT) || defined(KCOUNTER)
      tg->v[i] = new Vertex<T>();
#else
      tg->v[i] = new Vertex();
#endif

      tg->rsignal[i] = sem_open((appName + std::to_string(getpid()) + std::to_string(i)).c_str(), O_CREAT, 0600, 0);
    }
  }

  void ti_init(int threads) {
    for(int i = 0; i < threads; i++) {
      ti[i].thread_id = i;
      ti[i].tg = tg;
    }
  }
    

#if defined(KSET)
  void parallel_kcontainer_add_init(int threads)
#elif defined(KDICT) || defined(KCOUNTER)
  void parallel_kcontainer_add_init(int threads, const std::function<T(T&, T&)> &f)
#endif
  {
#if defined(KDICT)
    ti = (ThreadInfo<T>*) calloc(threads, sizeof(ThreadInfo<T>));
    tg = (ThreadGlobals<T>*) calloc(1, sizeof(ThreadGlobals<T>));
    tg->merge_func = new std::function<T(T&, T&)>(f);
#elif defined(KCOUNTER)
    ti = (ThreadInfo<T>*) calloc(threads, sizeof(ThreadInfo<T>));
    tg = (ThreadGlobals<T>*) calloc(1, sizeof(ThreadGlobals<T>));
    tg->merge_func = new std::function<T(T&, T&)>(f);
#else
    ti = (ThreadInfo*) calloc(threads, sizeof(ThreadInfo));
    tg = (ThreadGlobals*) calloc(1, sizeof(ThreadGlobals));
#endif

    tg_init(threads);
    ti_init(threads);


    for(int i = 0; i < threads; i++) {
      // NOTE: spin up worker threads
      pthread_create(&tg->p_threads[i], NULL, parallel_kcontainer_add_consumer, &ti[i]);
    }
  }

  void parallel_kcontainer_add_join() {
    //std::cout << "joining threads" << std::endl;

    for(int i = 0; i < tg->nthreads; i++) {
      // NOTE: we post twice to finish working.
      // once to finish any kmers in the pipe
      // and a second time to pass the empty kmer list to
      // the thread to signal work is done
      sem_post(tg->rsignal[i]);
      sem_post(tg->rsignal[i]);
    }

    //std::cout << "posted" << std::endl;

    int total_vs = 0;
    //uint16_t total_kmers = 0;
    for(int i = 0; i < tg->nthreads; i++) {
      pthread_join(tg->p_threads[i], NULL);
      //std::cout << "joined " << i << std::endl;
      //std::cout << "uc size: " << v[i]->uc.size << std::endl;
      total_vs += tg->v[i]->get_vs_size();
      sem_close(tg->rsignal[i]);
    }

#if defined(KDICT) || defined(KCOUNTER)
    v.set_vs(new Vertex<T>[total_vs]);
#else
    v.set_vs(new Vertex[total_vs]);
#endif

    v.set_vs_size(total_vs);

    // NOTE: clean up memory
    int idx = 0;
    for(int i = 0; i < tg->nthreads; i++) {
      //std::cout << "thread: " << i << "\t" << v[i]->vs_size << std::endl;

#if defined(KDICT) || defined(KCOUNTER)
      Vertex<T>* thread_vs = tg->v[i]->get_vs();
#else
      Vertex* thread_vs = tg->v[i]->get_vs();
#endif
      uint16_t thread_vs_size = tg->v[i]->get_vs_size();

      if(thread_vs != NULL) {
	for(uint16_t vs_idx = 0; vs_idx < thread_vs_size; vs_idx++) {
	  v.get_vs()[idx + vs_idx] = std::move(thread_vs[vs_idx]);
	}

	v.set_pref_pres(v.get_pref_pres() | tg->v[i]->get_pref_pres());
	idx += thread_vs_size;
	delete[] thread_vs;
#if defined(KDICT) || defined(KCOUNTER)
	tg->v[i]->set_vs((Vertex<T>*) NULL);
#else
	tg->v[i]->set_vs((Vertex*) NULL);
#endif

      }

      delete tg->v[i];
      free(tg->blocks[i]);
      (*tg->kmers)[i].clear();
    }

    free(tg->v);
    free(tg->signal_b);
    free(tg->rsignal);
    free(tg->p_threads);
    free(tg->wbin);
    free(tg->rbin);
    free(tg->blocks);
    tg->kmers->clear();
    delete tg->kmers;
#if defined(KDICT) || defined(KCOUNTER)
    delete tg->merge_func;
#endif

    free(tg);
    free(ti);

    tg = NULL;
    ti = NULL;
  }

  static void* parallel_kcontainer_add_consumer(void* bin_ptr) {
    //int bin = *((int*) bin_ptr);
#if defined(KDICT) || defined(KCOUNTER)
    ThreadInfo<T>* cti = (ThreadInfo<T>*) bin_ptr;
    ThreadGlobals<T>* ctg = cti->tg;
#else
    ThreadInfo* cti = (ThreadInfo*) bin_ptr;
    ThreadGlobals* ctg = cti->tg;
#endif

    int bin = cti->thread_id;
    int cur_rbin;

    while(true) {
      sem_wait(ctg->rsignal[bin]);
      cur_rbin = ctg->rbin[bin];
      pthread_mutex_lock(&ctg->blocks[bin][cur_rbin]);

      // NOTE: release mutex

      // NOTE: if the list is empty we're done
      if((*ctg->kmers)[bin][cur_rbin].size() == 0) {
	pthread_mutex_unlock(&ctg->blocks[bin][cur_rbin]);
	break;
      }

      // NOTE: insert all kmers
      int kmer_idx = 0;
      for(auto i : (*ctg->kmers)[bin][cur_rbin]) {
	kmer_idx++;
#if KSET
	ctg->v[bin]->vertex_insert(i, ctg->k);
	free(i);
#elif defined(KDICT) || defined(KCOUNTER)
	ctg->v[bin]->vertex_insert(i.first, ctg->k, i.second, *ctg->merge_func);
	free(i.first);
#endif
      }

      // NOTE: clear the vector
      (*ctg->kmers)[bin][cur_rbin].clear();
      pthread_mutex_unlock(&ctg->blocks[bin][cur_rbin]);
      ctg->rbin[bin]++;
      if(ctg->rbin[bin] == ctg->work_queues) {
	ctg->rbin[bin] = 0;
      }
    }
#if defined(KDICT) || defined(KCOUNTER)
    ctg->v[bin]->burst_uc(ctg->k, *ctg->merge_func);
#else
    ctg->v[bin]->burst_uc(ctg->k);
#endif

    return NULL;
  }

#if defined(KSET) || defined(KCOUNTER)
  void parallel_kcontainer_add(const char* kmer)
#elif defined(KDICT)
  void parallel_kcontainer_add(const char* kmer, T& value)
#endif
  {
    uint8_t* pbseq = (uint8_t*) calloc(tg->bk, sizeof(uint8_t));
    int ret = serialize_kmer(kmer, k, pbseq);
    if(ret != -1) {
      free(pbseq);
      throw std::invalid_argument("Parallel add op: Could not serialize kmer, ambiguity bases present.");
    }
#if defined(KSET)
    parallel_kcontainer_add_bseq(pbseq);
#elif defined(KDICT)
    parallel_kcontainer_add_bseq(pbseq, value);
#elif defined(KCOUNTER)
    int count = 1;
    parallel_kcontainer_add_bseq(pbseq, count);
#endif
  }

#if defined(KSET) || defined(KCOUNTER)
  void parallel_kcontainer_add_seq(const char* seq, uint32_t length)
#elif defined(KDICT)
#if defined(PYTHON)
  void parallel_kcontainer_add_seq(const char* seq, uint32_t length, py::iterable& values)
#else
  template <typename Iterable>
  void parallel_kcontainer_add_seq(const char* seq, uint32_t length, Iterable& values)
#endif
#endif
  {
    //std::cout << "parallel kcontainer add seq" << std::endl;
    int size64 = tg->k / 32;
    if(tg->k % 32 > 0) {
      size64++;
    }

    int i;
    uint64_t* bseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
    uint8_t* bseq8 = (uint8_t*) bseq64;
    uint64_t* bseq64_sub = (uint64_t*) calloc(size64, sizeof(uint64_t));
    uint8_t* bseq8_sub = (uint8_t*) bseq64_sub;
    uint8_t last_index = (tg->k - 1) % 4;

    // serialize the first kmer
    int ret = serialize_kmer(seq, k, bseq8);
    int start = 0;
    while(ret != -1 && start + ret + k < length) {
      start += ret + 1;
      ret = serialize_kmer(&seq[start], k, bseq8);
    }

    for(i = 0; i < size64; i++) {
      bseq64_sub[i] = bseq64[i];
    }
#if defined(KSET)
    parallel_kcontainer_add_bseq(bseq8_sub);
#elif defined(KCOUNTER)
    int count = 1;
    parallel_kcontainer_add_bseq(bseq8_sub, count);
#elif defined(KDICT)

#if defined(PYTHON)
    auto iter = py::iter(values);
    auto value = iter->cast<T>();
#else
    auto iter = values.begin();
    auto value = *iter;
#endif
    parallel_kcontainer_add_bseq(bseq8_sub, value);

#endif

    for(uint32_t j = start + tg->k; j < length; j++) {
      bseq64[0] >>= 2;
      for(int i = 1; i < size64; i++) {
        bseq64[i - 1] |= (bseq64[i] << 62);
        bseq64[i] >>= 2;
      }

      ret = serialize_position(j, tg->bk - 1, last_index, bseq8, seq);

      if(ret != -1) {
	start = j - k + 1;
	ret = k - 1;

	while(ret != -1 && start + ret + k < length) {
	  // NOTE: can't big shift part of it since a lot of could be wrong
	  start += ret + 1;
	  memset(bseq64, 0, sizeof(uint64_t) * size64);
	  ret = serialize_kmer(&seq[start], k, bseq8);
	}

	if(ret != -1) {
	  break;
	}

	j = start + k - 1;
      }
 
      //std::cout << "inserting: " << deserialize_kmer(k, bseq8) << std::endl;
      bseq64_sub = (uint64_t*) calloc(size64, sizeof(uint64_t));
      bseq8_sub = (uint8_t*) bseq64_sub;
      for(i = 0; i < size64; i++) {
	bseq64_sub[i] = bseq64[i];
      }

#if defined(KSET)
      //std::cout << "about to add 2" << std::endl;
      parallel_kcontainer_add_bseq(bseq8_sub);
#elif defined(KCOUNTER)
      int count = 1;
      parallel_kcontainer_add_bseq(bseq8_sub, count);
#elif defined(KDICT)
      //std::cout << j << "/" << length << "\tadding val 2: " << std::endl;; //<< std::string(py::str(*iter)) << std::endl;
      #if defined(PYTHON)
      py::gil_scoped_acquire acquire;
      #endif
      std::advance(iter, 1);
      #if defined(PYTHON)
      py::gil_scoped_release release;
      #endif

      //std::cout << j << " casting" << std::endl;
      #if defined(PYTHON)
      //handle_value = *iter;
      //value = handle_value.cast<T>();
      auto value = iter->cast<T>();
      #else
      value = *iter;
      #endif

      parallel_kcontainer_add_bseq(bseq8_sub, value);
#endif
    }
    //std::cout << "done adding stuff" << std::endl;

    free(bseq64);
}

#if defined(KSET)
  void parallel_kcontainer_add_bseq(uint8_t* bseq)
#elif defined(KDICT) || defined(KCOUNTER)
  void parallel_kcontainer_add_bseq(uint8_t* bseq, T& obj)
#endif
  {
    //std::cout << "adding bseq" << std::endl;
    uint idx = (unsigned) bseq[0];

    // NOTE: determine bin
    // hard coded for 4 threads right now
    uint bin = idx >> tg->bits_to_shift;
    int cur_wbin = tg->wbin[bin];
    //std::cout << "bin: " << bin << "\t" << cur_wbin << std::endl;

    // NOTE: check if producer has the mutex
    pthread_mutex_lock(&tg->blocks[bin][cur_wbin]);

    // NOTE: add to queuue
#if defined(KSET)
    (*tg->kmers)[bin][cur_wbin].push_back(bseq);
#elif defined(KDICT) || defined(KCOUNTER)
    std::pair<uint8_t*, T> data(bseq, obj);
    (*tg->kmers)[bin][cur_wbin].push_back(data);
#endif

    // NOTE: if there are enough items in the queue, release the mutex
    if((*tg->kmers)[bin][cur_wbin].size() == tg->MAX_BIN_SIZE) {
      // NOTE: move to next thread queue
      tg->wbin[bin]++;
      if(tg->wbin[bin] == tg->work_queues) {
	tg->wbin[bin] = 0;
      }

      // NOTE: signal to thread to work
      sem_post(tg->rsignal[bin]);
    }
    pthread_mutex_unlock(&tg->blocks[bin][cur_wbin]);
  }
private:
  int k;
  
#if defined(KDICT) || defined(KCOUNTER)
  Vertex<T> v;
  ThreadInfo<T>* ti;
  ThreadGlobals<T>* tg;
#else
  ThreadInfo* ti;
  ThreadGlobals* tg;
  Vertex v;
#endif
};


