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
#include <jemalloc/jemalloc.h>
#include <math.h>
#include <functional>

#if defined(PYTHON)
namespace py = pybind11;
#endif

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif

struct ThreadInfo{
  int thread_id;
#if defined(KDICT) || defined(KCOUNTER)
  std::function<T(T, T)>* merge_func;
#endif
};

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
  Vertex<T>** v;
#else
  Vertex** v;
#endif

  sem_t* signal_b;
  pthread_mutex_t** blocks;
  sem_t** rsignal;
  int bk, k, nthreads;
  pthread_t* p_threads;
  //int* bin_ids;
#if defined(KDICT) || defined(KCOUNTER)
  ThreadInfo<T>* bin_ids;
#else
  ThreadInfo* bin_ids;
#endif
  int* wbin;
  int* rbin;
  int work_queues;
  int bits_to_shift;
  size_t MAX_BIN_SIZE;
};


#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif
class Kcontainer {
private:
  int k;

#if defined(KDICT) || defined(KCOUNTER)
  Vertex<T>* v;
#else
  Vertex* v;
#endif

protected:
#if defined(KDICT) || defined(KCOUNTER)
  static ThreadGlobals<T>* tg;
#elif defined(KSET)
  static ThreadGlobals* tg;
#endif

public:
  Kcontainer(const int k) : k(k) {
#if defined(KDICT) || defined(KCOUNTER)
    v = new Vertex<T>();
#else
    v = new Vertex();
#endif
  }

  ~Kcontainer() {
    delete v;
  }

#if defined(KDICT) || defined(KCOUNTER)
  Vertex<T>* get_v() { return v; }
#else
  Vertex* get_v() { return v; }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  std::string kcontainer_get_child_suffix(Vertex<T>* v, int idx) {
#else
  std::string kcontainer_get_child_suffix(Vertex* v, int idx) {
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
    serialize_kmer(kmer, k, bseq);
    bool res = v->vertex_contains(bseq, k);
    free(bseq);
    return res;
  }

  void kcontainer_remove(const char* kmer)
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    serialize_kmer(kmer, k, bseq);
    v->vertex_remove(bseq, k);
    free(bseq);
  }

#if defined(KDICT) || defined(KCOUNTER)
  void kcontainer_add(const char* kmer, T obj, std::function<T(T, T)>& merge_func)
#elif KSET
  void kcontainer_add(const char* kmer)
#endif
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    serialize_kmer(kmer, k, bseq);
#if defined(KDICT) || defined(KCOUNTER)
    v->vertex_insert(bseq, k, obj, merge_func);
#elif KSET
    v->vertex_insert(bseq, k);
#endif
    free(bseq);
  }

#if defined KDICT || defined KCOUNTER
  T kcontainer_get(char* kmer)
  {
    uint8_t* bseq = (uint8_t*) calloc(k, sizeof(uint8_t));
    serialize_kmer(kmer, k, bseq);
#if defined(KDICT) || defined(KCOUNTER)
    T res = v->vertex_get(bseq, k);
#endif
    free(bseq);
    return res;
  }
#endif

  uint64_t kcontainer_size()
  {
    return v->get_vertex_size();
  }

#if defined(KDICT)
#if defined(PYTHON)
  void kcontainer_add_seq(const char* seq, uint32_t length, py::iterable& values, std::function<T(T, T)> &f)
#else
  template <typename Iterable>
  void kcontainer_add_seq(const char* seq, uint32_t length, Iterable& values, std::function<T(T, T)> &f)
#endif
#elif defined(KCOUNTER)
  void kcontainer_add_seq(const char* seq, uint32_t length, std::function<T(T, T)> &f)
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
    serialize_kmer(seq, k, bseq8);

#if KSET
    v->vertex_insert(bseq8, k);
#elif defined(KCOUNTER)
    v->vertex_insert(bseq8, k, 1, f);
#elif defined(KDICT)
#if defined(PYTHON)
    auto iter = py::iter(values);
    v->vertex_insert(bseq8, k, (*iter).cast<T>(), f);
#else
    auto iter = values.begin();
    v->vertex_insert(bseq8, k, *iter, f);
#endif
#endif

    //std::cout << strlen(seq) << std::endl;
    for(uint32_t j = k; j < length; j++) {
      //std::cout << j << "\t" << seq[j] << std::endl;
      // shift all the bits over
      //bseq8[0] <<= 2;
      bseq64[0] >>= 2;
      //std::cout << "shifting\t" << deserialize_kmer(k, bseq8) << std::endl;
      //for(int i = 1; i < bk; i++) {
      for(int i = 1; i < size64; i++) {
	bseq64[i - 1] |= (bseq64[i] << 62);
	bseq64[i] >>= 2;
	//std::cout << "shifting\t" << deserialize_kmer(k, bseq8) << std::endl;
      }

      serialize_position(j, bk - 1, last_index, bseq8, seq);
#if KSET
      v->vertex_insert(bseq8, k);
#elif defined(KCOUNTER)
      v->vertex_insert(bseq8, k, 1, f);
#elif defined(KDICT)
#if defined(PYTHON)
      std::advance(iter, 1);
      v->vertex_insert(bseq8, k, (*iter).cast<T>(), f);
#else
      v->vertex_insert(bseq8, k, *iter, f);
#endif
#endif
    }

    free(bseq64);
  }

#if defined(KSET)
  void parallel_kcontainer_add_init(int threads)
#elif defined(KDICT) || defined(KCOUNTER)
  void parallel_kcontainer_add_init(int threads, const std::function<T(T, T)> &f)
#endif
  {
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

#if defined(KDICT) || defined(KCOUNTER)
    tg->bin_ids = (ThreadInfo<T>*) calloc(threads, sizeof(ThreadInfo<T>));
#else
    tg->bin_ids = (ThreadInfo*) calloc(threads, sizeof(ThreadInfo));
#endif

    tg->wbin = (int*) calloc(threads, sizeof(int));
    tg->rbin = (int*) calloc(threads, sizeof(int));

    tg->blocks = (pthread_mutex_t**) calloc(threads, sizeof(pthread_mutex_t*));
    std::string appName = "kcollections";

#if defined(KDICT) || defined(KCOUNTER)
    tg->kmers = (std::vector<std::vector<std::vector<std::pair<uint8_t*, T>>>>*) calloc(threads, sizeof(std::vector<std::vector<std::vector<uint8_t*, T>>>*));
#elif defined(KSET)
    tg->kmers = (std::vector<std::vector<std::vector<uint8_t*>>>*) calloc(threads, sizeof(std::vector<std::vector<std::vector<uint8_t*>>>*));
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
      tg->bin_ids[i].thread_id = i;
#if defined(KDICT) || defined(KCOUNTER)
      tg->bin_ids[i].merge_func = new std::function<T(T, T)>(f);
#endif

      // NOTE: spin up worker threads
      pthread_create(&tg->p_threads[i], NULL, parallel_kcontainer_add_consumer, &tg->bin_ids[i]);
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
    v->set_vs((Vertex<T>**) calloc(total_vs, sizeof(Vertex<T>*)));
#else
    v->set_vs((Vertex**) calloc(total_vs, sizeof(Vertex*)));
#endif

    v->set_vs_size(total_vs);

    // NOTE: clean up memory
    int idx = 0;
    for(int i = 0; i < tg->nthreads; i++) {
      //std::cout << "thread: " << i << "\t" << v[i]->vs_size << std::endl;

#if defined(KDICT) || defined(KCOUNTER)
      Vertex<T>** v_vs = tg->v[i]->get_vs();
#else
      Vertex** v_vs = tg->v[i]->get_vs();
#endif

      if(v_vs != NULL) {
#if defined(KDICT) || defined(KCOUNTER)
	std::memmove(&v->get_vs()[idx], tg->v[i]->get_vs(), tg->v[i]->get_vs_size() * sizeof(Vertex<T>*));
#else
	std::memmove(&v->get_vs()[idx], tg->v[i]->get_vs(), tg->v[i]->get_vs_size() * sizeof(Vertex*));
#endif

	v->set_pref_pres(v->get_pref_pres() | tg->v[i]->get_pref_pres());
	idx += tg->v[i]->get_vs_size();
	free(v_vs);
#if defined(KDICT) || defined(KCOUNTER)
	tg->v[i]->set_vs((Vertex<T>**) NULL);
#else
	tg->v[i]->set_vs((Vertex**) NULL);
#endif

      }

      delete tg->v[i];
      free(tg->blocks[i]);
      (*tg->kmers)[i].clear();
#if defined(KDICT) || defined(KCOUNTER)
      delete tg->bin_ids[i].merge_func;
#endif
    }

    //std::cout << "cleaning" << std::endl;

    free(tg->v);
    free(tg->signal_b);
    free(tg->rsignal);
    free(tg->p_threads);
    free(tg->bin_ids);
    free(tg->wbin);
    free(tg->rbin);
    free(tg->blocks);

    tg->kmers->clear();

    // NOTE: merge all vertices together
  }

  static void* parallel_kcontainer_add_consumer(void* bin_ptr) {
    //int bin = *((int*) bin_ptr);
#if defined(KDICT) || defined(KCOUNTER)
    ThreadInfo<T>* ti = (ThreadInfo<T>*) bin_ptr;
#else
    ThreadInfo* ti = (ThreadInfo*) bin_ptr;
#endif

    int bin = ti->thread_id;
    int cur_rbin;

    while(true) {
      sem_wait(tg->rsignal[bin]);
      cur_rbin = tg->rbin[bin];
      pthread_mutex_lock(&tg->blocks[bin][cur_rbin]);

      // NOTE: release mutex

      // NOTE: if the list is empty we're done
      if((*tg->kmers)[bin][cur_rbin].size() == 0) {
	pthread_mutex_unlock(&tg->blocks[bin][cur_rbin]);
	break;
      }

      // NOTE: insert all kmers
      int kmer_idx = 0;
      for(auto i : (*tg->kmers)[bin][cur_rbin]) {
	kmer_idx++;
#if KSET
	tg->v[bin]->vertex_insert(i, tg->k);
	free(i);
#elif defined(KDICT) || defined(KCOUNTER)
	tg->v[bin]->vertex_insert(i.first, tg->k, i.second, *ti->merge_func);
	free(i.first);
#endif
      }

      // NOTE: clear the vector
      (*tg->kmers)[bin][cur_rbin].clear();
      pthread_mutex_unlock(&tg->blocks[bin][cur_rbin]);
      tg->rbin[bin]++;
      if(tg->rbin[bin] == tg->work_queues) {
	tg->rbin[bin] = 0;
      }
    }
#if defined(KDICT) || defined(KCOUNTER)
    tg->v[bin]->burst_uc(tg->k, *ti->merge_func);
#else
    tg->v[bin]->burst_uc(tg->k);
#endif
    //delete v[bin]->get_uc();

    return NULL;
  }

#if defined(KSET) || defined(KCOUNTER)
  void parallel_kcontainer_add(const char* kmer)
#elif defined(KDICT)
  void parallel_kcontainer_add(const char* kmer, T& value)
#endif
  {
    uint8_t* pbseq = (uint8_t*) calloc(tg->bk, sizeof(uint8_t));
    serialize_kmer(kmer, k, pbseq);
#if defined(KSET)
    parallel_kcontainer_add_bseq(pbseq);
#elif defined(KDICT)
    parallel_kcontainer_add_bseq(pbseq, value);
#elif defined(KCOUNTER)
    parallel_kcontainer_add_bseq(pbseq, 1);
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
    serialize_kmer(seq, tg->k, bseq8);
    for(i = 0; i < size64; i++) {
      bseq64_sub[i] = bseq64[i];
    }
#if defined(KSET)
    parallel_kcontainer_add_bseq(bseq8_sub);
#elif defined(KCOUNTER)
    parallel_kcontainer_add_bseq(bseq8_sub, 1);
#elif defined(KDICT)

#if defined(PYTHON)
    auto iter = py::iter(values);
    //std::cout << "casting" << std::endl;
    parallel_kcontainer_add_bseq(bseq8_sub, (*iter).cast<T>());
#else
    auto iter = values.begin();
    parallel_kcontainer_add_bseq(bseq8_sub, *iter);
#endif
    //std::cout << "after cast insert" << std::endl;
#endif

    for(uint32_t j = tg->k; j < length; j++) {
      bseq64[0] >>= 2;
      for(int i = 1; i < size64; i++) {
        bseq64[i - 1] |= (bseq64[i] << 62);
        bseq64[i] >>= 2;
      }

      serialize_position(j, tg->bk - 1, last_index, bseq8, seq);
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
      parallel_kcontainer_add_bseq(bseq8_sub, 1);
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
      parallel_kcontainer_add_bseq(bseq8_sub, iter->template cast<T>());
      //std::cout << "after cast insert" << std::endl;
      //#else
      //std::advance(iter, 1);
      //parallel_kcontainer_add_bseq(bseq8_sub, *iter);
      //#endif
#endif
    }
    //std::cout << "done adding stuff" << std::endl;

    free(bseq64);
}

#if defined(KSET)
  void parallel_kcontainer_add_bseq(uint8_t* bseq)
#elif defined(KDICT) || defined(KCOUNTER)
  void parallel_kcontainer_add_bseq(uint8_t* bseq, T obj)
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
    //std::cout << "creating pair" << std::endl;
    std::pair<uint8_t*, T> data(bseq, obj);
    //std::cout << "pair created" << std::endl;
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
};

#if defined(KDICT) || defined(KCOUNTER)
  //template<class T>
  //ThreadGlobals<T>* Kcontainer<T>::tg = (ThreadGlobals<T>*) calloc(1, sizeof(ThreadGlobals<T>));
#else
  //ThreadGlobals* Kcontainer::tg = (ThreadGlobals*) calloc(1, sizeof(ThreadGlobals));
#endif
