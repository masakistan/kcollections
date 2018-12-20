#include "Kcontainer.h"

#if defined(KSET) || defined(KCOUNTER) || defined(KCOLOR)
#if defined(KCOLOR)
std::vector<std::vector<std::vector<std::pair<uint8_t*, uint32_t>>>> kmers;
#else
std::vector<std::vector<std::vector<uint8_t*>>> kmers;
#endif

Vertex** v;
pthread_mutex_t** blocks;
sem_t** rsignal;
int bk, k, nthreads;
pthread_t* p_threads;
int* bin_ids;
int* wbin;
int* rbin;
int work_queues = 10;
int bits_to_shift;
int MAX_BIN_SIZE = 500;

void parallel_kcontainer_add_join(Kcontainer* kc) {
  for(int i = 0; i < nthreads; i++) {
    // NOTE: we post twice to finish working.
    // once to finish any kmers in the pipe
    // and a second time to pass the empty kmer list to
    // the thread to signal work is done
    sem_post(rsignal[i]);
    sem_post(rsignal[i]);
  }

  int total_vs = 0;
  uint16_t total_kmers = 0;
  for(int i = 0; i < nthreads; i++) {
    pthread_join(p_threads[i], NULL);
    //std::cout << "joined " << i << std::endl;
    total_vs += v[i]->vs_size;
    sem_close(rsignal[i]);
  }

  kc->v.vs = (Vertex*) calloc(total_vs, sizeof(Vertex));
  kc->v.vs_size = total_vs;

  // NOTE: clean up memory
  int idx = 0;
  for(int i = 0; i < nthreads; i++) {
    //std::cout << "placing vs at: " << idx << std::endl;
    std::memmove(&kc->v.vs[idx], v[i]->vs, v[i]->vs_size * sizeof(Vertex));
    kc->v.pref_pres |= v[i]->pref_pres;
    idx += v[i]->vs_size;
#if KCOLOR
    free_uc(&v[i]->uc, true);
#else
    free_uc(&v[i]->uc);
#endif
    free(v[i]->vs);
    free(v[i]);
    kmers[i].clear();
    free(blocks[i]);
  }

  kmers.clear();
  free(rsignal);
  free(v);
  free(p_threads);
  free(bin_ids);
  free(wbin);
  free(rbin);
  free(blocks);

  // NOTE: merge all vertices together
}

void parallel_kcontainer_add_init(Kcontainer* kd, int threads) {
    // NOTE: initialize variables
    k = kd->k;
    nthreads = threads;
    bits_to_shift = 8 - log2((double) threads);
    bk = calc_bk(kd->k);
    v = (Vertex**) calloc(threads, sizeof(Vertex*));
    rsignal = (sem_t**) calloc(threads, sizeof(sem_t*));
    p_threads = (pthread_t*) calloc(threads, sizeof(pthread_t));
    bin_ids = (int*) calloc(threads, sizeof(int));
    wbin = (int*) calloc(threads, sizeof(int));
    rbin = (int*) calloc(threads, sizeof(int));

    blocks = (pthread_mutex_t**) calloc(threads, sizeof(pthread_mutex_t*));
    std::string appName = "kcollections";

    // NOTE: initialize per thread variables
    // NOTE: initialize mutexes
    for(int i = 0; i < threads; i++) {
#if KCOLOR
      kmers.push_back(std::vector<std::vector<std::pair<uint8_t*, uint32_t>>>());
#else
      kmers.push_back(std::vector<std::vector<uint8_t*>>());
#endif
      blocks[i] = (pthread_mutex_t*) calloc(work_queues, sizeof(pthread_mutex_t));
      for(int j = 0; j < work_queues; j++) {
        pthread_mutex_init(&blocks[i][j], NULL);
#if KCOLOR
        kmers[i].push_back(std::vector<std::pair<uint8_t*, uint32_t>>());
#else
        kmers[i].push_back(std::vector<uint8_t*>());
#endif
      }
        v[i] = (Vertex*) calloc(1, sizeof(Vertex));
        init_vertex(v[i]);

        //ids.push_back(std::to_string(i).c_str());
        rsignal[i] = sem_open((appName + std::to_string(getpid()) + std::to_string(i)).c_str(), O_CREAT, 0600, 0);
        bin_ids[i] = i;
        // NOTE: spin up worker threads
        pthread_create(&p_threads[i], NULL, parallel_kcontainer_add_consumer, &bin_ids[i]);
    }
}

void* parallel_kcontainer_add_consumer(void* bin_ptr) {
    int bin = *((int*) bin_ptr);
    int cur_rbin;

#if KCOUNTER
    int count;
#endif

    while(true) {
        sem_wait(rsignal[bin]);
        cur_rbin = rbin[bin];
        pthread_mutex_lock(&blocks[bin][cur_rbin]);

        // NOTE: release mutex

        // NOTE: if the list is empty we're done
        if(kmers[bin][cur_rbin].size() == 0) {
          pthread_mutex_unlock(&blocks[bin][cur_rbin]);
          break;
        }

        // NOTE: insert all kmers
        for(auto i : kmers[bin][cur_rbin]) {
#if KSET
            vertex_insert(v[bin], i, k, 0);
            free(i);
#elif KCOLOR
            roaring_bitmap_t* r = roaring_bitmap_of(i.second);
            vertex_insert(v[bin], i.first, k, 0, r);
            free(i.first);
#elif KCOUNTER
            count = vertex_get_counter(v[bin], i, k, 0);
            vertex_insert(v[bin], i, k, 0, ++count);
            free(i);
#endif
        }

        // NOTE: clear the vector
        kmers[bin][cur_rbin].clear();
        pthread_mutex_unlock(&blocks[bin][cur_rbin]);
        rbin[bin]++;
        if(rbin[bin] == work_queues) {
          rbin[bin] = 0;
        }
    }
    burst_uc(v[bin], k, 0);

    return NULL;
}

#if defined(KCOLOR)
void parallel_kcontainer_add_bseq(Kcontainer* kd, uint8_t* bseq, uint32_t color)
#else
void parallel_kcontainer_add_bseq(Kcontainer* kd, uint8_t* bseq)
#endif
{
   uint idx = (unsigned) bseq[0];

  // NOTE: determine bin
  // hard coded for 4 threads right now
  uint bin = idx >> bits_to_shift;
  int cur_wbin = wbin[bin];

  // NOTE: check if producer has the mutex
  pthread_mutex_lock(&blocks[bin][cur_wbin]);

  // NOTE: add to queuue
#if defined(KCOLOR)
  kmers[bin][cur_wbin].push_back(std::make_pair(bseq, color));
#else
  kmers[bin][cur_wbin].push_back(bseq);
#endif

  // NOTE: if there are enough items in the queue, release the mutex
  if(kmers[bin][cur_wbin].size() == MAX_BIN_SIZE) {
    // NOTE: move to next thread queue
    wbin[bin]++;
    if(wbin[bin] == work_queues) {
      wbin[bin] = 0;
    }

    // NOTE: signal to thread to work
    sem_post(rsignal[bin]);
  }
  pthread_mutex_unlock(&blocks[bin][cur_wbin]);
}

#if defined(KCOLOR)
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer, uint32_t color)
#else
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer)
#endif
{
    // NOTE: start adding
  uint8_t* pbseq = ( uint8_t* ) calloc( bk, sizeof( uint8_t ) );
  serialize_kmer(kmer, kd->k, pbseq);
#if defined(KCOLOR)
  parallel_kcontainer_add_bseq(kd, pbseq, color);
#else
  parallel_kcontainer_add_bseq(kd, pbseq);
#endif
}

#if defined(KCOLOR)
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, uint32_t color)
#else
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length)
#endif
{
  int size64 = kd->k / 32;
  if(kd->k % 32 > 0) {
      size64++;
  }

  int i;
  uint64_t* bseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
  uint8_t* bseq8 = (uint8_t*) bseq64;
  uint64_t* bseq64_sub = (uint64_t*) calloc(size64, sizeof(uint64_t));
  uint8_t* bseq8_sub = (uint8_t*) bseq64_sub;

  uint bk = calc_bk(kd->k);
  int count = 0;
  uint8_t holder;
  uint8_t last_index = (kd->k - 1) % 4;

  // serialize the first kmer
  serialize_kmer(seq, kd->k, bseq8);
  for(i = 0; i < size64; i++) {
    bseq64_sub[i] = bseq64[i];
  }
#if defined(KCOLOR)
  parallel_kcontainer_add_bseq(kd, bseq8_sub, color);
#else
  parallel_kcontainer_add_bseq(kd, bseq8_sub);
#endif

  for(int j = kd->k; j < length; j++) {
    // shift all the bits over
    bseq64[0] >>= 2;
    //bseq8[0] <<= 2;
    //for(int i = 1; i < size64; i++) {
    for(int i = 1; i < size64; i++) {
        bseq64[i - 1] |= (bseq64[i] << 62);
        //bseq64[i - 1] |= (bseq64[i] >> 62);
        bseq64[i] >>= 2;
        //bseq64[i] <<= 2;
    }

    serialize_position(j, bk - 1, last_index, bseq8, seq);
    //std::cout << "inserting: " << deserialize_kmer(k, calc_bk(k), bseq8) << std::endl;

    bseq64_sub = (uint64_t*) calloc(size64, sizeof(uint64_t));
    bseq8_sub = (uint8_t*) bseq64_sub;
    for(i = 0; i < size64; i++) {
      bseq64_sub[i] = bseq64[i];
    }

#if defined(KCOLOR)
    parallel_kcontainer_add_bseq(kd, bseq8_sub, color);
#else
    parallel_kcontainer_add_bseq(kd, bseq8_sub);
#endif
  }

  free(bseq64);
}

#if defined(KCOLOR)
void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, uint32_t color)
#else
void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length)
#endif
{
    int size64 = kd->k / 32;
    if(kd->k % 32 > 0) {
        size64++;
    }

    uint64_t* bseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
    uint8_t* bseq8 = (uint8_t*) bseq64;

    uint bk = calc_bk(kd->k);
    uint8_t holder;
    uint8_t last_index = (kd->k - 1) % 4;

    // serialize the first kmer
    serialize_kmer(seq, kd->k, bseq8);

#if KSET
    vertex_insert(&(kd->v), bseq8, kd->k, 0);
#elif KCOLOR
    roaring_bitmap_t* r = roaring_bitmap_of(color);
    vertex_insert(&(kd->v), bseq8, kd->k, 0, r);
#elif KCOUNTER
    // get current count of bseq
    int count = vertex_get_counter(&(kd->v), bseq8, kd->k, 0);
    vertex_insert(&(kd->v), bseq8, kd->k, 0, ++count);
#endif

    std::cout << strlen(seq) << std::endl;
    for(uint32_t j = kd->k; j < length; j++) {
        //std::cout << j << "\t" << seq[j] << std::endl;
        // shift all the bits over
        //bseq8[0] <<= 2;
        bseq64[0] >>= 2;
        //std::cout << "shifting\t" << deserialize_kmer(kd->k, calc_bk(kd->k), bseq8) << std::endl;
        //for(int i = 1; i < bk; i++) {
        for(int i = 1; i < size64; i++) {
            bseq64[i - 1] |= (bseq64[i] << 62);
            bseq64[i] >>= 2;
            //std::cout << "shifting\t" << deserialize_kmer(kd->k, calc_bk(kd->k), bseq8) << std::endl;
        }

        serialize_position(j, bk - 1, last_index, bseq8, seq);
#if KSET
        vertex_insert(&(kd->v), bseq8, kd->k, 0);
#elif KCOLOR
        roaring_bitmap_t* r = roaring_bitmap_of(color);
        vertex_insert(&(kd->v), bseq8, kd->k, 0, r);
#elif KCOUNTER
        // get current count of bseq
        count = vertex_get_counter(&(kd->v), bseq8, kd->k, 0);
        vertex_insert(&(kd->v), bseq8, kd->k, 0, ++count);
#endif
    }

    free(bseq64);

}
#endif
