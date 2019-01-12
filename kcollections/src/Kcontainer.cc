#include "Kcontainer.h"

#if defined(KSET) || defined(KCOUNTER)
std::vector<std::vector<std::vector<std::tuple<uint16_t, uint32_t, uint8_t*>>>> kmers;
Vertex** v;
sem_t* signal_b;
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
    free_uc(&v[i]->uc);
    free(v[i]->vs);
    free(v[i]);
    kmers[i].clear();
  }

  kmers.clear();
  free(signal_b);
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
    signal_b = (sem_t*) calloc(threads, sizeof(sem_t));
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
      kmers.push_back(std::vector<std::vector<std::tuple<uint16_t, uint32_t, uint8_t*>>>());
      blocks[i] = (pthread_mutex_t*) calloc(work_queues, sizeof(pthread_mutex_t));
      for(int j = 0; j < work_queues; j++) {
        pthread_mutex_init(&blocks[i][j], NULL);
        kmers[i].push_back(std::vector<std::tuple<uint16_t, uint32_t, uint8_t*>>());
      }
        v[i] = (Vertex*) calloc(1, sizeof(Vertex));
        init_vertex(v[i]);

        //sem_init(&signal_b[i], 0, 0);
        //ids.push_back(std::to_string(i).c_str());
        rsignal[i] = sem_open((appName + std::to_string(getpid()) + std::to_string(i)).c_str(), O_CREAT, 0600, 0);
        //rsignal[i] = &signal_b[i];
        bin_ids[i] = i;
        // NOTE: spin up worker threads
        pthread_create(&p_threads[i], NULL, parallel_kcontainer_add_consumer, &bin_ids[i]);
    }
}

void* parallel_kcontainer_add_consumer(void* bin_ptr) {
    int bin = *((int*) bin_ptr);
    int cur_rbin;
    uint8_t* bseq;

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
            bseq = std::get<2>(i);
            vertex_insert(v[bin], bseq, k, 0, &i, false);
            free(bseq);
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
    filter_vertices(v[bin], calc_bk(k));

    return NULL;
}

void parallel_kcontainer_add_bseq(Kcontainer* kd, uint8_t* bseq, uint16_t gidx, uint32_t pos) {
   uint idx = (unsigned) bseq[0];

  // NOTE: determine bin
  // hard coded for 4 threads right now
  uint bin = idx >> bits_to_shift;
  int cur_wbin = wbin[bin];

  // NOTE: check if producer has the mutex
  pthread_mutex_lock(&blocks[bin][cur_wbin]);

  // NOTE: add to queuue
  //std::cout << gidx << "\t" << pos << std::endl;
  kmers[bin][cur_wbin].push_back(std::make_tuple(gidx, pos, bseq));

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

void parallel_kcontainer_add(Kcontainer* kd, const char* kmer, uint16_t gidx, uint32_t pos) {
    // NOTE: start adding
  uint8_t* pbseq = ( uint8_t* ) calloc( bk, sizeof( uint8_t ) );
  serialize_kmer(kmer, kd->k, pbseq);
  parallel_kcontainer_add_bseq(kd, pbseq, gidx, pos);
}

void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, uint16_t gidx, uint32_t offset) {
  int size64 = kd->k / 32;
  if(kd->k % 32 > 0) {
      size64++;
  }

  int i;
  uint64_t* cbseq64;

  uint64_t* fbseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
  uint8_t* fbseq8 = (uint8_t*) fbseq64;

  uint64_t* rbseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
  uint8_t* rbseq8 = (uint8_t*) rbseq64;

  uint64_t* bseq64_sub = (uint64_t*) calloc(size64, sizeof(uint64_t));
  uint8_t* bseq8_sub = (uint8_t*) bseq64_sub;

  uint bk = calc_bk(kd->k);
  int count = 0;
  uint8_t holder;
  uint8_t last_index = (kd->k - 1) % 4;
  uint8_t reverse_clear_bits_to_shift = last_index * 2;

  // serialize the first kmer
  serialize_kmer(seq, kd->k, fbseq8);
  serialize_kmer_rev(seq, kd->k, rbseq8);

  //std::cout << "for: " << deserialize_kmer(kd->k, bk, fbseq8) << std::endl;
  //std::cout << "rev: " << deserialize_kmer(kd->k, bk, rbseq8) << std::endl;

  cbseq64 = fbseq64;
  if(memcmp(fbseq64, rbseq64, sizeof(uint64_t) * size64) > 0) {
    cbseq64 = rbseq64;
  }

  for(i = 0; i < size64; i++) {
    bseq64_sub[i] = cbseq64[i];
  }
  parallel_kcontainer_add_bseq(kd, bseq8_sub, gidx, offset);

  for(int j = kd->k; j < length; j++) {
    offset++;

    // shift all the bits over
    fbseq64[0] >>= 2;

    rbseq8[bk - 1] &= ~(3 << reverse_clear_bits_to_shift);
    rbseq64[size64 - 1] <<= 2;


    for(int i = 1; i < size64; i++) {
        fbseq64[i - 1] |= (fbseq64[i] << 62);
        fbseq64[i] >>= 2;

        rbseq64[size64 - i] |= (rbseq64[size64 - i - 1] >> 62);
        rbseq64[size64 - i - 1] <<= 2;
    }

    serialize_position(j, bk - 1, last_index, fbseq8, seq);
    serialize_position_comp(j, 0, 0, rbseq8, seq);
    //std::cout << "inserting: " << deserialize_kmer(k, calc_bk(k), bseq8) << std::endl;
    //std::cout << "for: " << deserialize_kmer(bk * 4, bk, fbseq8) << std::endl;
    //std::cout << "rev: " << deserialize_kmer(bk * 4, bk, rbseq8) << std::endl;

    //std::cout << "******************" << std::endl;
    //std::cout << "for: " << deserialize_kmer(kd->k, bk, fbseq8) << std::endl;
    //std::cout << "rev: " << deserialize_kmer(kd->k, bk, rbseq8) << std::endl;

    cbseq64 = fbseq64;
    if(memcmp(fbseq64, rbseq64, sizeof(uint64_t) * size64) > 0) {
      cbseq64 = rbseq64;
    }

    bseq64_sub = (uint64_t*) calloc(size64, sizeof(uint64_t));
    bseq8_sub = (uint8_t*) bseq64_sub;
    for(i = 0; i < size64; i++) {
      bseq64_sub[i] = cbseq64[i];
    }

    parallel_kcontainer_add_bseq(kd, bseq8_sub, gidx, offset);
  }
}

void filter_vertices(Vertex* v, int bk) {
  // NOTE: iterate over suffixes
  std::list<uint8_t>::iterator countIter;
  std::list<uint32_t>::iterator coordIter;

  for(int i = 0; i < v->uc.size; i++) {
    PgData* kdata = &v->uc.data[i];
    int gidx = 0;
    countIter = kdata->counts->begin();
    coordIter = kdata->coords->begin();

    if(kdata->counts->size() == 1) {
      int suffix_idx = bk * i;
      v->uc.size -= 1;
      memmove(
              &v->uc.suffixes[suffix_idx],
              &v->uc.suffixes[suffix_idx + bk],
              (v->uc.size - i) * sizeof(uint8_t) * bk
              );
      v->uc.suffixes = (uint8_t*) realloc(v->uc.suffixes, sizeof(uint8_t) * v->uc.size * bk);

      free_pgdata(kdata);
      memmove(
              &v->uc.data[i],
              &v->uc.data[i + 1],
              (v->uc.size - i) * sizeof(PgData)
              );
      v->uc.data = (PgData*) realloc(v->uc.data, sizeof(PgData) * v->uc.size);
      i -= 1;

      continue;
    }

    while(countIter != kdata->counts->end()) {
      // NOTE: get gidx
      gidx = next_set_bit(&kdata->genomes, gidx, 32);

      if(*countIter == 2) {
        countIter = kdata->counts->erase(countIter);
        coordIter = kdata->coords->erase(coordIter);
        kdata->genomes &= ~(0x1 << gidx);
      } else {
        countIter++;
        coordIter++;
      }
      gidx++;
    }
  }

  for(int i = 0; i < v->vs_size; i++) {
    filter_vertices(&v->vs[i], bk - 1);
  }
}

#endif
