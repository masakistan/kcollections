#include "Kcontainer.h"

#if KSET
std::vector<std::vector<std::vector<uint8_t*>>> kmers;
Vertex** v;
sem_t* signal_b;
pthread_mutex_t** blocks;
sem_t** rsignal;
int bk, k, nthreads;
pthread_t* p_threads;
int* bin_ids;
int* wbin;
int* rbin;
int work_queues = 5;

void parallel_kcontainer_insert_join(Kcontainer* kc) {
  for(int i = 0; i < nthreads; i++) {
    //std::cout << "last post for bin " << i << std::endl;
    // NOTE: we post twice to finish working.
    // once to finish any kmers in the pipe
    // and a second time to pass the empty kmer list to
    // the thread to signal work is done
    sem_post(rsignal[i]);
    sem_post(rsignal[i]);
  }

  int total_ccs = 0;
  uint16_t total_kmers = 0;
  for(int i = 0; i < nthreads; i++) {
    pthread_join(p_threads[i], NULL);
    total_ccs += v[i]->cc_size;
  }

  kc->v.cc = (CC*) calloc(total_ccs, sizeof(CC));
  kc->v.cc_size = total_ccs;

  //std::cout << "total compressed containers: " << total_ccs << std::endl;

  // NOTE: clean up memory
  int idx = 0;
  for(int i = 0; i < nthreads; i++) {
    //std::cout << "placing CCs at: " << idx << std::endl;
    std::memmove(&kc->v.cc[idx], v[i]->cc, v[i]->cc_size * sizeof(CC));
    idx += v[i]->cc_size;
    free_uc(&v[i]->uc);
    free(v[i]->cc);
    free(v[i]);
  }

  free(signal_b);
  free(rsignal);
  free(v);
  free(p_threads);
  free(bin_ids);

  // NOTE: merge all vertices together
}

void parallel_kcontainer_insert_init(Kcontainer* kd, int threads) {
    // NOTE: initialize variables
    k = kd->k;
    nthreads = threads;
    bk = calc_bk(kd->k);
    v = (Vertex**) calloc(threads, sizeof(Vertex*));
    signal_b = (sem_t*) calloc(threads, sizeof(sem_t));
    rsignal = (sem_t**) calloc(threads, sizeof(sem_t*));
    p_threads = (pthread_t*) calloc(threads, sizeof(pthread_t));
    bin_ids = (int*) calloc(threads, sizeof(int));
    wbin = (int*) calloc(threads, sizeof(int));
    rbin = (int*) calloc(threads, sizeof(int));

    blocks = (pthread_mutex_t**) calloc(threads, sizeof(pthread_mutex_t*));

    // NOTE: initialize per thread variables
    // NOTE: initialize mutexes
    for(int i = 0; i < threads; i++) {
      kmers.push_back(std::vector<std::vector<uint8_t*>>());
      blocks[i] = (pthread_mutex_t*) calloc(work_queues, sizeof(pthread_mutex_t));
      for(int j = 0; j < work_queues; j++) {
        pthread_mutex_init(&blocks[i][j], NULL);
        kmers[i].push_back(std::vector<uint8_t*>());
      }
        v[i] = (Vertex*) calloc(1, sizeof(Vertex));
        init_vertex(v[i]);

        sem_init(&signal_b[i], 0, 0);
        rsignal[i] = &signal_b[i];
        bin_ids[i] = i;
        // NOTE: spin up worker threads
        pthread_create(&p_threads[i], NULL, parallel_insert_consumer, &bin_ids[i]);
    }
    //sleep(1);
}

void* parallel_insert_consumer(void* bin_ptr) {
    int bin = *((int*) bin_ptr);
    int cur_rbin;
    std::cout << "spinning up thread:\t" << bin << std::endl;
    while(true) {
        sem_wait(rsignal[bin]);
        cur_rbin = rbin[bin];
        pthread_mutex_lock(&blocks[bin][cur_rbin]);
        //std::cout << "\tworking on " << bin << ":" << cur_rbin << "\t" << kmers[bin][cur_rbin].size() << std::endl;

        // NOTE: release mutex

        // NOTE: if the list is empty we're done
        if(kmers[bin][cur_rbin].size() == 0) {
          //std::cout << "thread " << bin << " is done!" << std::endl;
          pthread_mutex_unlock(&blocks[bin][cur_rbin]);
          break;
        }

        // NOTE: insert all kmers
        for(auto i : kmers[bin][cur_rbin]) {
            vertex_insert(v[bin], i, k, 0);
            free(i);
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
}

void parallel_kcontainer_insert(Kcontainer* kd, const char* kmer) {
    // NOTE: start adding
  uint8_t* pbseq = ( uint8_t* ) calloc( bk, sizeof( uint8_t ) );
  serialize_kmer(kmer, kd->k, pbseq);
  uint idx = (unsigned) pbseq[0];

  // NOTE: determine bin
  // hard coded for 4 threads right now
  uint bin = idx >> 6;
  int cur_wbin = wbin[bin];

  // NOTE: check if producer has the mutex
  //std::cout << "trying to lock:\t" << bin << std::endl;
  pthread_mutex_lock(&blocks[bin][cur_wbin]);

  // NOTE: add to queuue
  kmers[bin][cur_wbin].push_back(pbseq);
  //std::cout << "added to " << bin << "(" << wkmers[bin]->size() << ") from index:" << idx << std::endl;


  // NOTE: if there are enough items in the queue, release the mutex
  if(kmers[bin][cur_wbin].size() == 5000) {
    //std::cout << "bin " << bin << ":" << cur_wbin << " is done" << std::endl;
    // NOTE: move to next thread queue
    wbin[bin]++;
    if(wbin[bin] == work_queues) {
      wbin[bin] = 0;
    }

    // NOTE: signal to thread to work
    sem_post(rsignal[bin]);
    //sleep(1);
  }
  pthread_mutex_unlock(&blocks[bin][cur_wbin]);
}

void release_all_locks() {

}
#endif
