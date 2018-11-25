#include "Kcontainer.h"

#if KSET
std::vector<uint8_t*>** kmers_a;
std::vector<uint8_t*>** kmers_b;
std::vector<uint8_t*>** wkmers;
std::vector<uint8_t*>** rkmers;
Vertex** v;
sem_t* signal_b;
pthread_mutex_t* locks;
sem_t** rsignal;
int bk, k, nthreads;
pthread_t* p_threads;
int* bin_ids;

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

  free(kmers_a);
  free(kmers_b);
  free(wkmers);
  free(rkmers);
  free(locks);
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
    kmers_a = (std::vector<uint8_t*>**) calloc(threads, sizeof(std::vector<uint8_t*>*));
    kmers_b = (std::vector<uint8_t*>**) calloc(threads, sizeof(std::vector<uint8_t*>*));
    wkmers = (std::vector<uint8_t*>**) calloc(threads, sizeof(std::vector<uint8_t*>*));
    rkmers = (std::vector<uint8_t*>**) calloc(threads, sizeof(std::vector<uint8_t*>*));
    locks = (pthread_mutex_t*) calloc(threads, sizeof(pthread_mutex_t));
    signal_b = (sem_t*) calloc(threads, sizeof(sem_t));
    rsignal = (sem_t**) calloc(threads, sizeof(sem_t*));
    p_threads = (pthread_t*) calloc(threads, sizeof(pthread_t));
    bin_ids = (int*) calloc(threads, sizeof(int));

    // NOTE: initialize per thread variables
    // NOTE: initialize mutexes
    for(int i = 0; i < threads; i++) {
        v[i] = (Vertex*) calloc(1, sizeof(Vertex));
        init_vertex(v[i]);
        kmers_a[i] = new std::vector<uint8_t*>();
        kmers_b[i] = new std::vector<uint8_t*>();
        wkmers[i] = kmers_a[i];
        rkmers[i] = kmers_b[i];

        pthread_mutex_init(&locks[i], NULL);

        sem_init(&signal_b[i], 0, 0);
        rsignal[i] = &signal_b[i];
        bin_ids[i] = i;
        // NOTE: spin up worker threads
        pthread_create(&p_threads[i], NULL, parallel_insert_consumer, &bin_ids[i]);
    }
}


void* parallel_insert_consumer(void* bin_ptr) {
    int bin = *((int*) bin_ptr);
    //std::cout << "spinning up thread:\t" << bin << std::endl;
    while(true) {
        sem_wait(rsignal[bin]);
        pthread_mutex_lock(&locks[bin]);
        // NOTE: swap lists
        std::vector<uint8_t*>* t_list = wkmers[bin];
        wkmers[bin] = rkmers[bin];
        rkmers[bin] = t_list;
        //std::cout << "\treading from: " << bin << "\t" << rkmers[bin] << std::endl;
        //std::cout << "\tstarting work on bin " << bin << "\t(" << rkmers[bin]->size() << " items)" << std::endl;

        // NOTE: release mutex
        pthread_mutex_unlock(&locks[bin]);

        // NOTE: if the list is empty we're done
        if(rkmers[bin]->size() == 0) {
          //std::cout << "bin " << bin << " is done!" << std::endl;
            break;
        }

        // NOTE: insert all kmers
        for(auto i : *rkmers[bin]) {
            vertex_insert(v[bin], i, k, 0);
            free(i);
        }

        // NOTE: clear the vector
        rkmers[bin]->clear();
    }

    delete kmers_a[bin];
    delete kmers_b[bin];
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

    // NOTE: check if producer has the mutex
    //std::cout << "trying to lock:\t" << bin << std::endl;
    pthread_mutex_lock(&locks[bin]);

    // NOTE: add to queuue
    wkmers[bin]->push_back(pbseq);
    //std::cout << "added to " << bin << "(" << wkmers[bin]->size() << ") from index:" << idx << std::endl;

    // NOTE: unlock the mutex

    // NOTE: if there are enough items in the queue, release the mutex
    if(wkmers[bin]->size() == 5000) {
        int value;
        sem_getvalue(rsignal[bin], &value);
        if(value != 1) {
          //std::cout << "writing to: " << bin << "\t" << wkmers[bin] << std::endl;
            sem_post(rsignal[bin]);
            //std::cout << "posting! " << bin << "\t" << value  << "\t" << wkmers[bin]->size() << std::endl;
        }
    }
    pthread_mutex_unlock(&locks[bin]);
}

void release_all_locks() {

}
#endif
