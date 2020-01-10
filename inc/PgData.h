#ifndef __PGDATA_H_
#define __PGDATA_H_

#include <list>
#include <vector>
#include "helper.h"

struct PgData{
  uint32_t vidx;
  std::vector<uint32_t>* coords;
  uint32_t genomes;
  uint32_t orientation;
  //std::list<uint8_t>* counts;
  std::vector<bool>* first;
  std::vector<bool>* second;
  uint8_t size;
};

inline void free_pgdata(PgData* d) {
  delete d->coords;
  //delete d->counts;
  delete d->first;
  delete d->second;
}


/*class GraphBuild
{
  private:
    char* mseq;
    uint32_t mlen;
  public:
    GraphBuild(const char* seq) {mseq = seq;}
    ~GraphBuild(){}
}*/

#endif // __PGDATA_H_
