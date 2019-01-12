#ifndef __PGDATA_H_
#define __PGDATA_H_

#include <list>

struct PgData{
  //uint32_t vidx;
  std::list<uint32_t>* coords;
  uint32_t genomes;
  std::list<uint8_t>* counts;
  uint8_t size;
};

#endif // __PGDATA_H_
