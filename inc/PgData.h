#ifndef __PGDATA_H_
#define __PGDATA_H_

struct PgData{
  //uint32_t vidx;
  std::vector<uint32_t>* coords;
  uint32_t genomes;
  std::vector<uint8_t>* counts;
  uint8_t size;
};

#endif // __PGDATA_H_
