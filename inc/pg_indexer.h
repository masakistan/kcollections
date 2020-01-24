#include "Kdict.h"


struct PgData{
  uint64_t vidx;
  std::list<uint32_t>* coords;
  uint64_t genomes;
  uint64_t orientation;
  std::list<uint8_t>* counts;
  uint8_t size;
};


PgData& merge_func(PgData& o, PgData& n) {
  return o;
}


class PgIndex {
 private:
  Kdict<PgData>* kd;
  int k;
  
 public:
  PgIndex(int k) {
    kd = new Kdict<PgData>(k);
    kd->set_merge_func(merge_func);
    this.k = k;
  }

  ~PgIndex() {
    delete kd;
  }
};
