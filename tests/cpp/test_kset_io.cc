#include <cassert>
#include <cstdio>
#include <cstring>

#include "Kset.h"

int main() {
  const char* seq = "AAACTGTCTTCAAACTGTCTTT";
  const int k = 11;
  const char* path = "kcollections_cpp_test.kc";

  Kset ks(k);
  ks.add_seq(seq, std::strlen(seq));
  const uint64_t expected = ks.size();

  ks.write(path);

  Kset loaded;
  loaded.read(path);
  std::remove(path);

  assert(loaded.get_k() == k);
  assert(loaded.size() == expected);
  assert(loaded.contains("AAACTGTCTTC"));
  assert(!loaded.contains("GGGGGGGGGGG"));

  return 0;
}
