#include <cassert>
#include <cstdio>
#include <cstring>

#include "Kcounter.h"

int main() {
  const char* seq = "AAACTGTCTTCAAACTGTCTTT";
  const int k = 11;
  const char* path = "kcollections_cpp_kcounter_test.kc";

  Kcounter kc(k);
  kc.add_seq(seq, std::strlen(seq));
  const uint64_t expected = kc.size();

  kc.write(path);

  Kcounter loaded;
  loaded.read(path);
  std::remove(path);

  assert(loaded.get_k() == k);
  assert(loaded.size() == expected);
  assert(loaded.get("AAACTGTCTTC") >= 1);

  return 0;
}
