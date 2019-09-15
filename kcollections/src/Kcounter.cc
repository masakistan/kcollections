#include "Kcounter.h"


Kcounter::Kcounter(int k)
{
  kc = new Kcontainer<int>(k);
  m_k = k;
}

Kcounter::~Kcounter()
{
  delete kc;
}

void Kcounter::clear()
{
  delete kc;
  kc = new Kcontainer<int>(m_k);
}

void Kcounter::insert(char* kmer, count_dtype count)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  kc->kcontainer_add(kmer, count, merge_func );
}

count_dtype Kcounter::get(char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  return kc->kcontainer_get(kmer);
}

bool Kcounter::contains(char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  return kc->kcontainer_contains(kmer);
}

uint64_t Kcounter::size()
{
  return kc->kcontainer_size();
}

void Kcounter::remove(char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  kc->kcontainer_remove(kmer);
}

void Kcounter::add_seq(const char* seq, uint32_t length)
{
  kc->kcontainer_add_seq(seq, length, merge_func);
}

void Kcounter::parallel_add_seq(const char* seq, uint32_t length) {
  kc->parallel_kcontainer_add_seq(seq, length);
}
