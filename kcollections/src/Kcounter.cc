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

void Kcounter::insert(const char* kmer, count_dtype count)
{
  //std::cout << "insert " << kmer << std::endl;
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  kc->kcontainer_add(kmer, count, overwrite_merge_func );
}

count_dtype Kcounter::get(const char* kmer)
{
  //std::cout << "get " << kmer << std::endl;
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  return kc->kcontainer_get(kmer);
}

bool Kcounter::contains(const char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  return kc->kcontainer_contains(kmer);
}

uint64_t Kcounter::size()
{
  return kc->kcontainer_size();
}

void Kcounter::remove(const char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kcounter");
  kc->kcontainer_remove(kmer);
}

void Kcounter::add_seq(const char* seq)
{
  kc->kcontainer_add_seq(seq, strlen(seq), merge_func);
}

void Kcounter::parallel_add_seq(const char* seq) {
  kc->parallel_kcontainer_add_seq(seq, strlen(seq));
}
