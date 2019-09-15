#include "Kdict.h"


Kdict::Kdict( int k )
{
  kc = new Kcontainer<int>(k);
  m_k = k;
}

Kdict::~Kdict()
{
  delete kc;
}

void Kdict::clear()
{
  delete kc;
  kc = new Kcontainer<int>(m_k);
}

void Kdict::add(char* kmer, int obj)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
  kc->kcontainer_add(kmer, obj, merge_func);
}

int Kdict::get(char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
  return kc->kcontainer_get(kmer);
}

bool Kdict::contains(char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
  return kc->kcontainer_contains(kmer);
}

uint64_t Kdict::size()
{
  return kc->kcontainer_size();
}

void Kdict::remove(char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
  kc->kcontainer_remove(kmer);
}

void Kdict::add_seq(const char* seq, uint32_t length, py::iterable values, std::function<int(int, int)> &f)
{
  kc->kcontainer_add_seq(seq, length, values, f);
}

