#include "Kset.h"

Kset::Kset( int k )
{
  kc = new Kcontainer(k);
  m_k = k;
}

Kset::~Kset()
{
  delete kc;
}

void Kset::clear()
{
  delete kc;
  kc = new Kcontainer(m_k);
}

void Kset::add(const char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kset");
  kc->kcontainer_add(kmer);
}

bool Kset::contains(const char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kset");
  return kc->kcontainer_contains(kmer);
}

uint64_t Kset::size()
{
  return kc->kcontainer_size();
}

void Kset::remove(const char* kmer)
{
  CHECK_KMER_LENGTH(kmer, m_k, "Kset");
  kc->kcontainer_remove(kmer);
}

void Kset::add_seq(const char* seq)
{
  kc->kcontainer_add_seq(seq, strlen(seq));
}

void Kset::parallel_add_seq(const char* seq) {
  kc->parallel_kcontainer_add_seq(seq, strlen(seq));
}

