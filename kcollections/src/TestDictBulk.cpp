#include <iostream>
#include "Kdict.h"
#include <fstream>
#include <string>

using namespace std;

class Indices {
private:
  size_t s;
  int* indices;
public:
  Indices() : s(0), indices(NULL) {}
  Indices(int start) : s(1) {
    indices = (int*) malloc(sizeof(int));
    indices[0] = start;
  }
  ~Indices() {
    if(s > 0) {
      free(indices);
    }
  }
  
  void add(int index) {
    if(s == 0) {
      indices = (int*) malloc(sizeof(int));
    } else {
      indices = (int*) realloc(indices, sizeof(int) * (s + 1));
    }

    indices[s] = index;
    s++;
  }

  size_t size() { return s; }
  int* get_indices() { return indices; }
  int get(int i) { return indices[i]; }
};

Indices& merge_func(Indices& o, Indices& n) {
  for(int i = 0; i < n.size(); i++) {
    o.add(n.get(i));
  }
  //delete &n;
  return o;
}

std::vector<int>& merge_func_vector(std::vector<int>& o, std::vector<int>& n) {
  //std::cout << "merging " << o.size() << "\t" << n.size() << std::endl;
  for(auto i : n) {
    o.push_back(i);
  }
  return o;
}

class Position {
private:
  int pos;
  //Indices temp;
public:
  typedef std::input_iterator_tag iterator_category;    // iterator category
  typedef std::vector<int> value_type;           // value type
  typedef std::ptrdiff_t difference_type;           // difference type
  typedef std::vector<int>* pointer;
  typedef std::vector<int>& reference;
  
  Position(int pos = 0) : pos(pos) {}

  Position& operator++() {
    pos++;
    return *this;
  }

  Position operator++(int) {
    Position& old = *this;
    pos++;
    return old;
  }

  std::vector<int> operator*() {
    std::vector<int> temp = {pos};
    //std::cout << "creating temp\t" << temp.size() << std::endl;
    return temp;
  }

  pointer operator->() {
    std::vector<int> temp = {pos};
    return &temp;
  }

  Position& begin() {
    return *this;
  }
};

int main(int argc, char* argv[]) {
  int k = atoi(argv[1]);
  char* file_path = argv[2];

  string line;
  string seq;
  ifstream fh(file_path);
  //Kdict<Indices>* kc = new Kdict<Indices>(k);
  Kdict<std::vector<int>>* kc = new Kdict<std::vector<int>>(k);
  kc->set_merge_func(merge_func_vector);
  Position pos;

  kc->parallel_add_init(4);

  int c = 0;
  if(fh.is_open()) {
    while(getline(fh, line)) {
      //std::cout << "inserting: " << line << std::endl;
      //kc->add(line.c_str());
      if(line.c_str()[0] == '>') {
        continue;
      }
      seq.append(line);
    }

    fh.close();
    std::cout << "seq size: " << seq.size() << std::endl;
    kc->parallel_add_seq(seq.c_str(), pos);

  } else {
    std::cout << "Could not open file: " << file_path << std::endl;
  }

  kc->parallel_add_join();
  
  std::cout << "Kmer set contains " << kc->size() << " kmers" << std::endl;

  int count = 0;
  int total = kc->size();
  
  std::cout << "Printing all kmers!" << std::endl;
  for(auto& it : *kc) {
    if(it.second->size() > 7) {
      std::cout << count << " / " << total << "\tkmer: " << it.first << ":" << std::endl;
      for(auto i : *it.second) {
	std::cout << "\t" << i << std::endl;
      }
    }
    count += 1;
  }
  
  /*fh.open(file_path);
  if(fh.is_open()) {
    while(getline(fh, line)) {
      kc->remove(line.c_str());
      //std::cout << "kcet size: " << kc->size() << std::endl;
    }
  }*/

  std::cout << "\tDone printing all kmers!" << std::endl;
  std::cout << "Kmer set contains " << kc->size() << " kmers" << std::endl;
  delete kc;

  return 0;
}
