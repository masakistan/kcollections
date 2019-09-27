#include <iostream>
#include "Kset.h"
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  int k = atoi(argv[1]);
  char* file_path = argv[2];

  string line;
  string seq;
  ifstream fh(file_path);
  Kset* kc = new Kset(k);
  //kc->parallel_add_init(4);

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
    kc->add_seq(seq.c_str());

  } else {
    std::cout << "Could not open file: " << file_path << std::endl;
  }

  //kc->parallel_add_join();
  
  std::cout << "Kmer set contains " << kc->size() << " kmers" << std::endl;

  int count = 0;
  int total = kc->size();
  
  std::cout << "Printing all kmers!" << std::endl;
  for(auto& it : *kc) {
    count += 1;
  }
  std::cout << "count: " << count << std::endl;
  
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
