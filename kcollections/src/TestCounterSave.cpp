#include <iostream>
#include "Kcounter.h"
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  int k = atoi(argv[1]);
  char* file_path = argv[2];

  string line;
  string seq;
  ifstream fh(file_path);
  Kcounter* kc = new Kcounter(k);
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

  int before_count = 0;
  int before_total = kc->size();
  
  for(auto& it : *kc) {
    std::cout << it.first << "\t" << *it.second << std::endl;
    before_count += 1;
  }

  if(before_count != kc->size()) {
    std::cout << "ERROR: iterating and traversed kmer counts are not equal" << std::endl;
  }

  
  std::cout << "\n\nSaving kset using boost..." << std::endl;
  kc->write("testsave.bs");
  std::cout << "Done using boost. Deleting kset..." << std::endl;
  delete kc;
  std::cout << "\tDone!" << std::endl;


  std::cout << "\n\nLoading kset using boost..." << std::endl;
  kc = new Kcounter();
  kc->read("testsave.bs");
  std::cout << "\tDone!" << std::endl;
  std::cout << "\n\nVerifying loaded kset..." << std::endl;
  std::cout << "\tk: " << kc->get_k() << std::endl;
  if(before_total == kc->size()){
    std::cout << "\tSUCCESS: current size matches previous size";
  } else {
    std::cout << "\tERROR: current size does not match previous size";
  }
  std::cout << " (current: " << kc->size() << ", previous: " << before_total << ")" << std::endl;
  std::cout << "\tDone!" << std::endl;

  int after_count = 0;
  for(auto& it : *kc) {
    std::cout << it.first << "\t" << *it.second << std::endl;
    after_count += 1;
  }

  std::cout << "\n\nTesting kmer iteration" << std::endl;
  if(after_count == before_count) {
    std::cout << "\tSUCCESS: kmer iteration mounts match";
  } else {
    std::cout << "\tERROR: kmer iteration counts failed";
  }

  std::cout << " (before: " << before_count << ", after: " << after_count << ")" << std::endl;

  /*fh.open(file_path);
  if(fh.is_open()) {
    while(getline(fh, line)) {
      kc->remove(line.c_str());
      //std::cout << "kcet size: " << kc->size() << std::endl;
    }
  }*/


  return 0;
}
