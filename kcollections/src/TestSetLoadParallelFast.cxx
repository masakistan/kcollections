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
  Kset* ks = new Kset(k);

  int c = 0;
  if(fh.is_open()) {
    while(getline(fh, line)) {
      //std::cout << "inserting: " << line << std::endl;
      //ks->add(line.c_str());
      if(line.c_str()[0] == '>') {
        continue;
      }
      seq.append(line);
    }

    fh.close();
    std::cout << "seq size: " << seq.size() << std::endl;
    ks->add_seq(seq.c_str(), seq.size());

  } else {
    std::cout << "Could not open file: " << file_path << std::endl;
  }
  std::cout << "Kmer set contains " << ks->size() << " kmers" << std::endl;

  // remove the kmers
  /*fh.open(file_path);
  if(fh.is_open()) {
    while(getline(fh, line)) {
      ks->remove(line.c_str());
      //std::cout << "kset size: " << ks->size() << std::endl;
    }
  }*/

  std::cout << "Kmer set contains " << ks->size() << " kmers" << std::endl;
  delete ks;

  return 0;
}
