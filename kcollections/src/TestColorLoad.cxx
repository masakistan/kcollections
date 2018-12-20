#include <iostream>
#include "Kcolor.h"
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  int k = atoi(argv[1]);
  char* file_path = argv[2];

  string line, seq;
  ifstream fh(file_path);
  Kcolor* kc = new Kcolor(k);
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
    kc->add_seq(seq.c_str(), seq.size(), 1);
    kc->parallel_add_join();

  } else {
    std::cout << "Could not open file: " << file_path << std::endl;
  }

  std::cout << "Kmer set contains " << kc->size() << " kmers" << std::endl;
  delete kc;

  return 0;
}
