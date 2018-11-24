#include <iostream>
#include "Kset.h"
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
  int k = atoi(argv[1]);
  char* file_path = argv[2];

  string line;
  ifstream fh(file_path);
  Kset* ks = new Kset(k);

  int c = 0;
  if(fh.is_open()) {
    while(getline(fh, line)) {
      //std::cout << "inserting: " << line << std::endl;
      ks->insert(line.c_str());
      c++;
    }
    fh.close();
  } else {
    std::cout << "Could not open file: " << file_path << std::endl;
  }

  std::cout << "Kmer set contains " << ks->size() << " kmers" << std::endl;
  delete ks;

  return 0;
}
