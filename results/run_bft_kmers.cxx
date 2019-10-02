extern "C" {
    #include <bft.h>
}
#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>

/*
The program expects 1 argument, the kmer size.
The program will read from stdin kmers and its count on every line. The count is ignored.
Pipe the kmer file into the program and provide the kmer size.
 */
int main(int argc, char* argv[]) {
    if(argc < 2) {
        std::cout << "Error: missing argument of kmer size." << std::endl;
        return 1;
    }

    int kmerSize = atoi(argv[1]);
    BFT* bft = create_cdbg(kmerSize, 1);
    char** kmerArray = new char*[1];
    int count;
    std::string kmer;
    std::cin >> kmer >> count;
    std::replace(kmer.begin(), kmer.end(), 'N', 'A');
    std::replace(kmer.begin(), kmer.end(), 'n', 'a');
    kmerArray[0] = (char*) kmer.c_str();
    insert_kmers_new_genome(1, kmerArray, "genome", bft);
    while(std::cin >> kmer >> count) {
        std::replace(kmer.begin(), kmer.end(), 'N', 'A');
        std::replace(kmer.begin(), kmer.end(), 'n', 'a');
        kmerArray[0] = (char*) kmer.c_str();
        insert_kmers_last_genome(1, kmerArray, bft);
    }

    free_cdbg(bft);
    delete[] kmerArray;

    std::cout << "Successfully finished!" << std::endl;

    return 0;
}
