extern "C" {
    #include <bft.h>
}
#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>

int main(int argc, char* argv[]) {
    if(argc < 2) {
        std::cout << "Error: missing arguments of kmer size." << std::endl;
        return 1;
    }

    int kmerSize = atoi(argv[1]);
    BFT* bft = create_cdbg(kmerSize, 1);
    for(std::string seq; std::getline(std::cin, seq);) {
        size_t numKmers = seq.length() - kmerSize;
        char** kmers = new char*[numKmers];
        for(int i = 0; i < numKmers; i++) {
            std::string kmer = seq.substr(i, kmerSize);
            std::replace(kmer.begin(), kmer.end(), 'N', 'A');
            // std::cout << kmer << std::endl;
            kmers[i] = new char[kmerSize];
            strcpy(kmers[i], kmer.c_str());
        }
        insert_kmers_new_genome(numKmers, kmers, "genome", bft);

        for(int i = 0; i < numKmers; i++) {
            delete[] kmers[i];
        }
        delete[] kmers;
    }

    free_cdbg(bft);

    return 0;
}
