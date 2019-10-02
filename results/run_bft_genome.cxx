extern "C" {
    #include <bft.h>
}
#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>

/*
The program expects 1 argument, the kmer size.
The kmers is passed in on stdin with each genome on a new line.
Given a genome in FASTA format, remove the headers and make each entry one line
with the following command: `awk '!/^>/ { printf "%s", $0; } /^>/ { print "" }' genome.fasta`
and pipe that into `./bft_genome KMER_SIZE`
 */
int main(int argc, char* argv[]) {
    if(argc < 2) {
        std::cout << "Error: missing argument of kmer size." << std::endl;
        return 1;
    }

    int kmerSize = atoi(argv[1]);
    BFT* bft = create_cdbg(kmerSize, 1);
    char** kmerArray = new char*[1];
    std::string seq;
    std::cin >> seq;
    std::replace(seq.begin(), seq.end(), 'N', 'A');
    std::replace(seq.begin(), seq.end(), 'n', 'a');
    kmerArray[0] = (char*) seq.c_str();
    insert_kmers_new_genome(1, kmerArray, "genome", bft);
    while(std::cin >> seq) {
        std::replace(seq.begin(), seq.end(), 'N', 'A');
        std::replace(seq.begin(), seq.end(), 'n', 'a');
        for(int i = 0; i < seq.length() - kmerSize; i++) {
            kmerArray[0] = (char*) seq.substr(i, kmerSize).c_str();
            insert_kmers_last_genome(1, kmerArray, bft);
        }
    }

    free_cdbg(bft);
    delete[] kmerArray;

    std::cout << "Successfully finished!" << std::endl;

    return 0;
}
