#include <fstream>
#include <iostream>
#include "Kdict.h"
#include <sstream>
#include <string>


int main( int argc, char* argv[] )
{
    int k = atoi( argv[ 1 ] );
    Kdict* kd = new Kdict( k );
    std::ifstream infile( argv[ 2 ] );

    int c = 0;
    std::string line;
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        //std::cout << line.substr( 0, k ) << std::endl;
        strcpy( kmer, line.substr( 0, k ).c_str() );

        kd->insert( kmer );
        std::cout << "inserted correctly? ";
        if( kd->contains( kmer ) )
            std::cout << "True!" << std::endl;
        else
            std::cout << "False!" << std::endl;

        c++;
    }
    std::cout << "Processed " << c << " kmers!" << std::endl;
}


