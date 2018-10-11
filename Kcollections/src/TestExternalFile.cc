#include <fstream>
#include <iostream>
#include "Kdict.h"
#include <sstream>
#include <string>


int main( int argc, char* argv[] )
{
    int c = 0, interval = 100000;
    int k = atoi( argv[ 1 ] );
    Kdict* kd = new Kdict( k );
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );
   
    std::cout << "Inserting kmers";
    std::ifstream infile( argv[ 2 ] );
    std::string line;
    while( std::getline( infile, line ) )
    {
        std::istringstream iss(line);
        strcpy( kmer, line.substr( 0, k ).c_str() );

        kd->insert( kmer );
        /*std::cout << "inserted correctly? ";
        if( kd->contains( kmer ) )
            std::cout << "True!" << std::endl;
        else
            std::cout << "False!" << std::endl;*/
        c++;
        if( c % interval == 0 )
        {
            std::cout << '.' << std::flush;
        }
    }
    infile.close();
    std::cout << "\nFinished inserting all kmers! Processed " << c << " kmers!" << std::endl << std::endl;;

    std::cout << "Verifing kmers were inserted into container";
    infile = std::ifstream( argv[ 2 ] );
    c = 0;
    while( std::getline( infile, line ) )
    {
        std::istringstream iss(line);
        strcpy( kmer, line.substr( 0, k ).c_str() );
        assert( kd->contains( kmer ) );

        c++;
        if( c % interval == 0 )
        {
            std::cout << '.' << std::flush;
        }
    }
    infile.close();
    std::cout << "\nFinished checking if container contains  all kmers!" << std::endl;
    delete kmer;
    delete kd;
}


