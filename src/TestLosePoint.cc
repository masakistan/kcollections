#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>

#include "set/Kcontainer.h"


int main( int argc, char* argv[] )
{
    int k = atoi( argv[ 1 ] );
    Kcontainer* kd = create_kcontainer( k );
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );

    uint64_t c = 0;
    int interval = 50000;
    std::cout << "Inserting kmers" << std::endl;
    std::ifstream infile( argv[ 2 ] );
    std::string line;
    while( std::getline( infile, line ) )
    {
        std::istringstream iss(line);
        //std::cout << "inserting kmer: " << line << "\t" << c << std::endl << std::flush;
        strcpy( kmer, line.substr( 0, k ).c_str() );

        kcontainer_insert( kd, kmer );
        //assert( contains( kd, kmer ) );

        c++;
        if( c % interval == 0 )
        {
            std::cout << '\r';
            std::cout << "\tinserted " << c << " kmers"<< std::flush;
        }
    }
    std::cout << "\nFinished inserting all kmers! Processed " << c << " kmers!" << std::endl << std::endl;;
    std::cout << "Total kmers in kcontainer: " << std::endl;
    infile.close();

    c = 0;
    std::cout << "Checking contains\n";
    infile = std::ifstream( argv[ 2 ] );
    while( std::getline( infile, line ) )
    {
        std::istringstream iss(line);
        strcpy( kmer, line.substr( 0, k ).c_str() );
        //std::cout << "searching for: " << line << std::endl;

        assert( kcontainer_contains( kd, kmer ) );
        /*std::cout << "inserted correctly? ";
        if( kd->contains( kmer ) )
            std::cout << "True!" << std::endl;
        else
            std::cout << "False!" << std::endl;*/
        c++;
        if( c % interval == 0 )
        {
            std::cout << '\r';
            std::cout << "\tchecked " << c << " kmers"<< std::flush;
        }
    }
    std::cout << "\tFinished checking contains all kmers! Processed " << c << " kmer!" << std::endl << std::endl << std::flush;
    infile.close();

    //print( &( kd->v.uc ), k, 0 );

    free( kmer );
    free_kcontainer( kd );

    return 0;
}


