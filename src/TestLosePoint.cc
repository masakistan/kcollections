#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>

#include "Kdict.h"


int main( int argc, char* argv[] )
{
    int k = atoi( argv[ 1 ] );
    Kdict* kd = create_kdict( k );
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );

    int c = 0, interval = 50000;
    std::cout << "Inserting kmers" << std::endl;
    std::ifstream infile( argv[ 2 ] );
    std::string line;
    while( std::getline( infile, line ) )
    {
        std::istringstream iss(line);
        //std::cout << "inserting kmer: " << line << "\t" << c << std::endl;
        strcpy( kmer, line.substr( 0, k ).c_str() );

        insert( kd, kmer );
        //assert( contains( kd, kmer ) );

        c++;
        if( c % interval == 0 )
        {
            std::cout << '\r';
            std::cout << "\tinserted " << c << " kmers"<< std::flush;
        }
    }
    std::cout << "\nFinished inserting all kmers! Processed " << c << " kmers!" << std::endl << std::endl;;
    std::cout << "Total kmers in kdict: " << std::endl;
    infile.close();

    c = 0;
    std::cout << "Checking contains\n";
    infile = std::ifstream( argv[ 2 ] );
    while( std::getline( infile, line ) )
    {
        std::istringstream iss(line);
        strcpy( kmer, line.substr( 0, k ).c_str() );
        //std::cout << "searching for: " << line << std::endl;

        assert( contains( kd, kmer ) );
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
    free_kdict( kd );

    return 0;
}


