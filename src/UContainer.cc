#include "UContainer.h"

UContainer::UContainer() : Container()
{
    bkmers = new std::set< Bkmer >();
}

void UContainer::insert( Bkmer bkmer )
{
    //bkmers->insert( bkmer );
}

bool UContainer::contains_kmer( Bkmer* bkmer )
{
    //if( bkmers->find( *bkmer ) != bkmers->end() )
    {
        return true;
    }

    return false;
}

bool UContainer::is_full()
{
    if( bkmers->size() == s_capacity )
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::set< Bkmer >* UContainer::get_bkmers()
{
    return bkmers;
}


