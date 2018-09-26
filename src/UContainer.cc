#include "UContainer.h"

UContainer::UContainer() : Container()
{
    m_bkmers = new std::set< Bkmer >();
}

UContainer::~UContainer()
{
    delete m_bkmers;
}

bool UContainer::contains( Bkmer* bkmer )
{
    std::set< Bkmer >::iterator index = m_bkmers->find( *bkmer );
    if( index == m_bkmers->end() )
    {
        return false;
    }
    return true;
}

void UContainer::insert( Bkmer* bkmer )
{
    std::set< Bkmer >::iterator index = m_bkmers->find( *bkmer );
    if( index == m_bkmers->end() )
    {
        m_bkmers->insert( *bkmer );
    }
}

bool UContainer::is_full()
{
    if( m_bkmers->size() == s_capacity )
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
    return m_bkmers;
}


