#include "UContainer.h"

UContainer::UContainer() : Container()
{
    m_bkmers = new std::vector< Bkmer* >();
}

UContainer::~UContainer()
{
    for( Bkmer* bkmer : *m_bkmers )
    {
        delete bkmer;
    }
    m_bkmers->clear();
    delete m_bkmers;
}

bool UContainer::contains( Bkmer* bkmer )
{
    //std::vector< Bkmer >::iterator index = m_bkmers->find( *bkmer );
    //std::cout << "searching uc: " << bkmer->get_seq() << "\t" << std::binary_search( m_bkmers->begin(), m_bkmers->end(), bkmer, compare_bkmer() ) << std::endl << std::flush;
    return std::binary_search( m_bkmers->begin(), m_bkmers->end(), bkmer, compare_bkmer() );
}

void UContainer::insert( Bkmer* bkmer )
{
    std::vector< Bkmer* >::iterator it =
        std::upper_bound( m_bkmers->begin(), m_bkmers->end(), bkmer, compare_bkmer() );
    
    m_bkmers->insert( it, new Bkmer( *bkmer ) );
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

std::vector< Bkmer* >* UContainer::get_bkmers()
{
    return m_bkmers;
}

void UContainer::remove( Bkmer* bkmer )
{
    //m_bkmers->erase( *bkmer );
    //std::cout << "checking for removal" << std::endl << std::flush;
    std::vector< Bkmer* >::iterator it =
        std::lower_bound( m_bkmers->begin(), m_bkmers->end(), bkmer, compare_bkmer() );
    if( it != m_bkmers->end() && **it == *bkmer )
    {
        //std::cout << "removing!" << std::endl << std::flush;
        delete *it;
        m_bkmers->erase( it );
    }
}


