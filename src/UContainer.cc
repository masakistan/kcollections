#include "UContainer.h"

UContainer::UContainer() : Container()
{
    m_bkmers = new std::vector< Bkmer >();
}

UContainer::~UContainer()
{
    delete m_bkmers;
}

bool UContainer::contains( Bkmer* bkmer )
{
    //std::vector< Bkmer >::iterator index = m_bkmers->find( *bkmer );
    std::cout << "searching uc: " << bkmer->get_seq() << "\t" << std::binary_search( m_bkmers->begin(), m_bkmers->end(), *bkmer, compare_bkmer() ) << std::endl << std::flush;
    return std::binary_search( m_bkmers->begin(), m_bkmers->end(), *bkmer, compare_bkmer() );
}

void UContainer::insert( Bkmer* bkmer )
{
    printf("uc insert:  %p\n", bkmer->get_bseq() );
    //m_bkmers->insert( *bkmer ); 
    std::vector< Bkmer >::iterator it =
        std::upper_bound( m_bkmers->begin(), m_bkmers->end(), *bkmer, compare_bkmer() );

    std::cout << "stuff before\t" << m_bkmers->size() << std::endl;
    for( Bkmer& tbkmer : *m_bkmers )
    {
        std::cout << tbkmer.get_seq() << std::endl << std::flush;
    }
    std::cout << "stuff after" << std::endl;

    printf("uc after ub search:  %p\n", bkmer->get_bseq() );
    std::cout << bkmer->get_seq() << "\tupb: " << ( it == m_bkmers->begin() ) << "\t" << ( it == m_bkmers->end() ) << std::endl << std::flush;

    printf("uc before insert:  %p\n", bkmer->get_bseq() );
    
    m_bkmers->insert( m_bkmers->begin(), *bkmer );
    printf("uc after insert:  %p\n", bkmer->get_bseq() );
    
    std::cout << "*************************" << std::endl << std::flush;

    for( Bkmer& tbkmer : *m_bkmers )
    {
        std::cout << tbkmer.get_seq() << std::endl << std::flush;
        printf("\taddr:  %p\n", tbkmer.get_bseq() );
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

std::vector< Bkmer >* UContainer::get_bkmers()
{
    return m_bkmers;
}

void UContainer::remove( Bkmer* bkmer )
{
    //m_bkmers->erase( *bkmer );
    std::vector< Bkmer >::iterator it =
        std::upper_bound( m_bkmers->begin(), m_bkmers->end(), *bkmer, compare_bkmer() );
    if( it != m_bkmers->end() )
    {
        m_bkmers->erase( it );
    }
}


