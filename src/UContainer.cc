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
    //printf("uc insert:  %p\n", bkmer->get_bseq() );
    //m_bkmers->insert( *bkmer ); 
    std::vector< Bkmer* >::iterator it =
        std::upper_bound( m_bkmers->begin(), m_bkmers->end(), bkmer, compare_bkmer() );

    /*std::cout << "stuff before\t" << m_bkmers->size() << std::endl;
    for( Bkmer& tbkmer : *m_bkmers )
    {
        std::cout << tbkmer.get_seq() << std::endl << std::flush;
    }
    std::cout << "stuff after" << std::endl;

    printf("uc after ub search:  %p\n", bkmer->get_bseq() );
    std::cout << bkmer->get_seq() << "\tupb: " << ( it == m_bkmers->begin() ) << "\t" << ( it == m_bkmers->end() ) << std::endl << std::flush;

    printf("uc before insert:  %p\n", bkmer->get_bseq() );*/
    
    m_bkmers->insert( it, new Bkmer( *bkmer ) );
    /*printf("uc after insert:  %p\n", bkmer->get_bseq() );
    
    std::cout << "*************************\t" << m_bkmers->size()<< std::endl << std::flush;

    for( int i = 0; i < m_bkmers->size(); i++ )
    {
        std::cout << m_bkmers->at( i ).get_bseq() << std::endl << std::flush;
    }

    for( Bkmer& tbkmer : *m_bkmers )
    {
        std::cout << tbkmer.get_seq() << std::endl << std::flush;
        printf("\taddr:  %p\n", tbkmer.get_bseq() );
    }*/
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


