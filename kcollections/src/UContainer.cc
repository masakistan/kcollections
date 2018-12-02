#include "UContainer.h"

void init_uc( UC* uc )
{
    //std::cout << "initializing uc" << std::endl;
    uc->suffixes = NULL;
#if KDICT
    uc->objs = NULL;
#endif
    uc->size = 0;
}

void print( UC* uc, int k, int depth )
{
    int len = calc_bk( k );
    int idx;
    for( int i = 0; i < uc->size; i++ )
    {
        idx = i * len;
        char* dseq = deserialize_kmer( k, len, &uc->suffixes[ idx ] );
        std::cout << "kmer: " << dseq << std::endl;
        free( dseq );
    }
}

#if KDICT
void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx, py::handle* obj )
#elif KSET
void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx )
#endif
{
    int len = calc_bk( k );
    if( uc->suffixes == NULL )
    {
        uc->suffixes = ( uint8_t* ) calloc( len , sizeof( uint8_t ) );
#if KDICT
        uc->objs = ( py::handle* ) calloc( len , sizeof( py::handle ) );
#endif
    }
    else
    {
        uc->suffixes = ( uint8_t* ) realloc(
                uc->suffixes,
                len * ( uc->size + 1 ) * sizeof( uint8_t )
                );
#if KDICT
        uc->objs = ( py::handle* ) realloc(
                uc->objs,
                ( uc->size + 1 ) * sizeof( py::handle )
                );
#endif
    }

    if( uc->size < CAPACITY )
    {
        int bytes_to_move = ( uc->size - idx ) * len;
        int suffix_idx = idx * len;
        if( bytes_to_move > 0 )
        {
            std::memmove(
                    &uc->suffixes[ suffix_idx + len ],
                    &uc->suffixes[ suffix_idx ],
                    bytes_to_move
                    );
            
#if KDICT
            bytes_to_move = ( uc->size - idx ) * sizeof( py::handle );
            std::memmove(
                    &uc->objs[ idx + 1 ],
                    &uc->objs[ idx ],
                    bytes_to_move
                    );
#endif
        }

#if KDICT
        std::memcpy( &uc->objs[ idx ], obj, sizeof( py::handle ) );
        uc->objs[ idx ].inc_ref();
#endif
        std::memcpy( &uc->suffixes[ suffix_idx ], bseq, len );
        uc->size++;
    }
}

void free_uc( UC* uc )
{
    //std::cout << "Free uc\t" << uc->size << std::endl;
    if( uc->suffixes != NULL )
    {
        //std::cout << "\tremoving suffixes\t" << uc->size << std::endl;
        free( uc->suffixes );
        //std::cout << "\tdone removing suffixes" << std::endl;
#if KDICT
        for( int i = 0; i < uc->size; i++ )
        {
            uc->objs[ i ].dec_ref();
        }
        free( uc->objs );
#endif
    }
}

void uc_remove( UC* uc, int bk, int idx )
{
    int suffix_idx = idx * bk;
    int bytes_to_move = ( uc->size - ( idx + 1 ) ) * bk;
    std::memmove(
            &uc->suffixes[ suffix_idx ],
            &uc->suffixes[ suffix_idx + bk ],
            bytes_to_move
            );
#if KDICT
    bytes_to_move = ( uc->size - ( idx + 1 ) ) * sizeof( py::handle );
    uc->objs[ idx ].dec_ref();
    std::memmove(
            &uc->objs[ idx ],
            &uc->objs[ idx + 1 ],
            bytes_to_move
            );
#endif

    uc->size--;
}

std::pair< bool, int > uc_find( UC* uc, int k, int depth, uint8_t* bseq )
{
    if( uc->suffixes == NULL )
    {
        return std::make_pair( false, uc->size );
    }

    return binary_search( uc->suffixes, uc->size, calc_bk( k ), bseq );
}


