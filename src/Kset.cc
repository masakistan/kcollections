#include "Kset.h"

Kset::Kset( int k )
{
    kc = create_kcontainer( k );
}

Kset::~Kset()
{
    free_kcontainer( kc );
}

void Kset::clear()
{
    free_kcontainer( kc );
    kc = create_kcontainer( k );
}

void Kset::insert( char* kmer )
{
    kcontainer_insert( kc, kmer );
}

bool Kset::contains( char* kmer )
{
    kcontainer_contains( kc, kmer );
}

uint64_t Kset::size()
{
    return kcontainer_size( kc );
}

void Kset::remove( char* kmer )
{
    kcontainer_remove( kc, kmer );
}

PYBIND11_MODULE( Kcollections, m )
{
    py::class_<Kset>(m, "Kset")
        .def(py::init<const int>())
        .def("add", &Kset::insert)
        .def("__contains__", &Kset::contains)
        .def("clear", &Kset::clear)
        .def("__len__", &Kset::size)
        .def("__delitem__", &Kset::remove);
}


