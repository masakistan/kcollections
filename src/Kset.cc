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
    return kcontainer_contains( kc, kmer );
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
    m.doc() = R"pbdoc(
        kcollections python bindings
        ----------------------------

        .. currentmodule:: Kcollections

        .. autosummary::
           :toctree: _generate

           Kset
    )pbdoc";

    py::class_<Kset>(m, "Kset")
        .def(py::init<const int>())
        .def("add", &Kset::insert)
        .def("__contains__", &Kset::contains)
        .def("clear", &Kset::clear)
        .def("__len__", &Kset::size)
        .def("__delitem__", &Kset::remove);
}


