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
        .def("add", &Kset::insert, R"pbdoc(
            Add a kmer to Kset

            Takes one argument, the kmer represented as a string
          )pbdoc")
        .def("__contains__", &Kset::contains, R"pbdoc(
            Checks if a kmer is in Kset

            Takes one argument, the kmer represented as a string and returns
            True if kmer is present in Kset or False if it is not present in Kset
          )pbdoc")
        .def("clear", &Kset::clear, R"pbdoc(
            Clears the Kset
          )pbdoc")
        .def("__len__", &Kset::size, R"pbdoc(
            Returns the number of kmers in Kset
          )pbdoc")
        .def("__delitem__", &Kset::remove, R"pbdoc(
            Removes a kmer from Kset

            Takes one argument, the kmer represented as a string and removes it 
            from the Kset
          )pbdoc");
}


