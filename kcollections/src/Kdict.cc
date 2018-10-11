#include "Kdict.h"

Kdict::Kdict( int k )
{
    this->k = k;
    kc = create_kcontainer( k );
}

Kdict::~Kdict()
{
    free_kcontainer( kc );
}

void Kdict::clear()
{
    free_kcontainer( kc );
    kc = create_kcontainer( k );
}

void Kdict::insert( char* kmer, py::object* obj )
{
    kcontainer_insert( kc, kmer, obj );
}

py::object* Kdict::get( char* kmer )
{
    return kcontainer_get( kc, kmer );
}

bool Kdict::contains( char* kmer )
{
    return kcontainer_contains( kc, kmer );
}

uint64_t Kdict::size()
{
    return kcontainer_size( kc );
}

void Kdict::remove( char* kmer )
{
    kcontainer_remove( kc, kmer );
}


PYBIND11_MODULE( _Kdict, m )
{
    m.doc() = R"pbdoc(
        kcollections python bindings
        ----------------------------

        .. currentmodule:: Kcollections

        .. autosummary::
           :toctree: _generate

           Kset
      )pbdoc";

    py::class_<Kdict>(m, "Kdict")
        .def(py::init<const int>())
        .def("__setitem__", &Kdict::insert, R"pbdoc(
            Add a kmer to Kdict

            Takes two arguments, the kmer represented as a string and the object to set it to.
          )pbdoc")
        .def("__getitem__", &Kdict::get, R"pbdoc()pbdoc")
        .def("__contains__", &Kdict::contains, R"pbdoc(
            Checks if a kmer is in Kdict

            Takes one argument, the kmer represented as a string and returns
            True if kmer is present in Kdict or False if it is not present in Kdict
          )pbdoc")
        .def("clear", &Kdict::clear, R"pbdoc(
            Clears the Kdict
          )pbdoc")
        .def("__len__", &Kdict::size, R"pbdoc(
            Returns the number of kmers in Kdict
          )pbdoc")
        .def("__delitem__", &Kdict::remove, R"pbdoc(
            Removes a kmer from Kdict

            Takes one argument, the kmer represented as a string and removes it 
            from the Kdict
          )pbdoc");
}


