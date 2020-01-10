#if KDICT
#include "Kdict.h"

PYBIND11_MODULE( _Kdict, m )
{
    m.doc() = R"pbdoc(
        kcollections python bindings
        ----------------------------

        .. currentmodule:: Kcollections

        .. autosummary::
           :toctree: _generate

           Kdict
      )pbdoc";

    py::class_<Kdict>(m, "Kdict")
        .def(py::init<const int>())
        .def("__setitem__", &Kdict::add, R"pbdoc(
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
          )pbdoc")
        .def("get_uc_kmer", &Kdict::get_uc_kmer)
        .def("get_uc_size", &Kdict::get_uc_size)
        .def("get_root", &Kdict::get_root, py::return_value_policy::reference)
        .def("get_vs_size", &Kdict::get_vs_size )
        .def("get_child_vertex", &Kdict::get_child_vertex, py::return_value_policy::reference )
        .def("get_child_suffix", &Kdict::get_child_suffix )
        .def_property_readonly("k", &Kdict::get_k);

    py::class_<Vertex>(m, "Vertex")
        .def_readonly( "uc", &Vertex::uc );

}
#elif KCOUNTER
#include "Kcounter.h"

PYBIND11_MODULE( _Kcounter, m )
{
    m.doc() = R"pbdoc(
        kcollections python bindings
        ----------------------------

        .. currentmodule:: Kcollections

        .. autosummary::
           :toctree: _generate

           Kcounter
      )pbdoc";

    py::class_<Kcounter>(m, "Kcounter")
        .def(py::init<const int>())
        .def("__setitem__", &Kcounter::insert, R"pbdoc(
            Add a kmer to Kcounter

            Takes two arguments, the kmer represented as a string and an int associated with the kmer.
          )pbdoc")
        .def("__getitem__", &Kcounter::get, R"pbdoc()pbdoc")
        .def("__contains__", &Kcounter::contains, R"pbdoc(
            Checks if a kmer is in Kcounter

            Takes one argument, the kmer represented as a string and returns
            True if kmer is present in Kcounter or False if it is not present in Kcounter
          )pbdoc")
        .def("clear", &Kcounter::clear, R"pbdoc(
            Clears the Kcounter
          )pbdoc")
        .def("__len__", &Kcounter::size, R"pbdoc(
            Returns the number of kmers in Kcounter
          )pbdoc")
        .def("__delitem__", &Kcounter::remove, R"pbdoc(
            Removes a kmer from Kcounter

            Takes one argument, the kmer represented as a string and removes it
            from the Kcounter
          )pbdoc")
        .def("get_uc_kmer", &Kcounter::get_uc_kmer)
        .def("get_uc_size", &Kcounter::get_uc_size)
        .def("get_root", &Kcounter::get_root, py::return_value_policy::reference)
        .def("get_vs_size", &Kcounter::get_vs_size )
        .def("get_child_vertex", &Kcounter::get_child_vertex, py::return_value_policy::reference )
        .def("get_child_suffix", &Kcounter::get_child_suffix )
        .def("add_seq", &Kcounter::add_seq)
        .def("parallel_add_init", &Kcounter::parallel_add_init)
        .def("parallel_add", &Kcounter::parallel_add)
        .def("parallel_add_seq", &Kcounter::parallel_add_seq)
        .def("parallel_add_join", &Kcounter::parallel_add_join)
        .def_property_readonly("k", &Kcounter::get_k);

    py::class_<Vertex>(m, "Vertex")
        .def_readonly( "vs_size", &Vertex::vs_size)
        .def_readonly( "uc", &Vertex::uc );
}
#elif KSET
#include "Kset.h"
#include <pybind11/stl.h>

PYBIND11_MODULE( _Kset, m )
{
    m.doc() = R"pbdoc(
        kcollections python bindings
        ----------------------------

        .. currentmodule:: Kcollections

        .. autosummary::
           :toctree: _generate

           Kset
      )pbdoc";

    m.def("testbit32", &testbit32);
    
    py::class_<PgData>(m, "PgData")
      .def_readwrite("coords", &PgData::coords, py::return_value_policy::reference)
      .def_readwrite("genomes", &PgData::genomes, py::return_value_policy::reference)
      //.def_readwrite("counts", &PgData::counts, py::return_value_policy::reference)
      .def_readwrite("vidx", &PgData::vidx, py::return_value_policy::reference)
      .def_readwrite("orientation", &PgData::orientation, py::return_value_policy::reference)
      .def_readwrite("size", &PgData::size, py::return_value_policy::reference);
    
    py::class_<Kset>(m, "Kset")
      .def(py::init<const int>())
      /*.def("add", &Kset::add, R"pbdoc(
	Add a kmer to Kset
	
	Takes one argument, the kmer represented as a string
	)pbdoc")*/
      .def("__getitem__", &Kset::get, R"pbdoc()pbdoc", py::return_value_policy::reference)
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
          )pbdoc")
      .def("get_uc_kmer", &Kset::get_uc_kmer)
      .def("get_uc_size", &Kset::get_uc_size)
      .def("get_root", &Kset::get_root, py::return_value_policy::reference)
      .def("get_vs_size", &Kset::get_vs_size )
      .def("get_child_vertex", &Kset::get_child_vertex, py::return_value_policy::reference )
      .def("get_child_suffix", &Kset::get_child_suffix )
      //.def("add_seq", &Kset::add_seq)
      .def("parallel_add_init", &Kset::parallel_add_init)
      .def("parallel_add", &Kset::parallel_add)
      .def("parallel_add_seq", &Kset::parallel_add_seq)
      .def("parallel_add_join", &Kset::parallel_add_join)
      .def_property_readonly("k", &Kset::get_k);
    
    py::class_<Vertex>(m, "Vertex")
      .def_readonly( "vs_size", &Vertex::vs_size)
      .def_readonly( "uc", &Vertex::uc );
}
#endif



