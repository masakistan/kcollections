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
          )pbdoc")
        .def("get_uc_kmer", &Kdict::get_uc_kmer)
        .def("get_uc_size", &Kdict::get_uc_size)
        .def("get_root", &Kdict::get_root, py::return_value_policy::reference)
        .def("get_cc_size", &Kdict::get_cc_size )
        .def("get_cc_child_size", &Kdict::get_cc_child_size )
        .def("get_cc_child_vertex", &Kdict::get_cc_child_vertex, py::return_value_policy::reference )
        .def("get_cc_child_suffix", &Kdict::get_cc_child_suffix )
        .def_property_readonly("k", &Kdict::get_k);

    /*py::class_<Kcontainer>(m, "DKcontainer")
        .def_readonly( "root", &Kcontainer::v );*/
    
    py::class_<Vertex>(m, "Vertex")
        .def_readonly( "cc_size", &Vertex::cc_size)
        .def_readonly( "uc", &Vertex::uc );
    
    /*py::class_<UC>(m, "DUC")
        .def_readonly("size", &UC::size);

    py::class_<CC>(m, "DCC")
        .def_readonly("size", &CC::size);
    
    m.def("get_kmer_from_uc", &get_uc_kmer);
    m.def("get_cc", &get_cc, py::return_value_policy::reference);*/
}
#elif KSET
#include "Kset.h"

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
          )pbdoc")
        .def("get_uc_kmer", &Kset::get_uc_kmer)
        .def("get_uc_size", &Kset::get_uc_size)
        .def("get_root", &Kset::get_root, py::return_value_policy::reference)
        .def("get_cc_size", &Kset::get_cc_size )
        .def("get_cc_child_size", &Kset::get_cc_child_size )
        .def("get_cc_child_vertex", &Kset::get_cc_child_vertex, py::return_value_policy::reference )
        .def("get_cc_child_suffix", &Kset::get_cc_child_suffix )
        .def_property_readonly("k", &Kset::get_k);

    py::class_<Vertex>(m, "Vertex")
        .def_readonly( "cc_size", &Vertex::cc_size)
        .def_readonly( "uc", &Vertex::uc );
    
    /*py::class_<UC>(m, "SUC")
        .def_readonly("size", &UC::size);

    py::class_<CC>(m, "SCC")
        .def_readonly("size", &CC::size);
    
    m.def("get_kmer_from_uc", &get_uc_kmer);
    m.def("get_cc", &get_cc, py::return_value_policy::reference);*/
}
#endif



