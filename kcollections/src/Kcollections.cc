#if KDICT
#include <vector>
#include <set>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include "Kdict.h"


template<typename T>
void declare_kdict_member(py::module &m, const std::string &typestr) {
  using CClass = Kdict<T>;
  using VClass = Vertex<T>;
  std::string pyclass_name = std::string("Kdict_") + typestr;
  
  m.doc() = R"pbdoc(
	kcollections python bindings
	----------------------------

	.. currentmodule:: Kcollections

	.. autosummary::
	   :toctree: _generate

	   Kdict
      )pbdoc";
  
  py::class_<CClass>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    .def(py::init<const int>())
    .def("write", &CClass::write, R"pbdoc(Write a kset to disk)pbdoc")
    .def("read", &CClass::read, R"pbdoc(Read a kset from disk)pbdoc")
    .def("__setitem__", &CClass::add, R"pbdoc(
	    Add a kmer to Kdict

	    Takes two arguments, the kmer represented as a string and the object to set it to.
	  )pbdoc")
    .def("__getitem__", &CClass::get, R"pbdoc()pbdoc")
    .def("__iter__", [](CClass& v) { return py::make_iterator(v.begin(), v.end()); })
    .def("__contains__", &CClass::contains, R"pbdoc(
	    Checks if a kmer is in Kdict

	    Takes one argument, the kmer represented as a string and returns
	    True if kmer is present in Kdict or False if it is not present in Kdict
	  )pbdoc")
    .def("clear", &CClass::clear, R"pbdoc(
	    Clears the Kdict
	  )pbdoc")
    .def("__len__", &CClass::size, R"pbdoc(
	    Returns the number of kmers in Kdict
	  )pbdoc")
    .def("__delitem__", &CClass::remove, R"pbdoc(
	    Removes a kmer from Kdict

	    Takes one argument, the kmer represented as a string and removes it
	    from the Kdict
	  )pbdoc")
    .def("get_uc_kmer", &CClass::get_uc_kmer)
    .def("get_uc_size", &CClass::get_uc_size)
    .def("get_root", &CClass::get_root, py::return_value_policy::reference)
    .def("get_vs_size", &CClass::get_vs_size )
    .def("get_child_vertex", &CClass::get_child_vertex, py::return_value_policy::reference )
    .def("get_child_suffix", &CClass::get_child_suffix )
    .def_property_readonly("k", &CClass::get_k)
    .def("parallel_add_init", &CClass::parallel_add_init, py::call_guard<py::gil_scoped_release>())
    .def("parallel_add", &CClass::parallel_add)
    .def("parallel_add_seq", &CClass::parallel_add_seq)
    .def("add_seq", &CClass::add_seq)
    .def("parallel_add_join", &CClass::parallel_add_join, py::call_guard<py::gil_scoped_release>())
    .def("set_merge_func", &CClass::set_merge_func);
  
  std::string vertex_pyclass_name = std::string("Vertex_") + typestr;
  py::class_<VClass, std::shared_ptr<VClass>>(m, vertex_pyclass_name.c_str(), py::module_local())
    .def( "uc", &VClass::get_uc );
}

template<typename T>
void declare_kdict(py::module& m, const std::string& name) {
  
  declare_kdict_member<T>(m, name);
  declare_kdict_member<std::vector<T>>(m, std::string("vector_") + name);
  declare_kdict_member<std::set<T>>(m, std::string("set_") + name);
  declare_kdict_member<std::list<T>>(m, std::string("list_") + name);
}

//PYBIND11_MAKE_OPAQUE(std::vector<int>);

PYBIND11_MODULE( _Kdict, m )
{
  declare_kdict<int>(m, "int");
  declare_kdict<float>(m, "float");
  // NOTE: we use char because std::vector<bool> does not return references to items
  // using char instead is a hack, should we change this?
  declare_kdict<char>(m, "bool");
  declare_kdict<std::string>(m, "string");
  declare_kdict<py::object>(m, "object");

  //declare_kdict_member<py::list>(m, "pylist");
  declare_kdict_member<std::vector<std::vector<int>>>(m, "list_list");
}

#elif KCOUNTER
#include "Kcounter.h"

PYBIND11_MODULE( _Kcounter, m ) {
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
      .def("write", &Kcounter::write, R"pbdoc(Write a kset to disk)pbdoc")
      .def("read", &Kcounter::read, R"pbdoc(Read a kset from disk)pbdoc")
      .def("__setitem__", &Kcounter::insert, R"pbdoc(
	    Add a kmer to Kcounter

	    Takes two arguments, the kmer represented as a string and an int associated with the kmer.
	  )pbdoc")
      .def("__getitem__", &Kcounter::get, R"pbdoc()pbdoc")
      .def("items", [](Kcounter& v) { return py::make_iterator(v.begin(), v.end()); })
      .def("iteritems", [](Kcounter& v) { return py::make_iterator(v.begin(), v.end()); })
      .def("__iter__", [](Kcounter& v) { return py::make_iterator(v.begin(), v.end()); })
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
    
    py::class_<Vertex<int>, std::shared_ptr<Vertex<int>>>(m, "Vertex", py::module_local())
      .def( "vs_size", &Vertex<int>::get_vs_size)
      .def( "uc", &Vertex<int>::get_uc );
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
      .def("write", &Kset::write, R"pbdoc(Write a kset to disk)pbdoc")
      .def("read", &Kset::read, R"pbdoc(Read a kset from disk)pbdoc")
      .def("__iter__", [](Kset& v) { return py::make_iterator(v.begin(), v.end()); })
      .def("add", &Kset::add, R"pbdoc(
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
      .def("get_vs_size", &Kset::get_vs_size )
      .def("get_child_vertex", &Kset::get_child_vertex, py::return_value_policy::reference )
      .def("get_child_suffix", &Kset::get_child_suffix )
      .def("add_seq", &Kset::add_seq)
      .def("parallel_add_init", &Kset::parallel_add_init)
      .def("parallel_add", &Kset::parallel_add)
      .def("parallel_add_seq", &Kset::parallel_add_seq)
      .def("parallel_add_join", &Kset::parallel_add_join)
      .def_property_readonly("k", &Kset::get_k);
    
    py::class_<Vertex, std::shared_ptr< Vertex>>(m, "Vertex", py::module_local())
      .def( "get_vs_size", &Vertex::get_vs_size)
      .def( "get_uc", &Vertex::get_uc );
  }

#endif



