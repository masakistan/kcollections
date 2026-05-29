#define PYTHON

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "Kdict.h"
#include "globals.h"

namespace py = pybind11;

template<typename T>
void declare_kdict_member(py::module& m, const std::string& typestr) {
  using CClass = Kdict<T>;
  std::string pyclass_name = std::string("Kdict_") + typestr;

  py::class_<CClass>(m, pyclass_name.c_str(), py::dynamic_attr())
      .def(py::init<>())
      .def(py::init<const int>())
      .def("write", &CClass::write)
      .def("read", &CClass::read)
      .def("__setitem__", &CClass::add)
      .def("__getitem__", &CClass::get, py::return_value_policy::reference)
      .def("__iter__", [](CClass& v) { return py::make_iterator(v.begin(), v.end()); },
           py::keep_alive<0, 1>())
      .def("__contains__", &CClass::contains)
      .def("clear", &CClass::clear)
      .def("__len__", &CClass::size)
      .def("__delitem__", &CClass::remove)
      .def_property_readonly("k", &CClass::get_k)
      .def("parallel_add_init", &CClass::parallel_add_init, py::call_guard<py::gil_scoped_release>())
      .def("parallel_add_seq", [](CClass& kd, const char* seq, py::iterable& iter) {
        kd.parallel_add_seq(seq, iter);
      })
      .def("add_seq", [](CClass& kd, const char* seq, py::iterable& iter) {
        kd.add_seq(seq, iter);
      })
      .def("parallel_add", &CClass::parallel_add)
      .def("parallel_add_join", &CClass::parallel_add_join, py::call_guard<py::gil_scoped_release>())
      .def("set_merge_func", &CClass::set_merge_func)
      .def("_trie_stats", [](CClass& kd) {
        auto* root = kd.get_root();
        py::dict stats;
        stats["root_vs_size"] = kd.get_vs_size(root);
        stats["root_uc_size"] = kd.get_uc_size(root);
        return stats;
      });
}

template<typename T>
void make_opaque(py::module& m, const std::string& name) {
  py::bind_vector<std::vector<T>>(m, ("ovector_" + name).c_str(), py::module_local());
}

template<typename T>
void declare_kdict(py::module& m, const std::string& name) {
  declare_kdict_member<T>(m, name);
  declare_kdict_member<std::vector<T>>(m, std::string("vector_") + name);
}

void bind_kdict(py::module& m) {
  make_opaque<int>(m, "int");
  make_opaque<float>(m, "float");
  make_opaque<char>(m, "bool");
  make_opaque<std::string>(m, "str");
  make_opaque<std::vector<int>>(m, "vector_int");
  make_opaque<std::vector<float>>(m, "vector_float");
  make_opaque<std::vector<char>>(m, "vector_bool");
  make_opaque<std::vector<std::string>>(m, "vector_str");

  declare_kdict<int>(m, "int");
  declare_kdict<float>(m, "float");
  declare_kdict<char>(m, "bool");
  declare_kdict<std::string>(m, "str");
}
