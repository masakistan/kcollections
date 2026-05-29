#define PYTHON

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Kset.h"

namespace py = pybind11;

static void add_seq_any(Kset& ks, py::object seq) {
  std::string owned;
  const char* data = nullptr;
  size_t len = 0;

  if (py::isinstance<py::str>(seq)) {
    owned = seq.cast<std::string>();
    data = owned.data();
    len = owned.size();
  } else if (py::isinstance<py::bytes>(seq) || py::isinstance<py::bytearray>(seq)) {
    owned = py::cast<std::string>(seq);
    data = owned.data();
    len = owned.size();
  } else if (PyObject_CheckBuffer(seq.ptr())) {
    Py_buffer view;
    if (PyObject_GetBuffer(seq.ptr(), &view, PyBUF_SIMPLE) != 0) {
      throw py::type_error("add_seq: could not obtain buffer view");
    }
    owned.assign(static_cast<const char*>(view.buf), static_cast<size_t>(view.len));
    PyBuffer_Release(&view);
    data = owned.data();
    len = owned.size();
  } else {
    throw py::type_error("add_seq expects str, bytes, bytearray, or buffer");
  }

  py::gil_scoped_release release;
  ks.add_seq(data, len);
}

void bind_kset(py::module& m) {
  py::class_<Kset>(m, "Kset")
      .def(py::init<>())
      .def(py::init<const int>())
      .def("write", &Kset::write)
      .def("read", &Kset::read)
      .def("__iter__", [](Kset& v) { return py::make_iterator(v.begin(), v.end()); },
           py::keep_alive<0, 1>())
      .def("add", &Kset::add)
      .def("__contains__", &Kset::contains)
      .def("clear", &Kset::clear)
      .def("__len__", &Kset::size)
      .def("__delitem__", &Kset::remove)
      .def("add_seq", &add_seq_any)
      .def("parallel_add_init", &Kset::parallel_add_init, py::call_guard<py::gil_scoped_release>())
      .def("parallel_add", &Kset::parallel_add)
      .def("parallel_add_seq", &Kset::parallel_add_seq, py::call_guard<py::gil_scoped_release>())
      .def("parallel_add_join", &Kset::parallel_add_join, py::call_guard<py::gil_scoped_release>())
      .def_property_readonly("k", &Kset::get_k);
}
