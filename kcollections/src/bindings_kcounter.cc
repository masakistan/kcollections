#define PYTHON

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Kcounter.h"

namespace py = pybind11;

static void kcounter_add_seq_any(Kcounter& kc, py::object seq) {
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
  kc.add_seq(data, len);
}

void bind_kcounter(py::module& m) {
  py::class_<Kcounter>(m, "Kcounter")
      .def(py::init<>())
      .def(py::init<const int>())
      .def("write", &Kcounter::write)
      .def("read", &Kcounter::read)
      .def("__setitem__", &Kcounter::insert)
      .def("__getitem__", &Kcounter::get)
      .def("__iter__", [](Kcounter& v) { return py::make_iterator(v.begin(), v.end()); },
           py::keep_alive<0, 1>())
      .def("__contains__", &Kcounter::contains)
      .def("clear", &Kcounter::clear)
      .def("__len__", &Kcounter::size)
      .def("__delitem__", &Kcounter::remove)
      .def("add_seq", &kcounter_add_seq_any)
      .def("parallel_add_init", &Kcounter::parallel_add_init)
      .def("parallel_add", &Kcounter::parallel_add)
      .def("parallel_add_seq", &Kcounter::parallel_add_seq, py::call_guard<py::gil_scoped_release>())
      .def("parallel_add_join", &Kcounter::parallel_add_join, py::call_guard<py::gil_scoped_release>())
      .def_property_readonly("k", &Kcounter::get_k);
}
