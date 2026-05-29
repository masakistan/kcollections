#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_kset(py::module& m);
void bind_kdict(py::module& m);
void bind_kcounter(py::module& m);

PYBIND11_MODULE(_kcollections, m) {
  m.doc() = "kcollections: memory-efficient k-mer sets, dicts, and counters";
  bind_kset(m);
  bind_kdict(m);
  bind_kcounter(m);
}
