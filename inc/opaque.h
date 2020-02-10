#if defined(KDICT) && defined(PYTHON)
#include <vector>
#include <set>

#include <pybind11/stl_bind.h>
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);
PYBIND11_MAKE_OPAQUE(std::vector<bool>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
PYBIND11_MAKE_OPAQUE(std::vector<pybind11::object>);

PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<bool>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<pybind11::object>>);

/*
PYBIND11_MAKE_OPAQUE(std::list<int>);
PYBIND11_MAKE_OPAQUE(std::list<float>);
PYBIND11_MAKE_OPAQUE(std::list<bool>);
PYBIND11_MAKE_OPAQUE(std::list<std::string>);
PYBIND11_MAKE_OPAQUE(std::list<py::object>);

PYBIND11_MAKE_OPAQUE(std::list<int>>);
PYBIND11_MAKE_OPAQUE(std::list<std::list<float>>);
PYBIND11_MAKE_OPAQUE(std::list<std::list<bool>>);
PYBIND11_MAKE_OPAQUE(std::list<std::list<std::string>>);
PYBIND11_MAKE_OPAQUE(std::list<std::list<py::object>>);
*/
/*
PYBIND11_MAKE_OPAQUE(std::set<int>);
PYBIND11_MAKE_OPAQUE(std::set<float>);
PYBIND11_MAKE_OPAQUE(std::set<bool>);
PYBIND11_MAKE_OPAQUE(std::set<std::string>);
//PYBIND11_MAKE_OPAQUE(std::set<pybind11::object>);
*/

#endif

