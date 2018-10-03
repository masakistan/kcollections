#include "Kdict.h"

Kdict::Kdict( int k ) : m_k( k )
{
    m_bk = calc_bk( k );
    root = new Vertex();
}

Kdict::~Kdict()
{
    delete root;
}

void Kdict::insert( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    root->insert( bkmer );
    delete bkmer;
}

bool Kdict::contains( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    bool res = root->contains( bkmer );
    delete bkmer;
    return res;
}

size_t Kdict::size()
{
    return root->size();
}

void Kdict::remove( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, m_bk, kmer );
    root->remove( bkmer );
    delete bkmer;
}

void Kdict::clear()
{
    delete root;
    root = new Vertex();
}

PYBIND11_MODULE(kc, m) {
    py::class_<Kdict>(m, "Kdict")
        .def(py::init<const int>())
        .def("insert", &Kdict::insert)
        .def("contains", &Kdict::contains)
        .def("remove", &Kdict::remove)
        .def("size", &Kdict::size)
        .def("clear", &Kdict::clear)
        .def("get_root", &Kdict::get_root, py::return_value_policy::reference);
    
    py::class_<Vertex>(m, "Vertex")
        .def("get_uc", &Vertex::get_uc, py::return_value_policy::reference)
        .def("get_ccs", &Vertex::get_ccs, py::return_value_policy::reference);

    py::class_<Bkmer>(m, "Bkmer")
        .def("get_seq", &Bkmer::get_seq, py::return_value_policy::take_ownership);
    
    py::class_<UContainer>(m, "UContainer")
        .def("get_bkmers", &UContainer::get_bkmers, py::return_value_policy::reference);

    py::class_<CContainer>(m, "CContainer")
        .def("get_suf_clust_data", &CContainer::get_suf_clust_data, py::return_value_policy::reference)
        //.def_readonly("scd", &CContainer::m_suf_clust_data)
        .def("prefix_from_clust", &CContainer::prefix_from_clust, py::return_value_policy::take_ownership);

    py::class_<SufClustData>(m, "SufClustData")
        .def("get_child_vertex", &SufClustData::get_child_vertex, py::return_value_policy::reference);
    
    m.def( "calc_bk", &calc_bk, "Calculate the number of bytes needed to store a kmer." );
}


