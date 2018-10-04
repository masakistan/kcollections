#include "Kdict.h"

Kdict::Kdict( int k ) : m_k( k )
{
    ibkmer = new Bkmer( k );
    m_bk = calc_bk( k );
    root = new Vertex();
}

Kdict::~Kdict()
{
    delete ibkmer;
    delete root;
}

void Kdict::insert( char* kmer )
{
    //Bkmer* bkmer = new Bkmer( m_k, kmer );
printf("this Value:  %p\n", ibkmer->get_bseq() );
    ibkmer->set_seq( kmer, m_k );
printf("this Value:  %p\n", ibkmer->get_bseq() );
    root->insert( ibkmer );
    //delete bkmer;
}

void Kdict::insert_bkmer( Bkmer* bkmer )
{
    root->insert( bkmer );
}

bool Kdict::contains( char* kmer )
{
    Bkmer* bkmer = new Bkmer( m_k, kmer );
    bool res = root->contains( bkmer );
    delete bkmer;
    return res;
}

bool Kdict::contains_bkmer( Bkmer* bkmer )
{
    return root->contains( bkmer );
}

size_t Kdict::size()
{
    return root->size();
}

void Kdict::remove( char* kmer )
{
    //Bkmer* bkmer = new Bkmer( m_k, kmer );
    ibkmer->set_seq( kmer, m_k );
    root->remove( ibkmer );
    //delete bkmer;
}

void Kdict::remove_bkmer( Bkmer* bkmer )
{
    root->remove( bkmer );
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
        .def("insert_bkmer", &Kdict::insert_bkmer)
        .def("contains", &Kdict::contains)
        .def("contains_bkmer", &Kdict::contains_bkmer)
        .def("remove", &Kdict::remove)
        .def("remove_bkmer", &Kdict::remove_bkmer)
        .def("size", &Kdict::size)
        .def("clear", &Kdict::clear)
        .def("get_root", &Kdict::get_root, py::return_value_policy::reference);
    
    py::class_<Vertex>(m, "Vertex")
        .def("get_uc", &Vertex::get_uc, py::return_value_policy::reference)
        .def("get_ccs", &Vertex::get_ccs, py::return_value_policy::reference);

    py::class_<Bkmer>(m, "Bkmer")
        .def(py::init<const int>())
        .def("get_seq", &Bkmer::get_seq, py::return_value_policy::take_ownership)
        .def("set_seq", &Bkmer::set_seq);
    
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


