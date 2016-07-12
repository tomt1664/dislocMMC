/****************************************************
**  dislocMMC
**
**  A Metropolis Monte Carlo algorithm to model the
**  evolution of structure in defective sp2 bonded
**  carbon systems
**
**  The code calls the LAMMPS atomistic simulation
**  program for structural optimisation with a
**  reactive force-field
**
**  T.Trevethan 2016
*****************************************************/

#include "polygon.h"

Polygon::Polygon(std::vector<Atom> crds, std::vector<long> indx)
{
    //copy the coordinates and indexes to the private vectors
    int psize = crds.size();
    for(int i=0; i < psize; i++)
    {
        m_coords.push_back(crds[i]);
        m_index.push_back(indx[i]);
    }
    m_n = psize;
}

int Polygon::coords(std::vector<Atom>& crds)
{
    crds.clear();
    //copy atom coordinates
    for(int i=0; i < m_n; i++)
    {
        crds.push_back(m_coords[i]);
    }
    return m_n;
}

int Polygon::index(std::vector<long>& indx)
{
    indx.clear();
    //copy atom coordinates
    for(int i=0; i < m_n; i++)
    {
        indx.push_back(m_index[i]);
    }
    return m_n;
}

long Polygon::sum()
{
    //sum of the polygon index numbers
    long psum = 0;
    for(int i=0; i< m_n; i++)
    {
        psum += m_index[i];
    }
    return psum;
}

Polygon::~Polygon()
{
    //dtor
}
