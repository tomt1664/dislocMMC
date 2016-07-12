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

#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include <atom.h>

//polygon class to store polygons (3,4,5 and 7 membered rings)

class Polygon
{
    public:
        Polygon(std::vector<Atom> crds, std::vector<long> indx); //constructor passed arrays containing the atom coordinates and index numbers
        virtual ~Polygon();
        int coords(std::vector<Atom>& crds); //return the number of atoms and an array of the atom coordinates
        int index(std::vector<long>& indx); //return the number of atoms and a list of the atom coordinates
        int ord() { return m_n; } //return order
        long sum(); //sum of atom indices
    protected:
    private:
        int m_n; //the order of the polygon (number of atoms)
        std::vector<Atom> m_coords;
        std::vector<long> m_index;
};

#endif // POLYGON_H
