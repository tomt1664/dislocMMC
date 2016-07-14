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

#ifndef ATOM_H
#define ATOM_H

//atom class: holds the atom cartesian coordinates

class Atom
{
    public:
        Atom();
        Atom(double x, double y, double z);
        virtual ~Atom();
        void setPos(double x, double y, double z);
        double getx() { return m_x; }
        double gety() { return m_y; }
        double getz() { return m_z; }

        double dist(Atom at, double bounds[5]); //method to return the distance to a second atom
    private:
        double m_x;
        double m_y;
        double m_z;
};

#endif // ATOM_H
