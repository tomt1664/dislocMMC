#ifndef POLYGON_H
#define POLYGON_H

#include <vector>

#include <atom.h>

class Polygon
{
    public:
        Polygon(std::vector<Atom> crds, std::vector<int> indx);
        virtual ~Polygon();
        int coords(std::vector<Atom>& crds);
        int index(std::vector<long>& indx);
    protected:
    private:
        int m_n; //the order of the polygon
        std::vector<Atom> m_coords;
        std::vector<long> m_index;
};

#endif // POLYGON_H
