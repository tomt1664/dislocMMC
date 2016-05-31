#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "config.h"

// input constructor
Config::Config(std::string infile)
{
    //open data input file for inital structure
    std::ifstream initstruct("data.in");
    if (!initstruct)
    {
        std::cerr << "Error: input structure could not be opened" << std::endl;
        exit(1);
    }

    //read from file (LAMMPS data format)
    std::string data;
    std::stringstream ss; //use stringstream
    long iat = 0;
    long indx1,attype;
    double x,y,z; //atom coordinates
    while(!initstruct.eof())
    {
        std::getline(initstruct,data);

        if(iat > 0)
        {
//            std::cout << data << std::endl;
            ss << data;
            ss >> indx1 >> attype >> x >> y >> z;
            ss.str(""); //clear stringsteam
            Atom atom(x,y,z);
            m_coords.push_back(atom);
            std::cout << x << " " << y << " " << z << std::endl;
        }

        if(data.find("Atoms") !=std::string::npos)
        {
            getline(initstruct,data);  //skip next line
            iat++;
        }
    }

    std::cout << " total " << m_coords.size() << std::endl;

}

Config::~Config()
{
    //dtor
}

