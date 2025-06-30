#include "ioUtils.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

//to write vectors to files //not used anywhere
int vectorOutput(std::vector <double> x,char * filename,int precision)
{
    std::ofstream vectorfile;
    vectorfile.open (filename,std::ios::out);
    if (!vectorfile.is_open())
    {   
        std::cout<<" Error!. Vector file cannot be opened\n";
        return 0;
    }
    vectorfile.seekp (0,std::ios::beg);
    vectorfile<<x.size()<<std::endl;
    
    for (long int i=0; i<x.size(); ++i)
    {
        vectorfile<< x.at(i)<<std::setprecision(precision)<<std::scientific<<std::endl;
    }
    vectorfile.close();
    return 1;
}

//to get vectors from files //not used anywhere
std::vector<double> vectorInput(char * filename)
{
    std::ifstream vectorfile;
    vectorfile.open (filename,std::ios::in);

    if (!vectorfile.is_open())
    {   
        std::cout<<" Error!. Vector file cannot be opened\n";
        std::vector<double> err(1,0);
        return err;
    }
    vectorfile.seekg (0,std::ios::beg);

    long int size;
    vectorfile>>size;

    std::vector <double> x(size,0);

    for (long int i=0; i<x.size(); ++i)
    {
        vectorfile>> x.at(i);
    }
    vectorfile.close();

    return x;   
}