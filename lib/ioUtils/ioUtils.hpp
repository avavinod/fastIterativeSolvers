#ifndef _IOUTILS_INCLUDE_
#define _IOUTILS_INCLUDE_

#include <vector>

int vectorOutput(std::vector<double> x,char * filename,int precision);
std::vector<double> vectorInput(char * filename);

#endif