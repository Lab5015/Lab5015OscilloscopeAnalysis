#include "interface/WFMReader.h"

#include <iostream>




int main(int argc, char** argv)
{
  std::string inFileName(argv[1]);
  
  WFMReader reader(inFileName,true);
}
