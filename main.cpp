#include <iostream>
#include "Plot_HEnergy_Voltage.h"

int main(int argc, char** argv) {


	
	Plot_HEnergy_Voltage Test_obj(argv[1], std::atof(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), argv[5]) ; 
	//if(argv.size() == 0)  Test_obj.Help();
	std::cout<<" Ending Routine "<<std::endl;
	return 0;	 
}
