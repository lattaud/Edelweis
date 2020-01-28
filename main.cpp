#include <iostream>
#include "Cryo_Run.h"
#include "Plot_HEnergy_Voltage.h"


int main(int argc, char** argv) {


	
	Plot_HEnergy_Voltage * Test_obj = new  Plot_HEnergy_Voltage(argv[1], std::atof(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), argv[5], std::stoi(argv[6]) ); 
	
	delete Test_obj ;
	std::cout<<" Ending Routine "<<std::endl;
	return 0;	 
}
