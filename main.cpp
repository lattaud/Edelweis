#include <iostream>
#include "Plot_HEnergy_Voltage.h"

//faire un vrai menu
int main(int argc, char** argv) {

//./Plot_Energy list_15_v_nbsi_fond.txt  20 0 0 TESTINGTIME 0 NbSi209 0
	TCLAP::CmdLine cmd("Skimming datafile to minimum", ' ', "0.1");
	TCLAP::ValueArg<std::string> inputfileList("", "inputList", "input files", true, "", "string",cmd);
	TCLAP::ValueArg<std::string> Detector_name("d", "detector", "Which detector", true, "", "string", cmd);
	TCLAP::ValueArg<std::string> Output_name("o", "output-dir", "Output Directory", true, "", "string", cmd);
	TCLAP::ValueArg<double> Temp("", "Temp", "Cryo temperature", true, 0, "double", cmd);
	TCLAP::ValueArg<double> IonCut("", "IonCut", "Ionization cut", true, 0, "double", cmd);
	TCLAP::ValueArg<unsigned int> Parity("", "Parity", "Which partition to process : 0 (all) , 1 (odd), 2,(pair)", true, 0, "unsigned int", cmd);
	TCLAP::SwitchArg RunOnRun("", "IsRun", "IsRun?", cmd);
	TCLAP::SwitchArg RunOnLocalList("", "IsLocalList", "IsLocalList?", cmd);
	TCLAP::SwitchArg RunOnProcessed("", "RunOnProcessed", "RunOnProcessed?", cmd);
	cmd.parse(argc, argv);
	try{
	    std::cout<<" Testing Arguments"<<std::endl;
	    if(Parity.getValue() > 2) throw 1;
	    if(Detector_name.getValue() != "RED30" && Detector_name.getValue() != "FID848" && Detector_name.getValue() != "NbSi209" ) throw 2 ;
	     
	}catch(int e){
	    if(e == 1) std::cerr<<" Parity mode not Handle "<<std::endl;
	    if(e == 2) std::cerr<<" Unknown Detector       "<<std::endl;
	    
	    return 0;		
	}
	
	std::cout<<"Launching instance with "<<inputfileList.getValue()<<" "<< Temp.getValue()<<" "<< RunOnRun.getValue()<<" "<< RunOnProcessed.getValue()<<" "<< Output_name.getValue()<<" "<< RunOnLocalList.getValue()<<" "<<  Detector_name.getValue()<<" "<<Parity.getValue()<<" Ion cut "<< IonCut.getValue()<<std::endl;
	Plot_HEnergy_Voltage * Test_obj = new  Plot_HEnergy_Voltage(inputfileList.getValue(), Temp.getValue(), RunOnRun.getValue(), RunOnProcessed.getValue(), Output_name.getValue(), RunOnLocalList.getValue(),  Detector_name.getValue(),Parity.getValue(),IonCut.getValue()); 
	
	delete Test_obj ;
	
	std::cout<<" Ending Routine "<<std::endl;
	return 0;	 
}
