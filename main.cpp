#include <iostream>
#include "Plot_HEnergy_Voltage.h"

//faire un vrai menu
int main(int argc, char** argv) {

//--------------------------------------------------------Menu to launch Selection and skimmed file production -----------------------------------------------------------
	TCLAP::CmdLine cmd("Skimming datafile to minimum", ' ', "0.1");
	TCLAP::ValueArg<std::string> inputfileList("", "inputList", "List of input list from the prod/list/ directory ", true, "", "string",cmd);
	TCLAP::ValueArg<std::string> Detector_name("d", "detector", "Which detector you want to analyze. Currently supported NbSi209, RED30, FID848", true, "", "string", cmd);
	TCLAP::ValueArg<std::string> Output_name("o", "output-dir", "Output Directory to store skimmed root files", true, "", "string", cmd);
	TCLAP::ValueArg<std::string> Prod_to_use("p", "prod", "Which prod du you want to use", true, "", "string", cmd);
	TCLAP::ValueArg<double> Temp("", "Temp", "Cryo temperature (Only used for labelling )", true, 0, "double", cmd);
	TCLAP::ValueArg<double> IonCut("", "IonCut", "Cut on the Ionization Energy (negative value stand for no cut)", true, 0, "double", cmd);
	TCLAP::ValueArg<unsigned int> Parity("", "Parity", "Which partition to process : 0 (all) , 1 (odd), 2,(pair)", true, 0, "unsigned int", cmd);
	TCLAP::SwitchArg RunOnRun("", "IsRun", "Do you want to loop on one run ?  (not fully supported)", cmd);
	TCLAP::SwitchArg RunOnLocalList("", "IsLocalList", "Do you run on a custom local list ?", cmd);
	TCLAP::SwitchArg RunOnProcessed("", "RunOnProcessed", "Do you only run on processed rootfile from NEPAL ", cmd);
	TCLAP::SwitchArg Rejected_IonCut("", "RejectedIoncut", "Get event rejected by ion cut?", cmd);
	cmd.parse(argc, argv);
	try{
	    std::cout<<" Testing Arguments"<<std::endl;
	    if(Parity.getValue() > 2) throw 1;
	    if(Detector_name.getValue() != "RED30" && Detector_name.getValue() != "FID848" && Detector_name.getValue() != "NbSi209" ) throw 2 ;
	    if(Prod_to_use.getValue() != "prodj" && Prod_to_use.getValue() != "prodk" && Prod_to_use.getValue() != "prodi" ) throw 3 ;
	     
	}catch(int e){
	    if(e == 1) std::cerr<<" Parity mode not Handle "<<std::endl;
	    if(e == 2) std::cerr<<" Unknown Detector       "<<std::endl;
	    if(e == 3) std::cerr<<" prod not supported     "<<std::endl;
	    
	    return 0;		
	}
	
	std::cout<<"Launching instance with "<<inputfileList.getValue()<<" "<< Temp.getValue()<<" "<< RunOnRun.getValue()<<" "<< RunOnProcessed.getValue()<<" "<< Output_name.getValue()<<" "<< RunOnLocalList.getValue()<<" "<<  Detector_name.getValue()<<" "<<Parity.getValue()<<" Ion cut "<< IonCut.getValue()<<" Prod "<<Prod_to_use.getValue()<<std::endl;
	Plot_HEnergy_Voltage * Skimmer_instance = new  Plot_HEnergy_Voltage(inputfileList.getValue(), Temp.getValue(), RunOnRun.getValue(), RunOnProcessed.getValue(), Output_name.getValue(), RunOnLocalList.getValue(),  Detector_name.getValue(),Parity.getValue(),IonCut.getValue(), Rejected_IonCut.getValue(),Prod_to_use.getValue()); 
	
	delete Skimmer_instance ;
	
	std::cout<<" Ending Routine "<<std::endl;
	return 0;	 
}
