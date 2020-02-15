#include "Cryo_Run.h"


using namespace std;


Cryo_Run::Cryo_Run(std::string Output_file){

	output_file = Output_file;
	Init();

}

Cryo_Run::~Cryo_Run(){



}

void Cryo_Run::Init(){

		g_reso_vs_time = new TGraph();
}

void Cryo_Run::fill_reso_time(Int_t npoint, Double_t reso, Double_t time){


	
	g_reso_vs_time-> SetPoint(npoint,reso,time);

}
void Cryo_Run::Write_to_file(std::string output_name){


	output = new TFile((output_file+".root").c_str(), "RECREATE");
	
	g_reso_vs_time->Write("Reso_vs_time");

}
