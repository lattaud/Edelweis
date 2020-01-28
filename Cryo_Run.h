#ifndef __CRYO_RUN_H__
#define __CRYO_RUN_H__


// standard libraries

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector> 
// Root related libraries
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TParameter.h"
#include "Plot_HEnergy_Voltage.h"


class Cryo_Run {


	private: 
	
	
		std::string output_file;
	
		TGraph* g_reso_vs_time;
		TFile * output;

	public:

		Cryo_Run(std::string output_file);
		~Cryo_Run();
		void fill_reso_time(Int_t npoint, Double_t reso, Double_t time);
		void Write_to_file(std::string output_name) ;
		void Init();
		friend class Plot_HEnergy_Voltage ;

};

#endif
