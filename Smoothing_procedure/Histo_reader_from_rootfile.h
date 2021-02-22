#ifndef __HISTO_READER_FROM_ROOTFILE__
#define __HISTO_READER_FROM_ROOTFILE__

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
#include "TH1D.h"
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
#include "tclap/CmdLine.h"
#include "TKDE.h"
#include "TParameter.h"

class Histo_reader_from_rootfile {
	public: 
	
		Histo_reader_from_rootfile( const std::string file_list, const std::string &list_histo , const std::string detector, const std::string outputname, bool verbosity,bool runontree, double voltage_user, std::string run, bool const & do_kernel, bool const & do_dru);		
		~Histo_reader_from_rootfile() = default;
		void Get_histo_from_File();
		void Open_and_store_from_list(const std::string & list, const std::string & type_obj);
		void Make_graph_from_vector(const std::string & name, std::vector<Double_t> data);
		void Init_gaussian_kernel();
		double * Fill_data_array(std::string name_variable);
		std::vector<double>  Fill_data_array_vec(std::string name_variable);
        void Optimise_smoothing  ( std::vector<double>  data);		
		double calculatechi2(TH1D *hdata,TH1* hkde,double emin,double emax);//courtoisie de Quentin
		int color_rainbow(int i,int nbval);
		void Open_inputfile(std::string const & name);
	private: 
	
		TFile *                  Input_file;
		std::string              Detector;
		Double_t                 Detector_Gemass;
		Double_t                 Time_exposition;
		std::vector<std::string> Vec_hist_name;
		std::vector<TGraph *>    Vec_paramam_vs_time_lowE ;
		std::vector<TGraph *>    Vec_paramam_vs_time_HightE ;
		std::vector<TH1D *>      Vec_Histo_to_fit;
		std::vector<std::vector<Double_t>>    Vec_parameters_lowE; 
		std::vector<std::vector<Double_t>>    Vec_parameters_hightE; 
		std::vector<double*>    Vec_Data_Hist;
		std::vector<double*>    Vec_Data_weight;
		TH1D*      Histo_from_tree;
		TF1*       Function_to_fit_lowE   ;
		TF1*       Function_to_fit_hightE ;
		std::string            Output_name; 
		std::string            RUN;
		bool IsVerbose; 
		TChain * event_chain;
		TTree  * event_tree;
		Double_t E_h;
		Double_t E_p;
		Double_t voltage;
		Double_t Voltage_user;
		Double_t weight;
		TKDE * kde_opti;
		TF1  * Analytical_bkg;
		bool Do_kernel;
		bool Do_DRU;

};


#endif
