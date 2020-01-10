#ifndef __PLOT_HERNERGY_VOLTAGE_H__
#define __PLOT_HERNERGY_VOLTAGE_H__


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


class Plot_HEnergy_Voltage{

	private:
		TChain* chain_voltage ;
		TChain* chain_index ;		
		TChain* chain_HeatEnergy ;
		TChain* chain_chi2A;
		TChain* chain_event_processed ;
		TChain* chain_voltage_pro;
		std::string file_name;
		std::string list_name;
		std::string Run_name;
		std::string Run_name_forlist;
		std::string OutputDir;
		TFile* Input_Files;
		TFile* Output_Files;
		
		Double_t Voltage;
		Int_t    Index_Calib;
		Double_t Eh;
		Double_t chi2_A; 
		Int_t Nb_voltage;
		Int_t Nb_index;
		Int_t Nb_HeatEnergy;
		Int_t Nb_chi2;
		TH1D * H_Eh;
		TH1D * H_Eh_lowres;
		TH1D * H_Ehee;
		TH2D * H2_Eh_chi2;
		TGraphErrors * PSD_plot;
		
		
		TH1D ** PSD_spectrum; 
		TH1D * PSD_spectrum_summed;
		Int_t count_line = 0;
		Double_t heat;
		Int_t Nlist;
		std::vector<Double_t> binning_vec ; 
		std::vector<Double_t> binning_vec_low_res ; 
		std::vector<Double_t> binning_vec_kevee;
		Int_t allRUN = 0 ;
		bool IS_PROCESSED;
		Double_t PDS_noise[1024][6];
		Double_t PDS_freq[1024];
		Double_t nVtoADU[6];
		Double_t polarion[4];
		Double_t cutofffreq;
		Int_t filter_order;
		Double_t Reso_cat;
		Double_t Reso_cat_buffer;
		
		Int_t micro_step;
		Int_t N_partition;
		Double_t f_max_heat;
		
		
	public:
		Plot_HEnergy_Voltage(std::string list_name, Double_t Heat, bool IsRun, bool On_processed, std::string outputdir );
		
		~Plot_HEnergy_Voltage();
		
		void RunOnly(Double_t Heat);
		void RunList(Double_t Heat, std::string list);
		void Help();
		void Init();
		void Parse_List();
		void Open_file( std::string file_name);
		void Load_chain( std::string tree_name);
		void Loop_over_Chain();
		void Loop_over_Chain_processed();
		void Write_histo_tofile(float temp, int voltage, std::string run_name);
		void SetTemp();
		void SetTemp(Double_t heat_);
		void SetRunname();
		void SetRunname(std::string runName);
		void Estimate_Run_ellapsed_time();
		void loop_over_generic_chain(TChain* chain);
		//void GetEntryChain(Int_t Nchain, std::string  chainnames),
		void cleaning();
		Double_t EpBinIndex(Double_t Ep, std::vector<Double_t> binning );
		Double_t Kevee_weight(Double_t Eh);
		//void merge_outputfiles(float temp, int voltage);
};

#endif
