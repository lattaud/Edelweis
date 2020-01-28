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
#include "Cryo_Run.h"

class Plot_HEnergy_Voltage {

	private:
		TChain* chain_voltage ;
		TChain* chain_index ;		
		TChain* chain_HeatEnergy ;
		TChain* chain_chi2A;
		TChain* chain_event_processed ;
		TChain* chain_voltage_pro ;
		TChain* chain_event_processed_fast ;
		
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
		Double_t Ei;
		Int_t Nb_voltage;
		Int_t Nb_index;
		Int_t Nb_HeatEnergy;
		Int_t Nb_chi2;
		
		
		
		
		TH1D * H_Eh;
		TH1D * H_Eh_lowres;
		TH1D * H_Ehee;
		TH2D * H2_Eh_chi2;
		TH1D * Time_per_voltage;
		TH2D * Ionration_vs_Ei;
		TH2D * Dchi2_vs_Ep_pass;
		TH2D * Dchi2_vs_Ep_fail;
		
		TH2D * chi2_cut_vs_Ep_pass;
		TH2D * chi2_cut_vs_Ep_fail;
		
		TGraph * reso_vs_time;
		
		TGraphErrors * PSD_plot;
		TGraphErrors * PSD_plot_reso;		
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
		bool local_list;
		
		Double_t chi2_norm;
		Double_t chi2_fast;
		Int_t point_time_reso = 0 ;
		long double Time_Crate ; 
		
	public:
		Plot_HEnergy_Voltage(const std::string &list_name , Double_t Heat , bool IsRun , bool On_processed , const std::string &outputdir, bool local_list );
		
		~Plot_HEnergy_Voltage() = default;
		
		void RunOnly(Double_t Heat);
		void RunList(Double_t Heat, const std::string &list);
		void Help();
		void Init();
		void Parse_List();
		void Open_file( const std::string &file_name);
		void Load_chain( const std::string &tree_name);
		void Loop_over_Chain();
		void Loop_over_Chain_processed();
		void Write_histo_tofile(float temp, int voltage, const std::string &run_name);
		void SetTemp();
		void SetTemp(Double_t heat_);
		void SetRunname();
		void SetRunname(const std::string &runName);
		void Estimate_Run_ellapsed_time();
		void loop_over_generic_chain( TChain* chain);
		//void GetEntryChain(Int_t Nchain, std::string  chainnames),
		void cleaning();
		Double_t EpBinIndex(Double_t Ep, std::vector<Double_t> binning );
		Double_t Kevee_weight(Double_t Eh);
		void Write_timed_reso(const std::string &name_output_reso);
		void Clean();
		//void merge_outputfiles(float temp, int voltage);
};

#endif
