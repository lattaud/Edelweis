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
#include "Smoothing_procedure/tclap/CmdLine.h"

class Plot_HEnergy_Voltage {

	private:
		TChain* chain_voltage ;
		TChain* chain_index ;		
		TChain* chain_HeatEnergy ;
		TChain* chain_chi2A;
		TChain* chain_event_processed ;
		TChain* chain_voltage_pro ;
		TChain* chain_event_processed_fast ;
		TChain* chain_event_processed_Slow ;
		TChain* chain_event_processed_NTD ;
		TChain* chain_event_Reso_processed ;
		
		std::string file_name;
		std::string list_name;
		std::string Run_name;
		std::string Run_name_forlist;
		std::string OutputDir;
		std::string Detector;
		TFile* Input_Files;
		TFile* Output_Files;
		
		Double_t Voltage;
		Int_t    Index_Calib;
		Double_t Eh;
		Double_t EhA;
		Double_t EhB;
		Double_t chi2_A; 
		Double_t EiA;
		Double_t EiB;
		Int_t Nb_voltage;
		Int_t Nb_index;
		Int_t Nb_HeatEnergy;
		Int_t Nb_chi2;
		
		TTree *outTree_ ;
		TParameter <double> * time_lengh ; 
		Double_t E_h_buf;
		Double_t E_p_buf;
		Double_t weight;
		Double_t voltage;
		std::vector<Double_t> *chi2_buf;
		
		
		
		TH1D * H_Eh;
		TH1D * H_EhB;
		TH1D * H_Ehtot;
		TH1D * H_EiA;
		TH1D * H_EiB;
		TH1D * H_EiFid;
		TH1D * H_Eh_lowres;
		TH1D * H_Ehee;
		TH1D * H_EhBee;
		TH1D * H_Ehtotee;
		
		TH1D * H_Eh_noweight;
		TH1D * H_EhB_noweight;
		TH1D * H_Ehtot_noweight;
		TH1D * H_Ehee_noweight;
		TH1D * H_EhBee_noweight;
		TH1D * H_Ehtotee_noweight;
		
		TH2D * H2_Eh_chi2;
		TGraph * G2_Eh_chi2;
		TH1D * Time_per_voltage;
		TH2D * Ionration_vs_Ei;
		TH2D * Dchi2_vs_Ep_pass;
		TH2D * Dchi2_vs_Ep_fail;
		TH2D * Dchi2Slow_vs_Ep_pass;
		TH2D * Dchi2Slow_vs_Ep_fail;
		TH2D * Dchi2NTD_vs_Ep_pass;
		TH2D * Dchi2NTD_vs_Ep_fail;		
		TH2D * chi2_cut_vs_Ep_pass;
		TH2D * chi2_cut_vs_Ep_fail;
		//TH2D * chi2ion_cut_vs_Ep_fail;
		TH2D * H2_Ei_chi2;
		
				
		TGraph * reso_vs_time;
		TGraph * EiFid_vs_Eh_passcut;
		TGraph * EiFid_vs_Eh_rejected;
		
		TGraph * EiFid_vs_chi2_passcut;
		TGraph * EiFid_vs_chi2_rejected;
				
		TGraphErrors * PSD_plot;
		TGraphErrors * PSD_plot_reso;		
		TH1D ** PSD_spectrum; 
		TH1D *  PSD_spectrum_summed;
		
		
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
		Double_t chi2_norm[6];
		Double_t chi2_fast[6];
		Double_t chi2_Slow[6];
		Double_t chi2_NTD[6];
		Double_t chi2_half[6];
		Double_t Energy_OF[6];
		Double_t chi2_i[6];
		Int_t Time_Crate;
		Int_t Mega_stp;	
		Int_t N_partitiontree;	
		Int_t point_time_reso = 0 ;
		Int_t ndof_chi2;
		Double_t reso_eV ;
		Double_t weight_detector ;
		unsigned int Pair_partition;
		Double_t IonCut;
		bool Neg_cutIon;
		std::string Prod;

		
		
	public:
		Plot_HEnergy_Voltage(std::string const &list_name , Double_t const & Heat , bool const & IsRun , bool const & On_processed , std::string const & outputdir, bool const & local_list, std::string const & detector, unsigned int const & runonpairpart, Double_t const & ionCut, bool const & cut_ion_rej, std::string const &prod );		
		~Plot_HEnergy_Voltage() = default;		
		void RunOnly(Double_t const & Heat);
		void RunList(Double_t const & Heat, std::string const & list);
		void Help();
		void Init();
		void Parse_List();
		void Open_file(  std::string const  &file_name);
		void Load_chain(  std::string const  &tree_name);
		void Loop_over_Chain();
		void Loop_over_Chain_processed();
		void Write_histo_tofile(float const & temp, int const & voltage, std::string const &run_name);
		void SetTemp();
		void SetTemp(Double_t const & heat_);
		void SetRunname();
		void SetRunname( std::string const &runName);
		void Estimate_Run_ellapsed_time();
		void loop_over_generic_chain( TChain* chain);
		void cleaning();
		Double_t EpBinIndex(Double_t const & Ep, std::vector<Double_t> const & binning );
		Double_t Kevee_weight(Double_t const & Eh);
		void Write_timed_reso( std::string const & name_output_reso);
		void Clean();
		bool Pass_chi2_cut(std::string const & detector, Double_t const & Ep , Double_t const & chi2A, Double_t const & chi2B);
		bool Pass_Deltachi2_cut(std::string const & detector, Double_t const & chi2_1, Double_t const & chi2_2);
};

#endif
