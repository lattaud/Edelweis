#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>


// Root related libraries
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2D.h"
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
#include "TString.h"
#include "TGraphErrors.h"
#include "/pbs/home/h/hlattaud/tclap/CmdLine.h"
using namespace std;

TGraph* Global_eff = new TGraph();
double Eff_function(double* x, double* par);


void get_eff_from_file(std::string const& nameeffFile){

    TFile * f_efficiency = new TFile(nameeffFile.c_str(),"UPDATE");
    Global_eff = (TGraph*) f_efficiency->Get("EfficiencyRED3078Vee");
}

double Eff_function(double* x, double* par){

    Double_t Voltage = par[0];
    Double_t Ehee = x[0];
	Double_t f =  Global_eff->Eval(Ehee *1000 ) ;
	return f;

}


void Make_Nice_Hist(TH1D* hist , Double_t tension, Double_t norm, std::string xaxisname, int rebin, std::string detector_name,bool scale){

	
	if(fabs(tension) == 0 )		hist->SetLineColor(kBlack );
	if(fabs(tension )== 8 )		hist->SetLineColor(kBlue + 1);
	if(fabs(tension )== 15 )	hist->SetLineColor(kMagenta + 1);
	if(fabs(tension )== 30 )	hist->SetLineColor(kRed + 1);
	if(fabs(tension) == 60 )	hist->SetLineColor(kYellow + 1);
	if(fabs(tension )== 66 )	hist->SetLineColor(kGreen + 1);	
	if(fabs(tension )== 70 )	hist->SetLineColor(kGreen - 7);
	if(fabs(tension) == 78 )	hist->SetLineColor(kYellow - 2);
	if(fabs(tension) == 40 )	hist->SetLineColor(kCyan - 2);
	if(fabs(tension) == 51 )	hist->SetLineColor(kCyan + 3);
	if(fabs(tension) == 57 )	hist->SetLineColor(kCyan );

	
	double detector_weight = 0.03 ;
	if(detector_name == "NbSi209") detector_weight = 0.2 ;
	if(detector_name == "FID848") detector_weight = 0.87 ;
	if(scale) hist->Scale(1./(((detector_weight*norm/(3600.*24.)))));
	hist->SetLineWidth(2);
	hist->GetXaxis()->SetTitle(xaxisname.c_str());
	hist->GetXaxis()->SetTitleOffset(1.25);
	hist->GetYaxis()->SetTitle("NEvent.keV ^{-1}.day^{-1}.kg^{-1}");
	hist->GetYaxis()->SetTitleOffset(1.3);
	//hist->GetYaxis()->SetRangeUser(10,10000000);
	if (rebin == 1 ){
	
		hist->RebinX(4);
		hist->GetXaxis()->SetRangeUser(10,10000000);
	}
	hist->SetTitle("");

}
void Apply_efficiency(TH1D* hist , Double_t tension){

    TF1 * funct_eff = new TF1("eff_function",Eff_function,0.001, 1.3,1 );
    funct_eff->SetParameter(0,tension);
    for(int i = 1 ; i < hist ->GetNbinsX(); i++){
    
        double E   (hist->GetBinCenter(i));
        double DRU (hist->GetBinContent(i));
        double Eff ((funct_eff->Eval(E/(1+tension/3.)))/100.);
        double DRUcor (DRU * 1./Eff);
        //std::cout << " bin content corrected "<<DRUcor<<" bin  "<<i<<" E "<<E<<" Eff "<<funct_eff->Eval(E/(1+tension/3.))<<std::endl;
        hist->SetBinContent(i,DRUcor);
        
    }

}

void Renorm_Hist_and_store(TH1D* hist , Double_t tension, Double_t norm, std::string xaxisname, int rebin,std::string name_hist, std::string detector_name){

	
	if(fabs(tension) == 0 )		hist->SetLineColor(kBlack );
	if(fabs(tension )== 8 )		hist->SetLineColor(kBlue + 1);
	if(fabs(tension )== 15 )	hist->SetLineColor(kMagenta + 1);
	if(fabs(tension )== 30 )	hist->SetLineColor(kRed + 1);
	if(fabs(tension) == 60 )	hist->SetLineColor(kYellow + 1);
	if(fabs(tension )== 66 )	hist->SetLineColor(kGreen + 1);	
	if(fabs(tension )== 70 )	hist->SetLineColor(kGreen - 7);
	if(fabs(tension) == 78 )	hist->SetLineColor(kYellow - 2);
	if(fabs(tension) == 40 )	hist->SetLineColor(kCyan - 2);
	if(fabs(tension) == 51 )	hist->SetLineColor(kCyan + 3);
	if(fabs(tension) == 57 )	hist->SetLineColor(kCyan );

	double detector_weight = 0.03 ;
	if(detector_name == "NbSi209") detector_weight = 0.2 ;
	if(detector_name == "FID848") detector_weight = 0.87 ;
	std::cout<<"weight "<< detector_weight<<" norme time "<<norm /(3600.*24.)<<std::endl;
	hist->Scale(1./(((detector_weight*norm/(3600.*24.)))));
	hist->SetLineWidth(2);
	hist->GetXaxis()->SetTitle(xaxisname.c_str());
	hist->GetXaxis()->SetTitleOffset(1.25);
	hist->GetYaxis()->SetTitle("NEvent.keV ^{-1}.day^{-1}.kg^{-1}");
	hist->GetYaxis()->SetTitleOffset(1.3);
	//hist->GetYaxis()->SetRangeUser(10,10000000);
	if (rebin == 1 ){
	
		hist->RebinX(4);
		hist->GetXaxis()->SetRangeUser(10,10000000);
	}
	hist->SetTitle("");
	
	TFile * Ouputfile = new TFile(Form("%s_%s.root",name_hist.c_str(),detector_name.c_str()), "RECREATE");
	hist->Write(Form("%s_%s_%2.0f",name_hist.c_str(),detector_name.c_str(),tension));
	Ouputfile->Close();

}
void Fit_and_store_HO(TH1D* hist){
      
      TF1 * exp_to_fit= new TF1("HO_fit","[0]*exp(-x*[1])",0.,10);   
    //  exp_to_fit->SetParameter(0,8000);
    //  exp_to_fit->SetParameter(1,10);
      hist-> Fit(exp_to_fit,"","",0.35,0.9);
      
      TFile* output = new TFile("Output_HOFit.root","RECREATE");
      hist->Write();
      exp_to_fit->Write();
      output->Close();
      
      
}

void Make_Nice_Hist_time(TH1D* hist , Double_t tension, Double_t norm, std::string xaxisname, int rebin,int n){

	
	hist->SetLineColor(n+1);	
	hist->Scale(1./(((norm/(3600.*24.)))* 0.03));
	hist->SetLineWidth(2);
	hist->GetXaxis()->SetTitle(xaxisname.c_str());
	hist->GetXaxis()->SetTitleOffset(1.25);
	hist->GetYaxis()->SetTitle("NEvent.keV ^{-1}.day^{-1}.kg^{-1}");
	hist->GetYaxis()->SetTitleOffset(1.3);
	hist->GetXaxis()->SetRangeUser(0.08,10.);
	if (rebin == 1 ){
	
		hist->RebinX(4);
		hist->GetXaxis()->SetRangeUser(10,10000000);
	}
	hist->SetTitle("");
	//std::cout<<" Integral hist "<<hist->Integral()<<std::endl;

}

void Per_tensionHist(TH1D * hist , TH1D * hist_2 , Double_t Tension,  Double_t Tension_2){

		
	if( Tension ==  Tension_2) hist -> Add(hist_2) ;

}

void Launch_plotting(std::string  Input_file, std::string Temp, std::string List_ofrunandtension , std::string resoCAT,bool plotlimit, std::string detector_name){


	std::string bias_applied, buffer; 
	std::vector <std::string> Run_Name;
	std::vector <std::string> Run_time_hists;
	std::vector<std::string> tension_ ;
	int iter = 0 ;
	ifstream Runlistname(List_ofrunandtension.c_str(),ios::in);
   while(std::getline(Runlistname, buffer))   {    
       
       TString RUN = buffer ;
       if(RUN.Contains("#", TString::kIgnoreCase)) continue ;      
       Run_Name.push_back("Ephonon_pos"+buffer);
       Run_time_hists.push_back("Ephonon_pos"+buffer);
       tension_.push_back(buffer);
       std::cout<<" Voltage to Print " << "Ephonon_pos"+buffer<<std::endl;

	}
	
	TH1D * Hist_vec           [Run_Name.size()];
	TH1D * Hist_vec_normtime  [Run_time_hists.size()];
	TH1D * Hist_vec_kee       [Run_Name.size()];	
	Double_t Tension[Run_Name.size()];
	Double_t Norm_Tension[Run_time_hists.size()];
	std::string Tension_name;
	Double_t Sign = 0 ; 

	TLegend *leg = new TLegend(0.6, 0.6, .9, .9);
	TLegend *legneg = new TLegend(0.7, 0.7, .9, .9);
	TFile * Inputs = TFile::Open(Input_file.c_str());
	Double_t temp_tension = -999;
	TH1D * Hist_per_tension [tension_.size()];
	TH1D * temp_Hist [tension_.size()];
	Double_t who_first [tension_.size()];


	for(int i = 0; i < Run_time_hists.size()  ; i++){
		std::string hist_name = "Ephonon_pos"+tension_.at(i);	
		if(!(Inputs->Get((hist_name+"_ellapsed_time").c_str())))continue;
		Hist_vec_normtime     [i] = (TH1D*) Inputs->Get((hist_name+"_ellapsed_time").c_str());
		TString namehist = Hist_vec_normtime [i]->GetName();
		if(namehist.Contains("pos", TString::kIgnoreCase) == 1) Sign = +1;
		
		for(int j = 0 ; j < tension_.size() ; j++){
		
			if(namehist.Contains((tension_ .at(j)).c_str(), TString::kIgnoreCase) == 1) {
				
				Norm_Tension [i] = Hist_vec_normtime[i]->Integral();
				std::cout<<" tension "<<  tension_.at(j) << " ellapsed time "<<Norm_Tension [i] <<std::endl;
				//Tension [i] *= Sign ;
				continue;
				
			}
			
		}
		
	}
	
	
	int stop_it [Run_Name.size()][Run_Name.size()];
	int IS_use [Run_Name.size()] ;
	for(int i = 0; i < Run_Name.size()  ; i++){
		if(!(Inputs->Get((Run_Name.at(i)+"_tot").c_str()) ) ) continue;	
		Hist_vec     [i] = (TH1D*) Inputs->Get((Run_Name.at(i)+"_tot").c_str());
		Hist_vec_kee [i] = (TH1D*) Inputs->Get((Run_Name.at(i)+"_keVee_tot").c_str());
		Sign = -1;
		IS_use [i]= 0;
		TString namehist = Hist_vec [i]->GetName();
		if(namehist.Contains("pos", TString::kIgnoreCase) == 1) Sign = +1;
		
		for(int j = 0 ; j < tension_.size() ; j++){
		
			if(namehist.Contains((tension_ [j]).c_str(), TString::kIgnoreCase) == 1) {
				
				Tension [i] = std::stod(tension_ [j]);
				//Tension [i] *= Sign ;
				continue;
				
			}
			
		}		
		
		Renorm_Hist_and_store(Hist_vec_kee [i], Tension [i], Norm_Tension [i], " E_{heat} (keVee)", 0, resoCAT, detector_name);
		std::string prefixednamed = resoCAT + "_phonon";
		Apply_efficiency(Hist_vec [i], Tension [i]);
		Renorm_Hist_and_store(Hist_vec [i], Tension [i], Norm_Tension [i], " E_{heat} (keV)", 0, prefixednamed, detector_name);
		Make_Nice_Hist(Hist_vec_kee [i], Tension [i], Norm_Tension [i], " E_{heat} (keVee)",0, detector_name,false);
		Make_Nice_Hist(Hist_vec [i], Tension [i], Norm_Tension [i], " E_{phonon} (keV)", 0, detector_name, false);
		leg->AddEntry(Hist_vec [i],(detector_name+" "+to_string(int (Tension [i]))+"V ;"+to_string(int(Norm_Tension [i]/(3600.*24.)))+" days").c_str(),"l");
		//if(namehist.Contains("Ephonon_pos15.000000_22.mk", TString::kIgnoreCase) == 1) Fit_and_store_HO(Hist_vec [i]);
	}
	bool compare_detector = false ;

	TFile *Migdal_file = TFile::Open("Spectrum_migdal_HV.root", "UPDATE");
	
	TGraphErrors * Nbsi_graph = (TGraphErrors*) Migdal_file->Get("nbsi209_spectrum");
	Nbsi_graph->SetMarkerColor(kBlue + 1);
	Nbsi_graph->SetMarkerSize(1);
	Nbsi_graph->SetMarkerStyle(34);
	if(compare_detector)	leg->AddEntry(Nbsi_graph, "NbSi209", "p");
	
	TGraphErrors * fid842_graph = (TGraphErrors*) Migdal_file->Get("fid842_spectrum");
	fid842_graph->SetMarkerColor(kCyan + 1);
	fid842_graph->SetMarkerSize(1);
	fid842_graph->SetMarkerStyle(34);
	if(compare_detector)	leg->AddEntry(fid842_graph, "FID842", "p");
	
	TGraphErrors * fid803_graph = (TGraphErrors*) Migdal_file->Get("fid803_spectrum");
	fid803_graph->SetMarkerColor(kRed + 1);
	fid803_graph->SetMarkerSize(1);
	fid803_graph->SetMarkerStyle(34);
	if(compare_detector)	leg->AddEntry(fid803_graph, "FID803", "p");
	
	
	TGraphErrors * fid848_graph = (TGraphErrors*) Migdal_file->Get("fid848_spectrum");
	fid848_graph->SetMarkerColor(kGreen + 1);
	fid848_graph->SetMarkerSize(1);
	fid848_graph->SetMarkerStyle(34);
	if(compare_detector)	leg->AddEntry(fid848_graph, "FID848", "p");


   bool compare_migdal = false ;
   double Neventsig [4] = {0.};
   TFile *Migdal_SPEC_file = TFile::Open("Spectrum_migdal_Ee.root", "UPDATE");
   
   
   TGraphErrors * g_10Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("10_spectrum");
	g_10Mev_graph->SetMarkerColor(kYellow + 1);
	g_10Mev_graph->SetMarkerSize(1);
	g_10Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_10Mev_graph, "Migal 10 MeV 1e^{-30} cm^{2}", "p"); 
	
	TGraphErrors * g_15Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("15_spectrum");
	g_15Mev_graph->SetMarkerColor(kCyan - 8);
	g_15Mev_graph->SetMarkerSize(1);
	g_15Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_15Mev_graph, "Migal 15 MeV 1e^{-30} cm^{2}", "p"); 
	
	TGraphErrors * g_20Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("20_spectrum");
	g_20Mev_graph->SetMarkerColor(kMagenta + 2);
	g_20Mev_graph->SetMarkerSize(1);
	g_20Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_20Mev_graph, "Migal 20 MeV 1e^{-30} cm^{2}", "p"); 
   
   TGraphErrors * g_50Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("50_spectrum");
	g_50Mev_graph->SetMarkerColor(kCyan + 3);
	g_50Mev_graph->SetMarkerSize(1);
	g_50Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_50Mev_graph, "Migal 50 MeV 1e^{-30} cm^{2}", "p"); 
	
	TGraphErrors * g_100Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("100_spectrum");
	g_100Mev_graph->SetMarkerColor(kGreen + 1);
	g_100Mev_graph->SetMarkerSize(1);
	g_100Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_100Mev_graph, "Migal 100 MeV 1e^{-30} cm^{2}", "p"); 
	
	
	TGraphErrors * g_500Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("500_spectrum");
	g_500Mev_graph->SetMarkerColor(kRed + 1);
	g_500Mev_graph->SetMarkerSize(1);
	g_500Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_500Mev_graph, "Migal 500 MeV 1e^{-30} cm^{2}", "p"); 
	
	
	
	TGraphErrors * g_1000Mev_graph = (TGraphErrors*) Migdal_SPEC_file->Get("1000_spectrum");
	g_1000Mev_graph->SetMarkerColor(kCyan + 1);
	g_1000Mev_graph->SetMarkerSize(1);
	g_1000Mev_graph->SetMarkerStyle(33);
	if(compare_migdal)	leg->AddEntry(g_1000Mev_graph, "Migal 1000 MeV 1e^{-30} cm^{2}", "p");   

   TFile *STD_detector = TFile::Open("Spectrum_DMstandard_Phonon.root", "UPDATE");
   
   TGraphErrors * Nbsi_graphSTD = (TGraphErrors*) STD_detector->Get("nbsi209_spectrum");
	Nbsi_graphSTD->SetMarkerColor(kBlue + 1);
	Nbsi_graphSTD->SetMarkerSize(1);
	Nbsi_graphSTD->SetMarkerStyle(34);
	if(compare_detector)	leg->AddEntry(Nbsi_graphSTD, "NbSi209", "p");
	
	TFile *STD_sig_0_2 = TFile::Open("78V_sig/Ouputfile_0.200000.root", "UPDATE");
   
   TGraph * g_200Mev_graphSTD = (TGraph*) STD_sig_0_2->Get("g1");
	g_200Mev_graphSTD->SetMarkerColor(kBlue + 1);
	g_200Mev_graphSTD->SetMarkerSize(1);
	g_200Mev_graphSTD->SetMarkerStyle(34);
	
   
   TFile *STD_sig_0_3 = TFile::Open("78V_sig/Ouputfile_0.300000.root", "UPDATE");
   
   TGraph * g_300Mev_graphSTD = (TGraph*) STD_sig_0_3->Get("g1");
	g_300Mev_graphSTD->SetMarkerColor(kBlue + 1);
	g_300Mev_graphSTD->SetMarkerSize(1);
	g_300Mev_graphSTD->SetMarkerStyle(34);
	
	
	TFile *STD_sig_0_5 = TFile::Open("78V_sig_xsec/Ouputfile_0.500000.root", "UPDATE");
   
   TGraph * g_500Mev_graphSTD = (TGraph*) STD_sig_0_5->Get("g1");
	g_500Mev_graphSTD->SetMarkerColor(kCyan + 1);
	g_500Mev_graphSTD->SetMarkerSize(1);
	g_500Mev_graphSTD->SetMarkerStyle(34);
	
	
	
	TFile *STD_sig_1 = TFile::Open("78V_sig_xsec/Ouputfile_1.000000.root", "UPDATE");
   
   TGraph * g_1000Mev_graphSTD = (TGraph*) STD_sig_1->Get("g1");
	g_1000Mev_graphSTD->SetMarkerColor(kGreen + 1);
	g_1000Mev_graphSTD->SetMarkerSize(1);
	g_1000Mev_graphSTD->SetMarkerStyle(33);
	
	TFile *STD_sig_2 = TFile::Open("78V_sig_xsec/Ouputfile_2.000000.root", "UPDATE");
   
   TGraph * g_2000Mev_graphSTD = (TGraph*) STD_sig_2->Get("g1");
	g_2000Mev_graphSTD->SetMarkerColor(kBlue + 1);
	g_2000Mev_graphSTD->SetMarkerSize(1);
	g_2000Mev_graphSTD->SetMarkerStyle(33);
	
	TFile *STD_sig_5 = TFile::Open("78V_sig_xsec/Ouputfile_5.000000.root", "UPDATE");
   
   TGraph * g_5000Mev_graphSTD = (TGraph*) STD_sig_5->Get("g1");
	g_5000Mev_graphSTD->SetMarkerColor(kRed + 1);
	g_5000Mev_graphSTD->SetMarkerSize(1);
	g_5000Mev_graphSTD->SetMarkerStyle(33);
	
	
	TFile *STD_sig_0_5nbsi = TFile::Open("66V_sig/Ouputfile_0.500000.root", "UPDATE");
   
   TGraph * g_500Mev_graphSTDnbsi = (TGraph*) STD_sig_0_5nbsi->Get("g1");
	g_500Mev_graphSTDnbsi->SetMarkerColor(kMagenta + 2);
	g_500Mev_graphSTDnbsi->SetMarkerSize(1);
	g_500Mev_graphSTDnbsi->SetMarkerStyle(33);
	
	
	
	TFile *STD_sig_1nbsi = TFile::Open("66V_sig/Ouputfile_1.000000.root", "UPDATE");
   
   TGraph * g_1000Mev_graphSTDnbsi = (TGraph*) STD_sig_1nbsi->Get("g1");
	g_1000Mev_graphSTDnbsi->SetMarkerColor(kBlue + 1);
	g_1000Mev_graphSTDnbsi->SetMarkerSize(1);
	g_1000Mev_graphSTDnbsi->SetMarkerStyle(34);
	
	TFile *STD_sig_2nbsi = TFile::Open("66V_sig/Ouputfile_2.000000.root", "UPDATE");
   
   TGraph * g_2000Mev_graphSTDnbsi = (TGraph*) STD_sig_2nbsi->Get("g1");
	g_2000Mev_graphSTDnbsi->SetMarkerColor(kBlue + 1);
	g_2000Mev_graphSTDnbsi->SetMarkerSize(1);
	g_2000Mev_graphSTDnbsi->SetMarkerStyle(34);
	
	TFile *STD_sig_5nbsi = TFile::Open("66V_sig/Ouputfile_5.000000.root", "UPDATE");
   
   TGraph * g_5000Mev_graphSTDnbsi = (TGraph*) STD_sig_5nbsi->Get("g1");
	g_5000Mev_graphSTDnbsi->SetMarkerColor(kBlue + 1);
	g_5000Mev_graphSTDnbsi->SetMarkerSize(1);
	g_5000Mev_graphSTDnbsi->SetMarkerStyle(34);
   
 
   for(int it = 108 ; it < g_5000Mev_graphSTD->GetN() ;it ++ ){ 
break;
      Double_t x,y ;
      g_5000Mev_graphSTD->GetPoint(it,x,y);
      cout<<"signal Bin "<<it<< "  center "<<x<<" value "<<y<<endl;
   }
   
   for(int it = 0 ; it < 100 ;it ++ ){ 
break;
      Double_t x,y ;
      Nbsi_graphSTD->GetPoint(it,x,y);
      cout<<"NBSI Bin "<<it<< "  center "<<x<<" value "<<y<<endl;
   }

      std::vector<TGraphErrors*> g_sig;
      g_sig.push_back(g_10Mev_graph);
      g_sig.push_back(g_15Mev_graph);
      g_sig.push_back(g_20Mev_graph);
      g_sig.push_back(g_50Mev_graph);
      g_sig.push_back(g_100Mev_graph);
      g_sig.push_back(g_500Mev_graph);
      g_sig.push_back(g_1000Mev_graph);
      
     
      
      std::vector<TGraph*> g_sigSTD;


     // g_sigSTD.push_back(g_200Mev_graphSTD);
     // g_sigSTD.push_back(g_300Mev_graphSTD);
      g_sigSTD.push_back(g_500Mev_graphSTD);
      g_sigSTD.push_back(g_1000Mev_graphSTD);
      g_sigSTD.push_back(g_2000Mev_graphSTD);
      g_sigSTD.push_back(g_5000Mev_graphSTD);
      
      
      std::vector<TGraph*> g_sigSTDnbsi;


     // g_sigSTD.push_back(g_200Mev_graphSTD);
     // g_sigSTD.push_back(g_300Mev_graphSTD);
      g_sigSTDnbsi.push_back(g_500Mev_graphSTDnbsi);
      g_sigSTDnbsi.push_back(g_1000Mev_graphSTDnbsi);
      g_sigSTDnbsi.push_back(g_2000Mev_graphSTDnbsi);
      g_sigSTDnbsi.push_back(g_5000Mev_graphSTDnbsi);
      double mass_array[7] = {10,15,20,50,100,500,1000};
      double mass_arraySTD[4] = {500,1000,2000,5000};
     // double crosssec[4] = {8.49e-32,2.19e-32,2.05e-33,6.05e-34};
     TGraph * Limit = new TGraph();
     TGraph * Limit_nbsi = new TGraph();
     TGraph * LimitSTD = new TGraph();
     TGraph * Limit_nbsiSTD = new TGraph();
     std::vector<std::multimap<double, Int_t>> vect_Ee_bin_sig;
     std::vector<std::multimap<double, Int_t>> vect_Ee_bin_sigSTD;
     std::vector<std::multimap<double, Int_t>> vect_Ee_bin_sigSTDnbsi;
     std::multimap<double, Int_t> Ee_bin_bkg;
     std::multimap<double, Int_t> Ee_bin_nbsi;
     
     std::multimap<double, Int_t> Ee_bin_bkg_phonon;
     std::multimap<double, Int_t> Ee_bin_nbsi_phonon;
     if(plotlimit) {
     for(int it_dat = 0 ;it_dat < Hist_vec_kee [0]-> GetNbinsX() ; it_dat++ ){
     
     
      Ee_bin_bkg.insert(std::make_pair(Hist_vec_kee [0]->GetBinCenter(it_dat),it_dat));
     
     
     }
     
     for(int it_dat = 0 ;it_dat < Hist_vec [0]-> GetNbinsX() ; it_dat++ ){
     
     
     
      Ee_bin_bkg_phonon.insert(std::make_pair(Hist_vec [0]->GetBinCenter(it_dat),it_dat));
     
     
     }
     for(int it_dat = 0 ;it_dat < Nbsi_graph->GetN() ; it_dat++ ){
            double x,y;
            Nbsi_graph->GetPoint(it_dat,x,y);
            Ee_bin_nbsi.insert(std::make_pair(x,it_dat));
          
     }
     
     for(int it_dat = 0 ;it_dat < Nbsi_graphSTD->GetN() ; it_dat++ ){
            double x,y;
            Nbsi_graphSTD->GetPoint(it_dat,x,y);
            Ee_bin_nbsi_phonon.insert(std::make_pair(x,it_dat));
          
     }
     
     
     for(int nsig = 0; nsig <7 ; nsig++){
     
        std::multimap<double, Int_t> Ee_bin_sig;
        
         for(int it_dat = 0 ;it_dat < g_sig.at(nsig)->GetN() ; it_dat++ ){
            double x,y;
            g_sig.at(nsig)->GetPoint(it_dat,x,y);
            Ee_bin_sig.insert(std::make_pair(x,it_dat));
          
        }  
         vect_Ee_bin_sig.push_back(Ee_bin_sig) ;
      }
      
      for(int nsig = 0; nsig <4; nsig++){
     
        std::multimap<double, Int_t> Ee_bin_sig;
        
         for(int it_dat = 0 ;it_dat < g_sigSTD.at(nsig)->GetN() ; it_dat++ ){
            double x,y;
            g_sigSTD.at(nsig)->GetPoint(it_dat,x,y);
            Ee_bin_sig.insert(std::make_pair(x,it_dat));
          
        }  
         vect_Ee_bin_sigSTD.push_back(Ee_bin_sig) ;
      }
      
      for(int nsig = 0; nsig <4; nsig++){
     
        std::multimap<double, Int_t> Ee_bin_sig;
        
         for(int it_dat = 0 ;it_dat < g_sigSTDnbsi.at(nsig)->GetN() ; it_dat++ ){
            double x,y;
            g_sigSTDnbsi.at(nsig)->GetPoint(it_dat,x,y);
            Ee_bin_sig.insert(std::make_pair(x,it_dat));
          
        }  
         vect_Ee_bin_sigSTDnbsi.push_back(Ee_bin_sig) ;
      }
      
      
      // limit Migdal
      for(int nsig = 1; nsig <7 ; nsig++){
        if(!plotlimit) break;

           bool isbest_max =false;
           bool isbest_min =false;
           double E_min = 0.035;
           
           int iterator = 0 ;
           double limit = 999.;
           while(!isbest_min  ){
              double limit_out = 999;
              double E_max = E_min + 0.01;
              int iterator_2 = 0 ;
              isbest_max =false;
              while(!isbest_max  ){
               
               Int_t bin_up_sig  = vect_Ee_bin_sig.at(nsig).lower_bound(E_max)->second;         
               Int_t bin_down_sig = vect_Ee_bin_sig.at(nsig).lower_bound(E_min)->second  ;
               Int_t bin_up_dat = Ee_bin_bkg.lower_bound(E_max)->second;         
               Int_t bin_down_dat = Ee_bin_bkg.lower_bound(E_min)->second ;
           
               
               double Nsignal = g_sig.at(nsig)->Integral(bin_down_sig,bin_up_sig);
               double NEvent_data = Hist_vec_kee [0]->Integral(bin_down_dat,bin_up_dat,"WIDTH") * 1./0.6 ;
                         
               double Npoissonup =NEvent_data + 1.28*sqrt(NEvent_data);
               
               if(NEvent_data < 2.3 && NEvent_data > 0.1){
                  Npoissonup = 3.89 ;
               
               }else if(NEvent_data < 0.1 ){
               
                  Npoissonup = 2.3 ;
               }
               double limit_in = Npoissonup * 1e-30/ Nsignal;
               
               if(limit_in < limit_out ) {
                 std::cout<<" tested  Emin"<<E_min<<" Emax "<<E_max<<" corresponding bin data "<<bin_up_dat<<" "<< bin_down_dat<<" sig "<<bin_up_sig<<" "<<bin_down_sig<<std::endl; 
                  cout<<"   Signal integral Nevent "<<mass_array[nsig]<<"MeV  "<<Nsignal<<"  NEvent Data "<<NEvent_data<<endl;
                  cout<<"   Signal           limit "<<mass_array[nsig]<<"MeV  "<<limit_in<<" out "<<limit_out<<" stock "<<limit<<" Sensibility "<<Nsignal/ (NEvent_data+ Nsignal)<<endl;
                  cout<< " "<<endl;
                  limit_out = limit_in ;
                  iterator_2 ++;                  
                  E_max += 0.01 ;
               }else if( limit_out <= limit_in  && E_max < 0.95 ) {
                  // cout<<" Emax "<< E_max << endl;
                  E_max += 0.01 ;
               }else {
                     isbest_max = true ;
                     break;
               }
              
               if( Nsignal/ (NEvent_data+ Nsignal) == 0 ) break;
            }
            if(limit_out < limit ) {
                  limit = limit_out ;
                  iterator++;
                 // cout<<"   Inf bound           limit "<<mass_array[nsig]<<"MeV  "<<limit_out<<" out "<<limit<<endl;
                  cout<< " "<<endl;
                  E_min += 0.01 ;
               }else if( limit <= limit_out && E_min < 0.95 ) {
                   
                  E_min += 0.01 ;
               }else  {
                  isbest_min = true ;
                  break;

               } 
            
          }
         std::cout<<mass_array[nsig]<<"MeV Migdal ROI optimised limits " <<  limit << std::endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         Limit->SetPoint(nsig,mass_array[nsig], limit )  ; 
         
         
      }
     for(int nsig = 0; nsig <7 ; nsig++){
     if(!plotlimit) break;
         bool isbest_max =false;
         bool  isbest_min =false;
         double E_min = 0.1; // safe 100 eV threshold efficiencies.
           
         int iterator = 0 ;
         double limit = 999.;
           while(!isbest_min){
              double limit_out = 999;
              double E_max = E_min + 0.01;
              int iterator_2 = 0 ;
              isbest_max =false;
              while(!isbest_max){
           
               Int_t bin_up_sig  = vect_Ee_bin_sig.at(nsig).lower_bound(E_max)->second;         
               Int_t bin_down_sig = vect_Ee_bin_sig.at(nsig).lower_bound(E_min)->second  ;
               Int_t bin_up_dat = Ee_bin_nbsi.lower_bound(E_max)->second;         
               Int_t bin_down_dat = Ee_bin_nbsi.lower_bound(E_min)->second ;
           
              
               double Nsignal = g_sig.at(nsig)->Integral(bin_down_sig,bin_up_sig+1);
               double NEvent_data = Nbsi_graph->Integral(bin_down_dat,bin_up_dat+1) * 1./0.6 ;
              
             //            
               double Npoissonup =NEvent_data + 1.28*sqrt(NEvent_data);
                if(NEvent_data < 2.3 && NEvent_data > 0.1){
                  Npoissonup = 3.89 ;
               
               }else if(NEvent_data < 0.1 ){
               
                  Npoissonup = 2.3 ;
               }
               double limit_in = Npoissonup * 1e-30/ Nsignal;
               if( Nsignal/ (NEvent_data+ Nsignal) == 0 ) limit_in = 999;
               
               
               
               if(limit_in < limit_out ) {
                  limit_out = limit_in ;
                  iterator_2 ++;
                   std::cout<<" tested  Emin"<<E_min<<" Emax "<<E_max<<" corresponding bin data "<<bin_up_dat<<" "<< bin_down_dat<<" sig "<<bin_up_sig<<" "<<bin_down_sig<<std::endl; 
                   cout<<"   Signal integral Nevent "<<mass_array[nsig]<<"MeV  "<<Nsignal<<"  NEvent Data "<<Npoissonup<<endl;
                   cout<<"   Signal           limit "<<mass_array[nsig]<<"MeV  "<<limit_in<<" Sensibility "<<Nsignal/ (NEvent_data+ Nsignal)<<endl;
                   cout<< " "<<endl;
                   E_max+= 0.01 ;
               }else if( limit_out <= limit_in &&  E_max < 0.95) {

                  E_max+= 0.01 ;
                  
               }else {
                     isbest_max = true ;
                     break;
               }
               
               
            }
            if(limit_out <= limit ) {
                  limit = limit_out ;
                  iterator++;
                  E_min += 0.01 ;
               }else if( limit <= limit_out &&  E_min < 0.95) {
                   
                  E_min+= 0.01 ;
               }else {
                  isbest_min = true ;
                  break;
               } 
            if(E_min > 0.95) break;
            
          }
         std::cout<<mass_array[nsig]<<"MeV Migdal NBSI ROI optimised limits " <<  limit << std::endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         Limit_nbsi->SetPoint(nsig,mass_array[nsig], limit )  ; 
     }
     
     
     
     // limit standard WIMP search 
     
     for(int nsig = 1; nsig < 4 ; nsig++){
        
      if(!plotlimit)break;
           bool isbest_max =false;
           bool isbest_min =false;
           double E_min = 0.300;
           
           int iterator = 0 ;
           double limit = 999.;
           while(!isbest_min){
              double limit_out = 999;
              double E_max = E_min + 0.040;
              int iterator_2 = 0 ;
              isbest_max =false;
              while(!isbest_max){
           
               Int_t bin_up_sig  = vect_Ee_bin_sigSTD.at(nsig).lower_bound(E_max*1000.)->second;         
               Int_t bin_down_sig = vect_Ee_bin_sigSTD.at(nsig).lower_bound(E_min*1000.)->second  ;
               Int_t bin_up_dat = Ee_bin_bkg_phonon.lower_bound(E_max)->second;         
               Int_t bin_down_dat = Ee_bin_bkg_phonon.lower_bound(E_min)->second ;
           
             //  
               double Nsignal = g_sigSTD.at(nsig)->Integral(bin_down_sig,bin_up_sig);
               double NEvent_data = Hist_vec [0]->Integral(bin_down_dat,bin_up_dat,"WIDTH") * 1./0.6 ;
             //  cout<<"   Signal integral Nevent "<<mass_arraySTD[nsig]<<"MeV  "<<Nsignal<<"  NEvent Data "<<NEvent_data<<endl;          
                double Npoissonup =NEvent_data + 1.28*sqrt(NEvent_data);
                if(NEvent_data < 2.3 && NEvent_data > 0.1){
                  Npoissonup = 3.89 ;
               
               }else if(NEvent_data < 0.1 ){
               
                  Npoissonup = 2.3 ;
               }
               double limit_in = Npoissonup * 1e-37/ Nsignal;
             //  cout<<"   Signal           limit "<<mass_arraySTD[nsig]<<"MeV  "<<limit_in<<" Sensibility "<<Nsignal/ (NEvent_data+ Nsignal)<<endl;
               if( Nsignal/ (NEvent_data+ Nsignal) == 0 ) limit_in = 999;
               if(limit_in < limit_out ) {
                  limit_out = limit_in ;
                  iterator_2 ++;
                  std::cout<<" tested  Emin"<<E_min<<" Emax "<<E_max<<" corresponding bin data "<<bin_up_dat<<" "<< bin_down_dat<<" sig "<<bin_up_sig<<" "<<bin_down_sig<<std::endl; 
                  cout<<"   Signal integral Nevent "<<mass_arraySTD[nsig]<<"MeV  "<<Nsignal<<"  NEvent Data "<<Npoissonup<<endl;
                  cout<<"   Signal           limit "<<mass_arraySTD[nsig]<<"MeV  "<<limit_in<<" Sensibility "<<Nsignal/ (NEvent_data+ Nsignal)<<endl;
                  cout<< " "<<endl;
                  E_max+= 0.04 ;                  
               }else if( limit_out <= limit_in && E_max < 10.) {
                   
                  E_max+= 0.04 ;
               }else {
                     isbest_max = true ;
                     break;
               }


            }
          //  cout<<" limit out "<< limit_out << endl;
            if(limit_out < limit ) {
                  limit = limit_out ;
                  iterator++;
                  E_min += 0.04 ;
               }else if(limit < limit_out && E_min < 10.){
                  E_min += 0.04 ;
               }else {
                  isbest_min = true ;
                  break;
                 // std::cout<<" ROI  Emin"<<E_min<<" Emax "<<E_max<<std::endl; 
               } 
            
          }
         std::cout<<mass_arraySTD[nsig]<<"MeV Standard ROI optimised limits " <<  limit << std::endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         LimitSTD->SetPoint(nsig,mass_arraySTD[nsig], limit )  ; 
         
         
      }
     for(int nsig = 0; nsig <4 ; nsig++){
      if(!plotlimit)break;
         bool isbest_max =false;
         bool  isbest_min =false;
         double E_min = 1.; // safe 100 eV threshold efficiencies.
           
         int iterator = 0 ;
         double limit = 999.;
           while(!isbest_min ){
              double limit_out = 999;
              double E_max = E_min + 0.120;
              int iterator_2 = 0 ;
              isbest_max =false;
              while(!isbest_max ){
           
               Int_t bin_up_sig  = vect_Ee_bin_sigSTDnbsi.at(nsig).lower_bound(E_max*1000.)->second;         
               Int_t bin_down_sig = vect_Ee_bin_sigSTDnbsi.at(nsig).lower_bound(E_min*1000.)->second  ;
               Int_t bin_up_dat = Ee_bin_nbsi_phonon.lower_bound(E_max)->second;         
               Int_t bin_down_dat = Ee_bin_nbsi_phonon.lower_bound(E_min)->second ;
           
              // std::cout<<" tested  Emin"<<E_min<<" Emax "<<E_max<<" corresponding bin data "<<bin_up_dat<<" "<< bin_down_dat<<" sig "<<bin_up_sig<<" "<<bin_down_sig<<std::endl; 
               double Nsignal = g_sigSTDnbsi.at(nsig)->Integral(bin_down_sig,bin_up_sig+1);
               double NEvent_data = Nbsi_graphSTD->Integral(bin_down_dat,bin_up_dat+1) * 1./0.6 ;
               if(NEvent_data < 1. && NEvent_data > 0.1){
                  NEvent_data = 3.89 ;
               
               }else if(NEvent_data < 0.1 ){
               
                  NEvent_data = 2.8 ;
               }
               //cout<<"   Signal integral Nevent "<<mass_arraySTD[nsig]<<"MeV  "<<Nsignal<<"  NEvent Data "<<NEvent_data<<endl;          
                double Npoissonup =NEvent_data + 1.28*sqrt(NEvent_data);
                if(NEvent_data < 2.3 && NEvent_data > 0.1){
                  Npoissonup = 3.89 ;
               
               }else if(NEvent_data < 0.1 ){
               
                  Npoissonup = 2.3 ;
               }
               double limit_in = Npoissonup * 1e-40/ Nsignal;
               if( Nsignal/ (NEvent_data+ Nsignal) == 0 ) limit_in = 999;
               
               if(limit_in < limit_out ) {
                  limit_out = limit_in ;
                  iterator_2 ++;
                  std::cout<<" tested  Emin"<<E_min<<" Emax "<<E_max<<" corresponding bin data "<<bin_up_dat<<" "<< bin_down_dat<<" sig "<<bin_up_sig<<" "<<bin_down_sig<<std::endl; 
                  cout<<"   Signal integral Nevent "<<mass_arraySTD[nsig]<<"MeV  "<<Nsignal<<"  NEvent Data "<<Npoissonup<<endl;
                  cout<<"   Signal           limit "<<mass_arraySTD[nsig]<<"MeV  "<<limit_in<<" Sensibility "<<Nsignal/ (NEvent_data+ Nsignal)<<endl;
                  cout<< " "<<endl;
                  E_max+= 0.120 ;
               }else if( limit_out <= limit_in && E_max < 10.) {
                   
                  E_max+= 0.012 ;
               }else {
                     isbest_max = true ;
                     break;
               }
               if(Nsignal/ (NEvent_data+ Nsignal) == 0 ) break;
               
            }
          //  cout<<" limit out "<< limit_out << endl;
            if(limit_out < limit ) {
                  limit = limit_out ;
                  iterator++;
                  E_min += 0.120 ;
               }else if(limit < limit_out && E_min < 10.){
               
                E_min += 0.120 ;
               }else  {
                  isbest_min = true ;
                  break;
                 // std::cout<<" ROI  Emin"<<E_min<<" Emax "<<E_max<<std::endl;

               } 
            
          }
         std::cout<<mass_arraySTD[nsig]<<"MeV Standard NBSI ROI optimised limits " <<  limit << std::endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         cout<< " "<<endl;
         Limit_nbsiSTD->SetPoint(nsig,mass_arraySTD[nsig], limit )  ; 
     }
     
     for(int nsig = 0; nsig <7 ; nsig++){
            if(!plotlimit) break;
           double Nsignal = g_sig.at(nsig)->Integral(144,199);
           double NEvent_data = Hist_vec_kee [0]->Integral(8,105,"WIDTH") * 1./0.6 ;
          // cout<<"  Signal integral Nevent "<<nsig<<"  " <<Nsignal<<endl;
           double SUM_SIG = 0 ;           
           
         for(int it = 144 ; it < 199 ;  it ++){
           double x,y,x2,y2,x3,y3;
           g_sig.at(nsig)->GetPoint(it+1,x,y);
           g_sig.at(nsig)->GetPoint(it,x2,y2);
           double Sig_step = x- x2 ;
            g_sig.at(nsig)->GetPoint(it,x3,y3);
            SUM_SIG += y3 * Sig_step;
         }
     //    cout<<"  Signal Sum Nevent "<<nsig<<"  " <<SUM_SIG<<" all "<<g_sig.at(nsig)->Integral(144,199)<<endl;
         
         
         double Npoissonup =NEvent_data + 1.28*sqrt(NEvent_data);
         double limit = Npoissonup * 1e-30/ Nsignal;
         cout<<"   Signal  Standard    limit "<<mass_array[nsig]<<"MeV  "<<limit<<" Sensibility "<<Nsignal/ (NEvent_data+ Nsignal)<<endl;
         //g_sig.at(nsig)->MovePoints(0,limit);
         Limit->SetPoint(nsig,mass_array[nsig], limit )  ; 
      }
     // cout<<"Data Nevent "<<Hist_vec_kee [0]->Integral(8,105,"WIDTH")<<endl;
      
      Limit->SetLineWidth(1.5);
      Limit->SetLineColor(kRed+2);
      Limit->SetMarkerStyle(33);
      Limit->SetMarkerSize(2);
      Limit->SetMarkerColor(kRed+2);
      Limit_nbsi->SetLineWidth(1.5);
      Limit_nbsi->SetLineColor(kGreen+2);
      Limit_nbsi->SetMarkerStyle(33);
      Limit_nbsi->SetMarkerSize(2);
      Limit_nbsi->SetMarkerColor(kGreen+2);
      
      LimitSTD->SetLineWidth(1.5);
      LimitSTD->SetLineColor(kRed+2);
      LimitSTD->SetMarkerStyle(34);
      LimitSTD->SetMarkerSize(2);
      LimitSTD->SetMarkerColor(kRed+2);
      Limit_nbsiSTD->SetLineWidth(1.5);
      Limit_nbsiSTD->SetLineColor(kGreen+2);
      Limit_nbsiSTD->SetMarkerStyle(34);
      Limit_nbsiSTD->SetMarkerSize(2);
      Limit_nbsiSTD->SetMarkerColor(kGreen+2);
      }//end limit computation
      std::vector<TGraph*> g_sigSTD_keV;


     // g_sigSTD.push_back(g_200Mev_graphSTD);
     // g_sigSTD.push_back(g_300Mev_graphSTD);
      TGraph * g_500Mev_graphSTDkeV = new TGraph();
      TGraph * g_1000Mev_graphSTDkeV = new TGraph();
      TGraph * g_2000Mev_graphSTDkeV = new TGraph();
      TGraph * g_5000Mev_graphSTDkeV = new TGraph();
      g_sigSTD_keV.push_back(g_500Mev_graphSTD);
      g_sigSTD_keV.push_back(g_1000Mev_graphSTD);
      g_sigSTD_keV.push_back(g_2000Mev_graphSTD);
      g_sigSTD_keV.push_back(g_5000Mev_graphSTD);
      
      for(int ite = 0 ; ite < g_sigSTD.size() ; ite ++){
         for(int it = 0 ; it < g_sigSTD.at(ite)->GetN() ;it ++ ){ 

         Double_t x,y,xlim = 0 ,ylim  =0 ;
         g_sigSTD.at(ite)->GetPoint(it,x,y);

         g_sigSTD_keV.at(ite)->SetPoint(it, x/1000., y * 1000. ) ;
         
         g_sigSTD_keV.at(ite)->SetMarkerColor(kRed + ite );
	      g_sigSTD_keV.at(ite)->SetMarkerSize(1);
	      g_sigSTD_keV.at(ite)->SetMarkerStyle(33);
	    //  std::cout<<ite<< " Energy "<<x/1000. <<" Amplitude "<<y*1000.<<std::endl;

     }
     std::cout<<" "<<std::endl;
      std::cout<<" "<<std::endl;
       std::cout<<" "<<std::endl;
      }
      

	TFile *red20_file = TFile::Open("output_red20.root", "UPDATE");
	
	
	TGraphErrors * Red_20_graph = (TGraphErrors*) red20_file->Get("Red20_results");
	Red_20_graph->SetMarkerColor(kBlack);
	Red_20_graph->SetMarkerSize(1);
	Red_20_graph->SetMarkerStyle(33);
	//leg->AddEntry(Red_20_graph, "RED20 surface spectrum", "p");
	
	
	TFile *DM_file = TFile::Open("DM_signal/Root_Files/DM_from_graph_V2.root", "UPDATE");
	
	
	TGraphErrors * DMlim = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_3e-38cm2_66.000000V");
	DMlim->SetMarkerColor(kBlue);
	DMlim->SetMarkerSize(1);
	DMlim->SetMarkerStyle(33);
	if(plotlimit)leg->AddEntry(DMlim, "WIMP 2 GeV 3e^{-38} cm ^{2}", "p");
	
	TGraphErrors * DMlim_2 = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_1e-38cm2_66.000000V");
	DMlim_2->SetMarkerColor(kBlue);
	DMlim_2->SetMarkerSize(0.7);
	DMlim_2->SetMarkerStyle(8);
	if(plotlimit)leg->AddEntry(DMlim_2, "WIMP 2 GeV 1e^{-38} cm ^{2}", "p");
	
	/*TGraphErrors * DMlim_3 = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_1e-37cm2_0.000000V");
	DMlim_3->SetMarkerColor(kBlue);
	DMlim_3->SetMarkerSize(0.7);
	DMlim_3->SetMarkerStyle(22);
	leg->AddEntry(DMlim_3, "WIMP 2 GeV 1e^{-37} cm ^{2}", "p");*/
	
	
	TGraphErrors * DMlim_test = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_4e-37cm2_0.000000V");
	DMlim_test->SetMarkerColor(kBlue);
	DMlim_test->SetMarkerSize(0.7);
	DMlim_test->SetMarkerStyle(22);
	if(plotlimit)leg->AddEntry(DMlim_test, "WIMP 2 GeV 4e^{-37} cm ^{2} (RED20)", "p");
	
	TGraphErrors * DMlim_4 = (TGraphErrors*) DM_file->Get("DM_Ephonon_1GeV_1e-37cm2_66.000000V");
	DMlim_4->SetMarkerColor(kRed);
	DMlim_4->SetMarkerSize(0.7);
	DMlim_4->SetMarkerStyle(33);
	if(plotlimit)leg->AddEntry(DMlim_4, "WIMP 1 GeV 1e^{-37} cm ^{2}", "p");
	
	TGraphErrors * DMlim_5 = (TGraphErrors*) DM_file->Get("DM_Ephonon_1GeV_1e-36cm2_66.000000V");
	DMlim_5->SetMarkerColor(kRed);
	DMlim_5->SetMarkerSize(0.7);
	DMlim_5->SetMarkerStyle(22);
	if(plotlimit)leg->AddEntry(DMlim_5, "WIMP 1 GeV 1e^{-36} cm ^{2}", "p");
	
	TCanvas *c = new TCanvas("c","c",800,800);    
	c->SetLogy();
	c->SetLogx();
	c->SetGridx();
	c->SetGridy();
	
	TH2D * axes = new TH2D("axes", "", 100, 0.1, 12., 100, 0.1, 10000000  );
	axes->GetXaxis()->SetTitle(" E_{phonon} (keV)");
	axes->GetXaxis()->SetTitleOffset(1.25);
	axes->GetYaxis()->SetTitle("NEvent.keV ^{-1}.day^{-1}.kg^{-1}");
	axes->GetYaxis()->SetTitleOffset(1.3);
	axes->GetYaxis()->SetRangeUser(1,10000000);
	
	axes->Draw("");
	axes->SetStats(kFALSE);
	//leg->AddEntry(Nbsi_graphSTD, "NbSi209", "p"); 
	
	//leg->AddEntry(g_sigSTD_keV.at(0), "DM 0.5 GeV #sigma = 1e^{-37}  cm^{2}", "p"); 
	//leg->AddEntry(g_sigSTD_keV.at(1), "DM 1 GeV #sigma = 1e^{-37}  cm^{2}", "p"); 
	//leg->AddEntry(g_sigSTD_keV.at(2), "DM 2 GeV #sigma = 1e^{-37}  cm^{2}", "p"); 
	//leg->AddEntry(g_sigSTD_keV.at(3), "DM 5 GeV #sigma = 1e^{-37}  cm^{2}", "p"); 
	
	//Red_20_graph->Draw("PSAME");
   //Nbsi_graphSTD->Draw("PSAME");
	
	bool plot_other = false ;
	
	 for(int nsig = 0; nsig <g_sigSTD_keV.size(); nsig++){

            //g_sigSTD_keV.at(nsig)->Draw("PSAME");
            
   }
	
	if(plotlimit){
	       // DMlim       ->Draw("PSAME");
	       // DMlim_2     ->Draw("PSAME");
	        //DMlim_3     ->Draw("PSAME");
	        //DMlim_4     ->Draw("PSAME");
	        //DMlim_5     ->Draw("PSAME");
	        //DMlim_test  ->Draw("PSAME");
	}
	for(int i = 0; i < Run_Name.size()  ; i++){
      if(!(Inputs->Get((Run_Name.at(i)+"_tot").c_str()) ) ) continue;
		Hist_vec [i]->SetStats(kFALSE);
		Hist_vec [i]->Draw("HIST E SAME");
		
		
	}
	leg->Draw();
	std::string name_plot = "Plot_output/Eh_vs_V_"+Temp+"mk_normalized"+resoCAT+".pdf" ;
	c->SaveAs(name_plot.c_str());

	
	
	TCanvas *cc = new TCanvas("cc","cc",1000,1000);    
	cc->SetLogy();
	cc->SetLogx();
	
	cc->SetGridy();
	cc->SetGridx();
	
	TH2D* axes2 = new TH2D("axes2", "", 100, 0.01, 2., 100, 1., 100000000  );
	axes2->GetXaxis()->SetTitle(" E_{heat} (keVee)");
	axes2->GetXaxis()->SetTitleOffset(1.25);

	axes2->GetYaxis()->SetTitle("NEvent.keVee ^{-1}.day^{-1}.kg^{-1}");
	axes2->GetYaxis()->SetTitleOffset(1.3);
	axes2->GetYaxis()->SetRangeUser(1,1000000000);
	
	axes2->Draw("");
	axes2->SetStats(kFALSE);
	
	/*
	Red_20_graph->Draw("PSAME");
	g_10Mev_graph->Draw("PSAME");
	g_15Mev_graph->Draw("PSAME");
	g_20Mev_graph->Draw("PSAME");
	g_50Mev_graph->Draw("PSAME");
	g_100Mev_graph->Draw("PSAME");
	g_500Mev_graph->Draw("PSAME");
	g_1000Mev_graph->Draw("PSAME");*/
	
	if(plot_other == true){
      Nbsi_graph->Draw("PSAME");
   	Nbsi_graph->Draw("PSAME");
   	fid842_graph->Draw("PSAME");
	   fid803_graph->Draw("PSAME");
	   fid848_graph	->Draw("PSAME");
	}
	for(int i = 0; i < Run_Name.size()  ; i++){

      if(!(Inputs->Get((Run_Name.at(i)+"_tot").c_str()) ) ) continue;
		Hist_vec_kee [i]->SetStats(kFALSE);
		Hist_vec_kee [i]->Draw("HIST E SAME");
	}
	leg->Draw();
	std::string name_plotneg = "Plot_output/Eh_kee_vs_V_"+Temp+"mk_normalized"+resoCAT+".pdf" ;
	cc->SaveAs(name_plotneg.c_str());
   TLegend *leglim = new TLegend(0.6, 0.8, .9, .9);
   //leglim->AddEntry(Limit, "RED30 Migdal", "p l");
   //leglim->AddEntry(Limit_nbsi, "NbSi209 Migdal", "p l");
   
   leglim->AddEntry(LimitSTD, "RED30 Standard", "p l");
   leglim->AddEntry(Limit_nbsiSTD, "NbSi209 Standard", "p l");
	axes2 = new TH2D("axes3", "", 100, 10, 6000, 100, 1e-41, 1e-27  );
	axes2->GetXaxis()->SetTitle(" Mass DM _{N} (MeV)");
	axes2->GetXaxis()->SetTitleOffset(1.25);
	axes2->GetYaxis()->SetTitle("#sigma_{DM_{N}} (cm^{2)}");
	axes2->GetYaxis()->SetTitleOffset(1.8);
	axes2->GetYaxis()->SetLabelSize(0.025);
	//axes2->GetYaxis()->SetRangeUser(1,1000000000);
   axes2->GetYaxis()->SetTitleSize(.025);
	axes2->Draw("");
	axes2->SetStats(kFALSE);
   /*Limit->Draw("PSAME ");
   Limit_nbsi->Draw("PSAME ");
   LimitSTD->Draw("PSAME ");
   Limit_nbsiSTD->Draw("PSAME ");*/
   
   
   leglim->Draw();
   //cc->SaveAs(("Limits_DM_"+resoCAT+".pdf").c_str());
}

void Launch_plotting_list(std::string  Input_list, std::string Temp, std::string List_ofrunandtension , std::string resoCAT,bool plotlimit){

        std::string name_plot = "Plot_output/Eh_vs_V_"+Temp+"mk_TIME"+resoCAT+".pdf" ;
	     std::string name_plotneg = "Plot_output/Eh_kee_vs_V_"+Temp+"mk_TIME"+resoCAT+".pdf" ;
	
	
	TCanvas *c = new TCanvas("c","c",800,800);    

        
        
   c->cd();
	c->SetLogy();
	c->SetLogx();
	c->SetGridx();
	c->SetGridy();
	
	TH2D * axes = new TH2D("axes", "", 100, 0.1, 12., 100, 10., 10000000  );
	axes->GetXaxis()->SetTitle(" E_{phonon} (keV)");
	axes->GetXaxis()->SetTitleOffset(1.25);
	axes->GetYaxis()->SetTitle("NEvent.keV ^{-1}.day^{-1}.kg^{-1}");
	axes->GetYaxis()->SetTitleOffset(1.3);
	axes->GetYaxis()->SetRangeUser(20,10000000);
	
	axes->Draw();
	axes->SetStats(kFALSE); 

        ifstream List_file(Input_list.c_str(),ios::in);
        std::string Input_file;
        
	 int nline = 0 ; 
        while(std::getline(List_file, Input_file)){
        
        std::cout<<" Opening "<<Input_file<<std::endl;
	std::string runname; 
	std::string runname2; 
	std::vector <std::string> Run_Name;
	std::vector <std::string> Run_time_hists;
	ifstream Runlistname(List_ofrunandtension.c_str(),ios::in);
        while(!Runlistname.eof() ) {
        
       		Runlistname>>runname>>runname2;
       		TString RUN = runname ;
       		if(RUN.Contains("#", TString::kIgnoreCase)) continue ;
       		if(!RUN.Contains("ellapsed_time", TString::kIgnoreCase)){
       			Run_Name.push_back(runname);
       		}else{
       			Run_time_hists.push_back(runname);
       		}
        	//std::cout<<" Run name " << runname<<" "<<runname2<<std::endl;
	       	 
	}
	
	TH1D * Hist_vec           [Run_Name.size()];
	TH1D * Hist_vec_normtime  [Run_time_hists.size()];
	TH1D * Hist_vec_kee       [Run_Name.size()];
	
	Double_t Tension[Run_Name.size()];
	Double_t Norm_Tension[Run_time_hists.size()];
	std::string Tension_name;
	Double_t Sign = 0 ; 
	
	std::string tension_ [11];
	
	tension_ [0] = "0.000000";
	tension_ [1] = "8.000000";
	tension_ [2] = "15.000000";
	tension_ [3] = "30.000000";
	tension_ [4] = "40.000000";
	tension_ [5] = "66.000000";
	tension_ [6] = "70.000000";
	tension_ [7] = "78.000000";
	tension_ [8] = "60.000000";
	tension_ [9] = "51.000000";
	tension_ [10] = "57.000000";
	TLegend *leg = new TLegend(0.6, 0.6, .9, .9);
	TLegend *legneg = new TLegend(0.7, 0.7, .9, .9);
	TFile * Inputs = TFile::Open(Input_file.c_str());
	Double_t temp_tension = -999;
	TH1D * Hist_per_tension [16];
	TH1D * temp_Hist [16];
	Double_t who_first [16];


	for(int i = 0; i < Run_time_hists.size()  ; i++){
			
		Hist_vec_normtime     [i] = (TH1D*) Inputs->Get((Run_time_hists.at(i)).c_str());
		TString namehist = Hist_vec_normtime [i]->GetName();
		if(namehist.Contains("pos", TString::kIgnoreCase) == 1) Sign = +1;
		
		for(int j = 0 ; j < 11 ; j++){
		
			if(namehist.Contains((tension_ [j]).c_str(), TString::kIgnoreCase) == 1) {
				
				Norm_Tension [i] = Hist_vec_normtime[i]->Integral();
				//std::cout<<" tension "<<  tension_ [j] << " ellapsed time "<<Norm_Tension [i] <<std::endl;
				//Tension [i] *= Sign ;
				continue;
				
			}
			
		}
		
	}
	
	
	int stop_it [Run_Name.size()][Run_Name.size()];
	int IS_use [Run_Name.size()] ;
	for(int i = 0; i < Run_Name.size()  ; i++){
			
		Hist_vec     [i] = (TH1D*) Inputs->Get((Run_Name.at(i)).c_str());
		Hist_vec_kee [i] = (TH1D*) Inputs->Get((Run_Name.at(i)+"_keVee").c_str());
		Sign = -1;
		IS_use [i]= 0;
		//if(Hist_vec [i]->GetName() == NULL) continue ; 
		TString namehist = Hist_vec [i]->GetName();
		if(namehist.Contains("pos", TString::kIgnoreCase) == 1) Sign = +1;
		
		for(int j = 0 ; j < 11 ; j++){
		
			if(namehist.Contains((tension_ [j]).c_str(), TString::kIgnoreCase) == 1) {
				
				Tension [i] = std::stod(tension_ [j]);
				//Tension [i] *= Sign ;
				continue;
				
			}
			
		}
		
		
		
		Make_Nice_Hist_time(Hist_vec [i], Tension [i], Norm_Tension [i], " E_{phonon} (keV)", 0,nline);
		leg->AddEntry(Hist_vec [i],(to_string(int (Tension [i]))+"V ;"+to_string(int(Norm_Tension [i]/(3600.*24.)))+" days").c_str(),"l");
		//std::cout<<" histo Integral "<<Hist_vec [i]->Integral()<<std::endl;
		TFile* output = new TFile(("Output_graph_"+Temp+".root").c_str(),"UPDATE");
		Hist_vec [i]->Write((to_string(nline)+"_"+to_string(int (Tension [i]))+"V").c_str());
		std::cout<< " Integral f(t) "<<Hist_vec [i]->Integral(0.03,1.,"WIDTH")<<"  "<< Hist_vec [i]-> Integral()<<std::endl;
		output->Close();
		
	}
	
	nline++;
	Inputs->Close();
	}
	
	
		
	//c->SaveAs(name_plot.c_str());
	
	


}

void Study_rate_HO(std::string list_file, std::string temp, std::string Volt){
   
      ifstream List_file(list_file.c_str(),ios::in);
      std::string Input_file;
      TString buffer ;   
	   int nline = 0 ; 
      while(std::getline(List_file, Input_file)){
         buffer = Input_file ;
         if(buffer.Contains("#", TString::kIgnoreCase))continue;
         TFile* file = TFile::Open(Input_file.c_str());
         TH2D * h2_chi2_E = (TH2D*) file->Get(("Ephonon_pos"+Volt+"_22.mk_vs_chi2").c_str());
         TH1D * h1_time_norm = (TH1D*) file->Get(("Ephonon_pos"+Volt+"_22.mk_ellapsed_time").c_str());
         h2_chi2_E->Scale(1./(((h1_time_norm->Integral()/(3600.*24.)))* 0.03));
         std::cout<<"zone A Run "<< Input_file<<" Rate  "<< h2_chi2_E->Integral(104,151,99,3999)<<std::endl;
         
         std::cout<<"zone B Run "<< Input_file<<" Rate "<< h2_chi2_E->Integral(104,151,999,9999)<<std::endl;
         
         std::cout<<"zone C Run "<< Input_file<<" Rate "<< h2_chi2_E->Integral(219,249,199,999)<<std::endl;
         
         std::cout<<"zone Analysis Run "<< Input_file<<" Rate  "<< h2_chi2_E->Integral(14,37,4,12)<<std::endl;
                  
      }
        
       
       




}



int main(int argc, char** argv) {


	TCLAP::CmdLine cmd("Plotting instance for DM edelweiss searches", ' ', "0.1");
	TCLAP::ValueArg<std::string> Input_file("i", "input-file", "input file", true,"", "string");
   TCLAP::ValueArg<std::string> inputList("", "input-list", "Text file containing input files", true, "input.list", "string");
   cmd.xorAdd(Input_file, inputList);
   TCLAP::ValueArg<std::string> Detector_name("d", "detector", "Which detector", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Output_name("o", "output-name", "Output name", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Temperature("t", "temp", "temperature", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Tension("", "tension", " list of tension to study", true, "", "string", cmd);
   TCLAP::SwitchArg IsSpectrum("", "spectrum", "spectrum?", cmd);
   TCLAP::SwitchArg Islimit("", "lim", "limits?", cmd);
   TCLAP::SwitchArg Istimedep("", "time-dep", "time studies?", cmd);
   TCLAP::SwitchArg IsrateHO("", "rateHO", "HO rate studies?", cmd);
   cmd.parse(argc, argv);
   get_eff_from_file("EfficiencyRED30.root");
   if(Istimedep.getValue())  Launch_plotting_list(inputList.getValue(),Temperature.getValue(), Tension.getValue(), Output_name.getValue(),Islimit.getValue());
   if(IsSpectrum.getValue() && Input_file.isSet()) Launch_plotting(Input_file.getValue(), Temperature.getValue(), Tension.getValue(), Output_name.getValue(),Islimit.getValue(),Detector_name.getValue()); 
   if(IsrateHO.getValue())  Study_rate_HO(inputList.getValue(),Temperature.getValue(),Tension.getValue());
std::cout<<" Ending Routine "<<std::endl;
    delete Global_eff;
	return 0;	 
}

