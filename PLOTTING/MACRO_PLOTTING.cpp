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

using namespace std;


void Make_Nice_Hist(TH1D* hist , Double_t tension, Double_t norm, std::string xaxisname, int rebin){

	
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

	
	hist->Scale(1./(((norm/(3600.*24.)))* 0.03));
	hist->SetLineWidth(1);
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

void Launch_plotting(std::string  Input_file, std::string Temp, std::string List_ofrunandtension , std::string resoCAT,std::string plotlimit){


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
        	std::cout<<" Run name " << runname<<" "<<runname2<<std::endl;
	       	 
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
				std::cout<<" tension "<<  tension_ [j] << " ellapsed time "<<Norm_Tension [i] <<std::endl;
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
		
		
		
		Make_Nice_Hist(Hist_vec [i], Tension [i], Norm_Tension [i], " E_{phonon} (keV)", 0);
		Make_Nice_Hist(Hist_vec_kee [i], Tension [i], Norm_Tension [i], " E_{heat} (keVee)",0);
		leg->AddEntry(Hist_vec [i],(to_string(int (Tension [i]))+"V ;"+to_string(int(Norm_Tension [i]/(3600.*24.)))+" days").c_str(),"l");
		if(namehist.Contains("Ephonon_pos15.000000_22.000000mk", TString::kIgnoreCase) == 1) Fit_and_store_HO(Hist_vec [i]);
	}
	

	TFile *red20_file = TFile::Open("output_red20.root", "UPDATE");
	
	
	TGraphErrors * Red_20_graph = (TGraphErrors*) red20_file->Get("Red20_results");
	Red_20_graph->SetMarkerColor(kBlack);
	Red_20_graph->SetMarkerSize(1);
	Red_20_graph->SetMarkerStyle(34);
	leg->AddEntry(Red_20_graph, "RED20 surface spectrum", "p");
	
	
	TFile *DM_file = TFile::Open("DM_signal/Root_Files/DM_from_graph_V2.root", "UPDATE");
	
	
	TGraphErrors * DMlim = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_3e-38cm2_66.000000V");
	DMlim->SetMarkerColor(kBlue);
	DMlim->SetMarkerSize(1);
	DMlim->SetMarkerStyle(33);
	if(plotlimit == "yes")leg->AddEntry(DMlim, "WIMP 2 GeV 3e^{-38} cm ^{2}", "p");
	
	TGraphErrors * DMlim_2 = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_1e-38cm2_66.000000V");
	DMlim_2->SetMarkerColor(kBlue);
	DMlim_2->SetMarkerSize(0.7);
	DMlim_2->SetMarkerStyle(8);
	if(plotlimit == "yes")leg->AddEntry(DMlim_2, "WIMP 2 GeV 1e^{-38} cm ^{2}", "p");
	
	/*TGraphErrors * DMlim_3 = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_1e-37cm2_0.000000V");
	DMlim_3->SetMarkerColor(kBlue);
	DMlim_3->SetMarkerSize(0.7);
	DMlim_3->SetMarkerStyle(22);
	leg->AddEntry(DMlim_3, "WIMP 2 GeV 1e^{-37} cm ^{2}", "p");*/
	
	
	TGraphErrors * DMlim_test = (TGraphErrors*) DM_file->Get("DM_Ephonon_2GeV_4e-37cm2_0.000000V");
	DMlim_test->SetMarkerColor(kBlue);
	DMlim_test->SetMarkerSize(0.7);
	DMlim_test->SetMarkerStyle(22);
	if(plotlimit == "yes")leg->AddEntry(DMlim_test, "WIMP 2 GeV 4e^{-37} cm ^{2} (RED20)", "p");
	
	TGraphErrors * DMlim_4 = (TGraphErrors*) DM_file->Get("DM_Ephonon_1GeV_1e-37cm2_66.000000V");
	DMlim_4->SetMarkerColor(kRed);
	DMlim_4->SetMarkerSize(0.7);
	DMlim_4->SetMarkerStyle(33);
	if(plotlimit == "yes")leg->AddEntry(DMlim_4, "WIMP 1 GeV 1e^{-37} cm ^{2}", "p");
	
	TGraphErrors * DMlim_5 = (TGraphErrors*) DM_file->Get("DM_Ephonon_1GeV_1e-36cm2_66.000000V");
	DMlim_5->SetMarkerColor(kRed);
	DMlim_5->SetMarkerSize(0.7);
	DMlim_5->SetMarkerStyle(22);
	if(plotlimit == "yes")leg->AddEntry(DMlim_5, "WIMP 1 GeV 1e^{-36} cm ^{2}", "p");
	
	TCanvas *c = new TCanvas("c","c",800,800);    
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
	
	axes->Draw("");
	axes->SetStats(kFALSE);
	
	
	Red_20_graph->Draw("PSAME");
	if(plotlimit == "yes"){
	        DMlim       ->Draw("PSAME");
	        DMlim_2     ->Draw("PSAME");
	        //DMlim_3     ->Draw("PSAME");
	        DMlim_4     ->Draw("PSAME");
	        DMlim_5     ->Draw("PSAME");
	        DMlim_test  ->Draw("PSAME");
	}
	for(int i = 0; i < Run_Name.size()  ; i++){

		Hist_vec [i]->SetStats(kFALSE);
		Hist_vec [i]->Draw("HIST E SAME");
		
		
	}
	leg->Draw();
	std::string name_plot = "Plot_output/Eh_vs_V_"+Temp+"mk_normalized"+resoCAT+".pdf" ;
	c->SaveAs(name_plot.c_str());

	
	
	TCanvas *cc = new TCanvas("cc","cc",800,800);    
	cc->SetLogy();
	cc->SetLogx();
	
	cc->SetGridy();
	cc->SetGridx();
	
	TH2D* axes2 = new TH2D("axes2", "", 100, 0.01, 12., 100, 100., 100000000  );
	axes2->GetXaxis()->SetTitle(" E_{heat} (keVee)");
	axes2->GetXaxis()->SetTitleOffset(1.25);
	axes2->GetYaxis()->SetTitle("NEvent.keV ^{-1}.day^{-1}.kg^{-1}");
	axes2->GetYaxis()->SetTitleOffset(1.3);
	axes2->GetYaxis()->SetRangeUser(100,1000000000);
	
	axes2->Draw("");
	axes2->SetStats(kFALSE);
	
	
	Red_20_graph->Draw("PSAME");
	
	for(int i = 0; i < Run_Name.size()  ; i++){


		Hist_vec_kee [i]->SetStats(kFALSE);
		Hist_vec_kee [i]->Draw("HIST E SAME");
	}
	leg->Draw();
	std::string name_plotneg = "Plot_output/Eh_kee_vs_V_"+Temp+"mk_normalized"+resoCAT+".pdf" ;
	cc->SaveAs(name_plotneg.c_str());
	
	


}

void Launch_plotting_list(std::string  Input_list, std::string Temp, std::string List_ofrunandtension , std::string resoCAT,std::string plotlimit){

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

void Study_rate_HO(std::string list_file, std::string temp, double Volt){
   
      ifstream List_file(list_file.c_str(),ios::in);
      std::string Input_file;
      TString buffer ;   
	   int nline = 0 ; 
      while(std::getline(List_file, Input_file)){
         buffer = Input_file ;
         if(buffer.Contains("#", TString::kIgnoreCase))continue;
         TFile* file = TFile::Open(Input_file.c_str());
         TH2D * h2_chi2_E = (TH2D*) file->Get(("Ephonon_pos"+to_string(Volt)+"_22.000000mk_vs_chi2").c_str());
         TH1D * h1_time_norm = (TH1D*) file->Get(("Ephonon_pos"+to_string(Volt)+"_22.000000mk_ellapsed_time").c_str());
         h2_chi2_E->Scale(1./(((h1_time_norm->Integral()/(3600.*24.)))* 0.03));
         std::cout<<"zone A Run "<< Input_file<<" Rate  "<< h2_chi2_E->Integral(104,151,99,3999)<<std::endl;
         
         std::cout<<"zone B Run "<< Input_file<<" Rate "<< h2_chi2_E->Integral(104,151,999,9999)<<std::endl;
         
         std::cout<<"zone C Run "<< Input_file<<" Rate "<< h2_chi2_E->Integral(219,249,199,999)<<std::endl;
         
         std::cout<<"zone Analysis Run "<< Input_file<<" Rate  "<< h2_chi2_E->Integral(14,37,4,12)<<std::endl;
                  
      }
        
       
       




}



int main(int argc, char** argv) {


	std::cout<<"test "<<argv[6]<<std::endl;
	std::string arg = argv[6] ;
        if(arg ==  "timeDep"){
        
                std::cout<<"test "<<argv[6]<<std::endl;
                Launch_plotting_list(argv[7], argv[2], argv[3], argv[4],argv[5]);
        }else if (arg ==  "Spectrum"){
           std::cout<<"test "<<std::endl;
	        Launch_plotting(argv[1], argv[2], argv[3], argv[4],argv[5]);
	}else{
	      
	      Study_rate_HO(argv[7],argv[2],stod(argv[6]));
	
	}
	
	std::cout<<" Ending Routine "<<std::endl;
	return 0;	 
}

