
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

using namespace std;

void Make_Nice_Hist(TH1D* hist , Double_t tension, std::string xaxisname){

	
	if(tension == 0 )	hist->SetLineColor(kBlack );
	if(tension == 8 )	hist->SetLineColor(kBlue + 1);
	if(tension == 15 )	hist->SetLineColor(kMagenta + 1);
	if(tension == 30 )	hist->SetLineColor(kRed + 1);
	if(tension == 60 )	hist->SetLineColor(kYellow + 1);
	if(tension == 66 )	hist->SetLineColor(kGreen + 1);	
	if(tension == 70 )	hist->SetLineColor(kGreen - 7);
	if(tension == 78 )	hist->SetLineColor(kYellow + 3);

	if(tension == -0 )	hist->SetLineColor(kBlack );
	if(tension == -8 )	hist->SetLineColor(kBlue + 1);
	if(tension == -15 )	hist->SetLineColor(kMagenta + 1);
	if(tension == -30 )	hist->SetLineColor(kRed + 1);
	if(tension == -60 )	hist->SetLineColor(kYellow + 1);
	if(tension == -66 )	hist->SetLineColor(kGreen + 1);	
	if(tension == -70 )	hist->SetLineColor(kGreen - 7);
	if(tension == -78 )	hist->SetLineColor(kYellow + 3);




	hist->GetXaxis()->SetTitle(xaxisname.c_str());
	hist->GetYaxis()->SetTitle("NEvent.keV ^{-1}");
	hist->GetYaxis()->SetTitleOffset(1.5);
	hist->SetTitle("");

}

void Per_tensionHist(TH1D * hist , TH1D * hist_2 , Double_t Tension,  Double_t Tension_2){

		
	if( Tension ==  Tension_2) hist -> Add(hist_2) ;

}

void Launch_plotting(std::string  Input_file, Double_t Temp, std::string List_ofrunandtension , std::string resoCAT){


	std::string runname; 
	std::string runname2; 
	std::vector <std::string> Run_Name;
	ifstream Runlistname(List_ofrunandtension.c_str(),ios::in);
        while(!Runlistname.eof() ) {
        
       		Runlistname>>runname>>runname2;
       		Run_Name.push_back(runname);
        	std::cout<<" Run name " << runname<<" "<<runname2<<std::endl;
	       	 
	}
	
	TH1D * Hist_vec [Run_Name.size()];
	TH1D * Hist_vec_kee [Run_Name.size()];
	
	Double_t Tension[Run_Name.size()];
	std::string Tension_name;
	Double_t Sign = 0 ; 
	
	std::string tension_ [9];
	
	tension_ [0] = "0.000000";
	tension_ [1] = "8.000000";
	tension_ [2] = "15.000000";
	tension_ [3] = "30.000000";
	tension_ [4] = "40.000000";
	tension_ [5] = "66.000000";
	tension_ [6] = "70.000000";
	tension_ [7] = "78.000000";
	tension_ [8] = "60.000000";
	TLegend *leg = new TLegend(0.7, 0.7, .9, .9);
	TLegend *legneg = new TLegend(0.7, 0.7, .9, .9);
	TFile * Inputs = TFile::Open(Input_file.c_str());
	Double_t temp_tension = -999;
	TH1D * Hist_per_tension [16];
	TH1D * temp_Hist [16];
	Double_t who_first [16];

	
	
	
	
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
		
		for(int j = 0 ; j < 9 ; j++){
		
			if(namehist.Contains((tension_ [j]).c_str(), TString::kIgnoreCase) == 1) {
				
				Tension [i] = std::stod(tension_ [j]);
				Tension [i] *= Sign ;
				continue;
				
			}
			
		}
		
		
		
		Make_Nice_Hist(Hist_vec [i], Tension [i], " E_{phonon} (keV)");
		Make_Nice_Hist(Hist_vec_kee [i], Tension [i], " E_{heat} (keVee)");
		
	}
	
	for(int i = 0 ; i < 16; i++){
	
		if(i < 16) temp_Hist [i] = new TH1D ( ("h"+to_string(i)).c_str(),("h"+to_string(i)).c_str(), 10 , 0, 1000);
	
	}
	
	
	
	

	//Drawing Part

	temp_Hist [0]->SetLineColor(kBlack );
	temp_Hist [1]->SetLineColor(kBlue + 1);
	temp_Hist [2]->SetLineColor(kMagenta + 1);
	temp_Hist [3]->SetLineColor(kRed + 1);
	temp_Hist [4]->SetLineColor(kYellow + 1);
	temp_Hist [5]->SetLineColor(kGreen + 1);
	temp_Hist [6]->SetLineColor(kGreen - 7);
	temp_Hist [7]->SetLineColor(kYellow + 3);
	temp_Hist [8]->SetLineColor(kBlue + 1);
	temp_Hist [9]->SetLineColor(kMagenta + 1);
	temp_Hist [10]->SetLineColor(kRed + 1);
	temp_Hist [11]->SetLineColor(kYellow + 1);
	temp_Hist [12]->SetLineColor(kGreen + 1);
	temp_Hist [13]->SetLineColor(kGreen - 7);
	temp_Hist [14]->SetLineColor(kYellow + 3);
	temp_Hist [15]->SetLineColor(kYellow + 5);
	
	
	leg->AddEntry(temp_Hist [0]," 0V ","l");
	leg->AddEntry(temp_Hist [1]," 8V ","l");
	leg->AddEntry(temp_Hist [2]," 15V ","l");
	//leg->AddEntry(temp_Hist [3]," 30V ","l");
	//leg->AddEntry(temp_Hist [4]," 40V ","l");
	leg->AddEntry(temp_Hist [5]," 66V ","l");
	leg->AddEntry(temp_Hist [15]," 60V ","l");
	leg->AddEntry(temp_Hist [6]," 70V ","l");
	leg->AddEntry(temp_Hist [7]," 78V ","l");


	legneg->AddEntry(temp_Hist [8]," -8V ","l");
	legneg->AddEntry(temp_Hist [9]," -15V ","l");
	//legneg->AddEntry(temp_Hist [10]," -30V ","l");
	//legneg->AddEntry(temp_Hist [11]," -40V ","l");
	legneg->AddEntry(temp_Hist [12]," -66V ","l");
	legneg->AddEntry(temp_Hist [13]," -70V ","l");
	legneg->AddEntry(temp_Hist [14]," -78V ","l");
	legneg->AddEntry(temp_Hist [15]," -60V ","l");
	
	
	TCanvas *c = new TCanvas("c","c",800,800);    
	c->SetLogy();
	c->SetLogx();
	for(int i = 0; i < Run_Name.size()  ; i++){

		if ( Tension [i] < 0) continue;
		if(i == 0 ){
			Hist_vec [i]->GetXaxis()->SetRangeUser(0.098, 12. );
			Hist_vec [i]->SetStats(kFALSE);
			Hist_vec [i]->DrawNormalized("HIST E ");
			
			
		}
		//std::cout<<"Integral "<<Hist_vec [i]->Integral()<<std::endl;
		Hist_vec [i]->GetXaxis()->SetRangeUser(0.098, 12. );
		Hist_vec [i]->SetStats(kFALSE);
		Hist_vec [i]->DrawNormalized("HIST E SAME");
	}
	leg->Draw();
	std::string name_plot = "Eh_vs_posV_"+to_string(Temp)+"mk_normalized"+resoCAT+".pdf" ;
	c->SaveAs(name_plot.c_str());

	
	for(int i = 0; i < Run_Name.size()  ; i++){

		if ( Tension [i] > 0) continue;

		if(i == 0 ){
			Hist_vec [i]->GetXaxis()->SetRangeUser(0.098, 12. );
			Hist_vec [i]->SetStats(kFALSE);
			Hist_vec [i]->DrawNormalized("HIST E ");
			
			
			
		}	
		Hist_vec [i]->GetXaxis()->SetRangeUser(0.098, 12. );
		Hist_vec [i]->SetStats(kFALSE);
		Hist_vec [i]->DrawNormalized("HIST E SAME");
	}
	legneg->Draw();
	std::string name_plotneg = "Eh_vs_negV_"+to_string(Temp)+"mk_normalized"+resoCAT+".pdf" ;
	c->SaveAs(name_plotneg.c_str());
	
	
	
	TCanvas *cc = new TCanvas("cc","cc",800,800);    
	cc->SetLogy();
	cc->SetLogx();
	for(int i = 0; i < Run_Name.size()  ; i++){

		if ( Tension [i] >= 0) continue;
		if(i == 0 ){
			Hist_vec_kee [i]->DrawNormalized("HIST E ");
			Hist_vec_kee [i]->SetStats(kFALSE);
			Hist_vec_kee [i]->GetXaxis()->SetRangeUser(0.0001, 0.1 );
		}
		Hist_vec_kee [i]->SetStats(kFALSE);
		Hist_vec_kee [i]->DrawNormalized("HIST E SAME");
	}
	legneg->Draw();
	name_plotneg = "Eh_kee_vs_negV_"+to_string(Temp)+"mk_normalized"+resoCAT+".pdf" ;
	cc->SaveAs(name_plotneg.c_str());
	
	for(int i = 0; i < Run_Name.size()  ; i++){

		if ( Tension [i] < 0) continue;
		if(i == 0 ){
		
			Hist_vec_kee [i]->DrawNormalized("HIST E ");
			Hist_vec_kee [i]->SetStats(kFALSE);
			Hist_vec_kee [i]->GetXaxis()->SetRangeUser(0.0001, 0.1 );
		}
		Hist_vec_kee [i]->SetStats(kFALSE);
		Hist_vec_kee [i]->DrawNormalized("HIST E SAME");
	}
	leg->Draw();
	name_plotneg = "Eh_kee_vs_posV_"+to_string(Temp)+"mk_normalized"+resoCAT+".pdf" ;
	cc->SaveAs(name_plotneg.c_str());
	
	
	
	


}

int main(int argc, char** argv) {


	

	Launch_plotting(argv[1], std::stod(argv[2]), argv[3], argv[4]);
	std::cout<<" Ending Routine "<<std::endl;
	return 0;	 
}

