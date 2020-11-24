
#include "Histo_reader_from_rootfile.h"
#include "tclap/CmdLine.h"

using namespace std;

Histo_reader_from_rootfile::Histo_reader_from_rootfile( const std::string file_list, const std::string &list_histo , const std::string detector, const std::string outputname, bool verbosity, bool runontree, double voltage_user, std::string run_arg )	{
        Vec_hist_name.clear();
	Vec_paramam_vs_time_lowE.clear() ;   
	Vec_paramam_vs_time_HightE.clear() ;                      
	Vec_Histo_to_fit.clear();
	Vec_parameters_lowE.clear(); 
	Vec_parameters_hightE.clear(); 
	Detector       = detector;
	Output_name    = outputname;
	IsVerbose      = verbosity;
	Voltage_user   = voltage_user;
	RUN            = run_arg;
	Histo_from_tree = new TH1D("Ephonon_spectrum_from_tree","Ephonon_spectrum_from_tree",1000,0.01,15.);	
	if(Detector == "RED30"){
	    Detector_Gemass = 0.034; 
	}else if(Detector == "FID848"){
	    Detector_Gemass = 0.87; 
	}else if(Detector == "NbSi209"){ 
	   Detector_Gemass = 0.2;
	}
	std::string type_to_extract = "TH1D" ;
	Open_and_store_from_list(list_histo, type_to_extract.c_str() );
	type_to_extract = "TFile";
	Open_and_store_from_list(file_list , type_to_extract.c_str());		
	Get_histo_from_File();
	std::vector<double> data_to_smooth =  Fill_data_array_vec("Ep");
	std::vector<double> data_weight    =  Fill_data_array_vec("weight");
	//smooth_spectrums(data_to_smooth,data_weight);
	Optimise_smoothing(data_to_smooth);	  
	for(int i = 0 ; i < 4 ; i++){	
	        TGraph * temp_graph = new TGraph();
	        Vec_paramam_vs_time_lowE.push_back(temp_graph);
	        temp_graph = new TGraph();
	        Vec_paramam_vs_time_HightE.push_back(temp_graph);
	}
}

void Histo_reader_from_rootfile::Get_histo_from_File(){
		for(unsigned int i = 0 ; i < Vec_hist_name.size() ; i ++ ){
			if(!(Input_file->Get((Vec_hist_name.at(i)).c_str()))) continue;
			TH1D * temp_Hist  = (TH1D*) Input_file->Get((Vec_hist_name.at(i)).c_str());
			TH1D * timed_Hist = (TH1D*)Input_file->Get(Form("Ephonon_pos%2.0f_ellapsed_time",Voltage_user));
			TParameter<double> * ellapsed_time = (TParameter<double>*)Input_file->Get(Form("timed_lenght_pos%2.0fV",Voltage_user));
			//Time_exposition = ellapsed_time->GetVal() /(3600.*24.);
			temp_Hist->Scale(1./ellapsed_time->GetVal());
			Vec_Histo_to_fit.push_back(temp_Hist);
			//Vec_Data_Hist.push_back(Fill_data_array("E_p"));
			//Vec_Data_weight.push_back(Fill_data_array("weight"));
			
			if(IsVerbose){
			        std::cout<<" getting  histo "<<	Vec_hist_name.at(i)<<" integral histo "<<Vec_Histo_to_fit.at(i)->Integral()<<std::endl;
			        std::cout<<"time "<<timed_Hist->Integral()<<std::endl;
			        std::cout<<" "<<std::endl;
			}
		}
}

void Histo_reader_from_rootfile::Fit_histo(){		
	Function_to_fit_lowE   = new TF1("HO_fit_lowE","[0]*exp(x*[1])+[2]*exp(x*[3])",0.3,0.9);;
	Function_to_fit_hightE = new TF1("HO_fit_HightE","[0]*exp(x*[1])+[2]*exp(x*[3])",1.1,2.5);;	
	for(unsigned int i = 0 ; i <Vec_Histo_to_fit.size(); i++){		
	        Vec_Histo_to_fit.at(i)->Fit(Function_to_fit_lowE,"R");
	        Vec_Histo_to_fit.at(i)->Fit(Function_to_fit_hightE,"R+");
	        std::vector<Double_t> temp_vect_1 ;
	        std::vector<Double_t> temp_vect_2 ;
	        for(unsigned int j = 0 ; j < 4 ; j++){
	                Vec_paramam_vs_time_lowE.at(j)->SetPoint(i,i,Function_to_fit_lowE->GetParameter(j));
	                Vec_paramam_vs_time_HightE.at(j)->SetPoint(i,i,Function_to_fit_hightE->GetParameter(j));
	                temp_vect_1.push_back(Function_to_fit_lowE->GetParameter(j));
	                temp_vect_2.push_back(Function_to_fit_hightE->GetParameter(j));
	                if(IsVerbose){ 
	                        std::cout<<"Low   E HO Fit parameter ["<<j<<"] = "<< Function_to_fit_lowE->GetParameter(j) <<std::endl;
	                        std::cout<<"Hight E HO Fit parameter ["<<j<<"] = "<< Function_to_fit_hightE->GetParameter(j) <<std::endl;
	                        std::cout<<" "<<std::endl;
	                }
	        }
	        Vec_parameters_lowE.push_back(temp_vect_1);
	        Vec_parameters_hightE.push_back(temp_vect_2);
	}			
}
void Histo_reader_from_rootfile::Open_and_store_from_list(const std::string &list, const std::string & type_obj){	
	std::string line_list_IN ="";
	std::cout << " opening list : "<<list<<std::endl;
	ifstream Listfile(list.c_str());
	std::vector<std::string> buffer_line ;
	if(!Listfile.fail()){
		while(  std::getline(Listfile, line_list_IN) )
		{
			TString temp_str = line_list_IN;
		       	if(temp_str.Contains("#", TString::kIgnoreCase)) continue ;
		       	buffer_line.push_back(line_list_IN);
		       	if(type_obj == "TFile"){
		       		TFile * temp_file = new TFile(line_list_IN.c_str());
		       		Input_file = temp_file;
		        }		        
		}
	}	
	if(type_obj == "TH1D")  { 
		Vec_hist_name = buffer_line	;
	}
	if(type_obj == "TFile") { 	
		std::cout<<Input_file->GetName()<<std::endl;
		event_tree = (TTree*)Input_file->Get("selected_events");
		event_tree->SetBranchAddress("Eh",&E_h);
		event_tree->SetBranchAddress("Ep",&E_p);
		event_tree->SetBranchAddress("voltage",&voltage);
		event_tree->SetBranchAddress("weight",&weight);
		TParameter<double> * ellapsed_time = (TParameter<double>*)Input_file->Get(Form("timed_lenght_pos%2.0fV",Voltage_user));
		Time_exposition = ellapsed_time->GetVal() ;
		TH1D * timed_Hist = (TH1D*)Input_file->Get(Form("Ephonon_pos%2.0f_ellapsed_time",Voltage_user));
		std::cout<<" time from  hist "<<timed_Hist->Integral() / (24.*3600.)<<" time from parameter "<< Time_exposition <<std::endl;
		
	}
	
}



void Histo_reader_from_rootfile::Store_to_output(){
        TFile * Output_file = new TFile((Output_name+"_"+Detector+".root").c_str(),"RECREATE");
        for(int i = 0 ; i < 4 ; i++){	

	        Vec_paramam_vs_time_lowE.at(i)->Write(("parameter_"+to_string(i)+"_vs_time_lowE").c_str());
	        Vec_paramam_vs_time_HightE.at(i)->Write(("parameter_"+to_string(i)+"_vs_time_HighE").c_str());
	}
        for(int i = 0 ; i < Vec_Histo_to_fit.size() ; i++){	
	        Vec_Histo_to_fit.at(i)->Write(("Histo_mois_"+to_string(i)).c_str());	        
	}	
        Output_file->Close();
}

 double * Histo_reader_from_rootfile::Fill_data_array(std::string name_variable){
        double  data [event_tree -> GetEntries()];
        if(IsVerbose) std::cout<< " Storing data to array  "<<std::endl;
        for(unsigned int i = 0 ; i < event_tree -> GetEntries() ; i++){     
                event_tree -> GetEntry(i);
                if(voltage != -66 ) continue;
                if(name_variable == "Ep") data[i] = E_p;
                if(name_variable == "weight") data[i] = weight;
                if(IsVerbose){ 
	                std::cout<<"Pushing back data bin ["<<i<<"] = "<< data[i]<<std::endl;
                   std::cout<<" "<<std::endl;
	        }        
        }
        return data;
}


std::vector<double>  Histo_reader_from_rootfile::Fill_data_array_vec(std::string name_variable){
        std::vector<double>  data (event_tree -> GetEntries());
        if(IsVerbose) std::cout<< " Storing data to array  "<<std::endl;
        for(unsigned int i = 0 ; i < event_tree -> GetEntries() ; i++){     
                event_tree -> GetEntry(i);
                if(fabs(voltage) != Voltage_user ) continue;
                if(name_variable == "Ep" && E_p < 15. && E_p > 0.01 ){ 
                    data[i] = E_p;
                    Histo_from_tree->Fill(E_p);                  
                }
                if(name_variable == "weight") data[i] = weight;       
        }
        double binwidth = (15. - 0.01)/1000.;
        //if(name_variable == "Ep")Histo_from_tree->Scale(1./(Time_exposition *Detector_Gemass*binwidth));
        std::cout <<" integral Hist from tree "<<Histo_from_tree->Integral()<<std::endl;
        return data;
}
void Histo_reader_from_rootfile::smooth_spectrums( double* data,  double* weight){
        if(IsVerbose) std::cout<< " smoothing spectrum  by TKDE"<<std::endl;
        double rho = 5.0; //default value
        TKDE * kde = new TKDE(event_tree -> GetEntries(), &data[0], 0.04,500., "", rho);
        TCanvas * c = new TCanvas("c","c", 1000,1000);
        TF1 * hk = kde->GetFunction(100000);
        TH1D* hist_kernel =(TH1D*) kde->GetFunction(100000)->GetHistogram();
        hist_kernel->Scale(1./hist_kernel->Integral("WIDTH"));
        hist_kernel->Scale(Vec_Histo_to_fit.at(0)->Integral());
        //hk->Draw();      
         Vec_Histo_to_fit.at(0)->Draw("");
       //hist_kernel->Draw("SAME");
        c->SaveAs("temp_kernel.pdf");
        if(IsVerbose) std::cout<< " smoothing spectrum  by smooth from TH1"<<std::endl;   
        for(unsigned int i = 0 ; i < Vec_Histo_to_fit.size() ; i++){  
              //  Vec_Histo_to_fit.at(i)->Smooth(); 
        }             
}


void Histo_reader_from_rootfile::smooth_spectrums( std::vector<double>  data,  std::vector<double>  weight){
      if(IsVerbose) std::cout<< " smoothing spectrum  by TKDE"<<std::endl;
        double rho = 10.0; //default value
        unsigned int N_data = event_tree -> GetEntries();
        TKDE * kde = new TKDE(N_data, &data[0], 0.01,500., "", rho);
        TCanvas * c = new TCanvas("c","c", 1000,1000);
        c->SetLogx();
        c->SetLogy();
       // TF1 * hk = kde->GetFunction(100000);
        TH1D* hist_kernel =(TH1D*) kde->GetFunction(300000)->GetHistogram();
        std::cout<<"Normalisation "<<hist_kernel->  Integral("WIDTH")<< "  "<<Vec_Histo_to_fit.at(0)->Integral()<<std::endl;
        hist_kernel->Scale(1./hist_kernel->Integral("WIDTH"));
        hist_kernel->Scale(Vec_Histo_to_fit.at(0)->Integral("WIDTH"));
        std::cout<<"Normalisation "<<hist_kernel->  Integral("WIDTH")<< "  "<<Vec_Histo_to_fit.at(0)->Integral("WIDTH")<<std::endl;
        hist_kernel->SetFillStyle(0);
        hist_kernel->SetLineColor(kRed);
        //hk->Draw();      
        Vec_Histo_to_fit.at(0)->GetXaxis()->SetRangeUser(0.08,15.);
       // Vec_Histo_to_fit.at(0)->GetYaxis()->SetRangeUser(0.0001,1e7);
        Vec_Histo_to_fit.at(0)->Draw("HIST");
        hist_kernel->Draw("HISTSAME");
        c->SaveAs("temp_kernel.pdf");
        std::multimap<double, Int_t> Ee_bin_bkg_phonon;
        for(int it_dat = 0 ;it_dat < hist_kernel-> GetNbinsX() ; it_dat++ ){    
                Ee_bin_bkg_phonon.insert(std::make_pair(hist_kernel->GetBinCenter(it_dat),it_dat));
         }
        double chi2_estimator= 0 ;
        for(unsigned int i = 1 ; i < /*Vec_Histo_to_fit.at(0)->GetNbinsX()*/ 170 ; i++){
                Int_t bin_up_dat = Ee_bin_bkg_phonon.lower_bound(Vec_Histo_to_fit.at(0)->GetBinCenter(i))->second; 
                double data          = (Vec_Histo_to_fit.at(0)->GetBinContent(i)/Vec_Histo_to_fit.at(0)->GetBinWidth(i));
                double kernel        = hist_kernel->GetBinContent(bin_up_dat)/hist_kernel->GetBinWidth(bin_up_dat);
                double sigma_kernel  = sqrt(kernel);
                std::cout<<"  data  "<< data  << " kernel "<<kernel<< " error "<<  sigma_kernel<<std::endl;
                chi2_estimator += pow((data - kernel)/sigma_kernel,2);
        }
        delete kde;
        delete hist_kernel;
        Vec_Histo_to_fit.at(0)->GetXaxis()->SetRangeUser(0.0,500.); 
}


double Histo_reader_from_rootfile::Optimise_smoothing( std::vector<double>  data){
      if(IsVerbose) std::cout<< " smoothing spectrum  by TKDE"<<std::endl;
      TCanvas * c = new TCanvas("c","c", 1000,1000);
      c->SetLogx();
      c->SetLogy();
      Histo_from_tree->SetStats(kFALSE);
      Histo_from_tree->GetYaxis()->SetRangeUser(1e-3,1e6);
      double binwidth = (15. - 0.01)/1000.;
      //Histo_from_tree->Scale(1./(Time_exposition *Detector_Gemass*binwidth));
      Histo_from_tree->Draw("HIST");
      TLegend * leg = new TLegend(0.1,0.1,0.5,0.5);
      leg->AddEntry(Histo_from_tree,"Data","l");
      int color = 1 ;
      std::multimap<double, double> MAP_rho_chi2;
      
      /*for(double iter = 0.1; iter < 3. ;  iter+= 0.3 ){  
            unsigned int N_data = event_tree -> GetEntries();
            TKDE * kde = new TKDE(N_data, &data[0], 0.01,15., "", iter);            
            TH1D* hist_kernel =(TH1D*) kde->GetFunction(100000)->GetHistogram();
            hist_kernel->Scale(1./hist_kernel->Integral("WIDTH"));
            hist_kernel->Scale(Histo_from_tree->Integral("WIDTH"));
            hist_kernel->SetFillStyle(0);           
            hist_kernel->SetLineColor(color); 
            color++;                       
            hist_kernel->Draw("HISTSAME");
            string detector = Detector ;             
            if(IsVerbose) std::cout<< " smoothing spectrum  by smooth from TH1"<<std::endl;          
            // std::cout<<" chi2 /ndof "<<calculatechi2(Histo_from_tree,hist_kernel,0.08,1.5)  <<std::endl;
            leg->AddEntry(hist_kernel,Form("KDE rho=%1.2f #chi^{2} / ndof=%2.2f ",iter,calculatechi2(Histo_from_tree,hist_kernel,0.1,2.)),"l");
            MAP_rho_chi2.insert(std::make_pair(calculatechi2(Histo_from_tree,hist_kernel,0.1,2.),iter)); 
            std::cout<<" chi  "<<calculatechi2(Histo_from_tree,hist_kernel,0.1,2.) <<" rho "<<  iter <<std::endl;   
        }*/
        double chi2buffer = 0.;
        double iter = 0.5 ;
        do{
            unsigned int N_data = event_tree -> GetEntries();
            TKDE * kde = new TKDE(N_data, &data[0], 0.01,15., "", iter);            
            TH1D* hist_kernel =(TH1D*) kde->GetFunction(100000)->GetHistogram();
            hist_kernel->Scale(1./hist_kernel->Integral("WIDTH"));
            hist_kernel->Scale(Histo_from_tree->Integral("WIDTH"));
            hist_kernel->SetFillStyle(0);           
            hist_kernel->SetLineColor(color); 
            color++;                       
            hist_kernel->Draw("HISTSAME");
            string detector = Detector ;             
            if(IsVerbose) std::cout<< " smoothing spectrum  by smooth from TH1"<<std::endl;          
            // std::cout<<" chi2 /ndof "<<calculatechi2(Histo_from_tree,hist_kernel,0.08,1.5)  <<std::endl;
            leg->AddEntry(hist_kernel,Form("KDE rho=%1.2f #chi^{2} / ndof=%2.2f ",iter,calculatechi2(Histo_from_tree,hist_kernel,0.1,5.)),"l");
            MAP_rho_chi2.insert(std::make_pair(calculatechi2(Histo_from_tree,hist_kernel,0.1,5.),iter)); 
            std::cout<<" chi  "<<calculatechi2(Histo_from_tree,hist_kernel,0.1,5.) <<" rho "<<  iter <<std::endl;
            chi2buffer= calculatechi2(Histo_from_tree,hist_kernel,0.1,5.);
            iter+=0.1;
        
        }while(chi2buffer < 1.1);
        
        leg->Draw();
        c->SaveAs(Form("Smoothed_kernel_%s.pdf",Detector.c_str()));// 
        
        double opti_rho     =  MAP_rho_chi2.lower_bound(1.)->second;
        double final_chi2   = MAP_rho_chi2.lower_bound(1.)->first;
        unsigned int N_data = event_tree -> GetEntries();
        kde_opti = new TKDE(N_data, &data[0], 0.01,15., "", opti_rho);
        TH1D* hist_buff_kernel =(TH1D*) kde_opti->GetFunction(1000)->GetHistogram();
        hist_buff_kernel->Scale(1./hist_buff_kernel->Integral("WIDTH"));
        Histo_from_tree->Scale(1./(/*Time_exposition *Detector_Gemass*/binwidth));
        hist_buff_kernel->Scale(Histo_from_tree->Integral("WIDTH"));
        Histo_from_tree->GetXaxis()->SetRangeUser(0.01,15.);
        hist_buff_kernel->SetName("Smooth_spectrum");
        hist_buff_kernel->SetFillStyle(0);
        TFile * Input_weighted_spectrum = new TFile("../PLOTTING/FEVRIER_TEST.root");
        TH1D  * tempt_Hist_wgt = (TH1D*)Input_weighted_spectrum->Get("Fevrier_phonon_NbSi209_51");      
        TFile * ouput_smooth_spectrum = new TFile(Form("Smooth_%s_%s.root",Detector.c_str(),RUN.c_str()),"RECREATE");
        hist_buff_kernel->Write();     
        TH1D * tempt_Hist = (TH1D*)Input_file->Get(Form("Ephonon_pos%2.0f_tot_noweight",Voltage_user));
        //tempt_Hist->Scale(1./(Detector_Gemass));
        tempt_Hist->SetLineColor(kRed);
        //tempt_Hist_wgt->Draw("HIST");
        //tempt_Hist->Draw("HISTSAME");
        Histo_from_tree->Draw("HIST");
        hist_buff_kernel->Draw("HISTSAME");
        tempt_Hist->Write("no_weight_spectrum");
        //tempt_Hist_wgt->Write("Weighted_spectrum");
        TParameter<double>* exposition = new TParameter<double>("Renorm_to_DRU",Time_exposition *Detector_Gemass);
        exposition->Write();
        leg = new TLegend(0.1,0.1,0.5,0.5);
        leg->AddEntry(Histo_from_tree,"Data","l");
        //tempt_Hist_wgt->Draw("HISTSAME");
        leg->AddEntry(hist_buff_kernel,Form("KDE rho=%1.2f #chi^{2} / ndof=%2.2f ",opti_rho,final_chi2),"l");
        leg->Draw();
        
        c->SaveAs(Form("Smoothed_picked_kernel_%s.pdf",Detector.c_str()));// 
        ouput_smooth_spectrum->Close();
        Input_weighted_spectrum->Close();
         
}
int Histo_reader_from_rootfile::color_rainbow(int i,int nbval){
	int color;
    //100,798,845,64,52
    // 100 red
    // 64 bleu clair
    // 798 orane
    // 845 vert
    // (2 vert
    if(nbval==1){return 1;}
    double tempcolor = 51+ (100-51.)*i/(nbval-1);
	color = (int)tempcolor;		
   if(nbval==3){
    	if(i==0){color=100;}
    	if(i==1){color=52;} 
    	if(i==2){color=64;}	
    }	    
     if(nbval==4){
      if(i==0){color=100;}
      if(i==1){color=52;} 
    	if(i==2){color=64;}	
    	if(i==3){color=798;}
    }	
	
    return color;
} 

double Histo_reader_from_rootfile::calculatechi2(TH1D *hdata,TH1* hkde,double emin,double emax)
{
     int binstart=hdata->FindBin(emin);
     int binend=hdata->FindBin(emax);
     int ntotbinshist=hdata->GetNbinsX();
     if(binend>ntotbinshist){binend=ntotbinshist;}
     cout<<"emin ="<<emin<<" bin = "<<binstart<<endl;
     cout<<"emax ="<<emax<<"    bin = "<<binend<<endl;
     cout<<"Nbinsmax ="<<hdata->GetNbinsX()<<endl;

     double chi2=0;
     int ndf=0;
     for(int bin=binstart;bin<=binend;bin++){
         double e=hdata->GetBinCenter(bin);
         double obs=(double)hdata->GetBinContent(bin);

         // for expected integrate over bin width of data
         int Nint=20;
         double newmin=hdata->GetBinLowEdge(bin);
         double newmax=newmin+hdata->GetBinWidth(bin);
         double newdelta=(newmax-newmin)/(double)Nint;
         double expected=0;
         for(int i=0;i<Nint;i++){
             double min=newmin+i*newdelta;
             double max=newmin+(i+1.)*newdelta;
             double emoy=(min+max)/2.;

expected+=hkde->GetBinContent(hkde->FindBin(emoy))/(double)Nint;
         }

         //double expected=hkde->GetBinContent(hkde->FindBin(e));
         if(expected>0.){
         chi2+=TMath::Power(obs-expected,2)/expected;
         ndf++;
         }
     }
     cout<<chi2<<"    "<<ndf<<" "<<chi2/ndf<<endl;
     //for(int bin=)

     return chi2/ndf;

}

int main(int argc, char** argv){  
   TCLAP::CmdLine cmd("Fitting bkg in time", ' ', "0.1");
   TCLAP::ValueArg<std::string> inputfileList("", "inputfile-list", "Text file containing input files", true, "", "string",cmd);
   TCLAP::ValueArg<std::string> inputplotList("", "inputplot-list", "Text file containing input files", true, "", "string",cmd);
   TCLAP::ValueArg<std::string> Detector_name("d", "detector", "Which detector", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Output_name  ("o", "output-name", "Output name", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> run_arg      ("r", "Run-name", "Run name", true, "", "string", cmd);
   TCLAP::ValueArg<double> Voltage           ("", "Voltage", "Voltage Voltage", true, 0, "double", cmd);
   TCLAP::SwitchArg verboseArg               ("", "v", "verbose?", cmd);
   TCLAP::SwitchArg treeArg                  ("", "tree", "tree/hist?", cmd);
   cmd.parse(argc, argv);
   Histo_reader_from_rootfile * reader_and_fitter_instance = new Histo_reader_from_rootfile(inputfileList.getValue(),inputplotList.getValue(),Detector_name.getValue() ,Output_name.getValue(),verboseArg.getValue(),treeArg.getValue(),Voltage.getValue(),run_arg.getValue());
  // reader_and_fitter_instance->Get_histo_from_File();
   //reader_and_fitter_instance->Fit_histo();
   //reader_and_fitter_instance->Store_to_output();
   return 0;
}
