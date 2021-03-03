#include "Histo_reader_from_rootfile.h"
#include "tclap/CmdLine.h"

using namespace std;
//------------------------------------- Class constructor-------------------------------------------------------------------------------------------
Histo_reader_from_rootfile::Histo_reader_from_rootfile( const std::string file_list, const std::string &list_histo , const std::string detector, const std::string outputname, bool verbosity, bool runontree, double voltage_user, std::string run_arg, bool const & do_kernel, bool const & do_dru )	{
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
	Do_kernel      = do_kernel;
	Do_DRU         = do_dru ;
//------------------------------------- Definition of histo for modeling choice of binning and boundaries are important !------------------------------------------------------------------------------------
	Histo_from_tree = new TH1D("Ephonon_spectrum_from_tree","Ephonon_spectrum_from_tree",4000,0.01,300.);
	Histo_from_tree->Sumw2();	
	if(Detector == "RED30"){
	    Detector_Gemass = 0.034; 
	}else if(Detector == "FID848"){
	    Detector_Gemass = 0.87; 
	}else if(Detector == "NbSi209"){ 
	   Detector_Gemass = 0.2;
	}
//------------------------------------- getting the inputs data-------------------------------------------------------------------------------------------
	std::string type_to_extract = "TH1D" ;
	Open_and_store_from_list(list_histo, type_to_extract.c_str() );
	type_to_extract = "TFile";
	Open_inputfile(file_list);		
	Get_histo_from_File();
	std::vector<double> data_to_smooth =  Fill_data_array_vec("Ep");
	std::vector<double> data_weight    =  Fill_data_array_vec("weight");
//------------------------------------- Modeling method kernel and/or Analytics-------------------------------------------------------------------------------------------
	Optimise_smoothing(data_to_smooth);	  
	
}

//------------------------------------- method to get input histogram for comparison sake-------------------------------------------------------------------------------------------
void Histo_reader_from_rootfile::Get_histo_from_File(){
		for(unsigned int i = 0 ; i < Vec_hist_name.size() ; i ++ ){
			if(!(Input_file->Get((Vec_hist_name.at(i)).c_str()))) continue;
			TH1D * temp_Hist  = (TH1D*) Input_file->Get((Vec_hist_name.at(i)).c_str());
			TH1D * timed_Hist = (TH1D*)Input_file->Get(Form("Ephonon_pos%2.0f_ellapsed_time",Voltage_user));
			TParameter<double> * ellapsed_time = (TParameter<double>*)Input_file->Get(Form("timed_lenght_pos%2.0fV",Voltage_user));
			temp_Hist->Scale(1./ellapsed_time->GetVal());
			Vec_Histo_to_fit.push_back(temp_Hist);
			
			if(IsVerbose){
			        std::cout<<" getting  histo "<<	Vec_hist_name.at(i)<<" integral histo "<<Vec_Histo_to_fit.at(i)->Integral()<<std::endl;
			        std::cout<<"time "<<timed_Hist->Integral()<<std::endl;
			        std::cout<<" "<<std::endl;
			}
		}
}

//------------------------------------- method to get input  file and set address to desired root trees -------------------------------------------------------------------------------------------
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
//------------------------------------- method to get input  file and set address to desired root trees -------------------------------------------------------------------------------------------
void Histo_reader_from_rootfile::Open_inputfile(std::string const & name){

    TFile * temp_file = new TFile(Form("%s",name.c_str()));
	Input_file = temp_file;
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
//------------------------------------- method to fill static table for gaussian kernel computation-------------------------------------------------------------------------------------------
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
//------------------------------------- method to fill array table for gaussian kernel computation and filling histo for analytical-------------------------------------------------------------------------
std::vector<double>  Histo_reader_from_rootfile::Fill_data_array_vec(std::string name_variable){
        std::vector<double>  data (event_tree -> GetEntries());
        if(IsVerbose) std::cout<< " Storing data to array  "<<std::endl;
        double binwidth = (300. - 0.01)/4000.;
        for(unsigned int i = 0 ; i < event_tree -> GetEntries() ; i++){     
                event_tree -> GetEntry(i);
                if(fabs(voltage) != Voltage_user ) continue;
                if(name_variable == "Ep" && E_p < 300. && E_p > 0.01 ){ 
                    data[i] = E_p;
                    if(Do_DRU){
                        Histo_from_tree->Fill(E_p,1./(Time_exposition *Detector_Gemass*binwidth));
                    }else{
                        Histo_from_tree->Fill(E_p,1./(binwidth));
                    }
                }
                if(name_variable == "weight") data[i] = weight;       
        }
        std::cout <<" integral Hist from tree "<<Histo_from_tree->Integral()<<std::endl;
        return data;
}
//------------------------------------- method to fit data -------------------------------------------------------------------------------------------
void Histo_reader_from_rootfile::Optimise_smoothing( std::vector<double>  data){
      if(IsVerbose) std::cout<< " smoothing spectrum  by TKDE"<<std::endl;
      TCanvas * c = new TCanvas("c","c", 1000,1000);
      c->SetLogx();
      c->SetLogy();
      TH2D * axes =new TH2D("","",1000,1e-2,300,1000,1e-1,1e8);
      axes->SetStats(kFALSE);
      axes->Draw();
      Histo_from_tree->SetStats(kFALSE);
      double binwidth = (60. - 0.01)/4000.;
      Histo_from_tree->Draw("HISTSAME");
      TLegend * leg = new TLegend(0.1,0.1,0.5,0.5);
      leg->AddEntry(Histo_from_tree,"Data","l");
      int color = 1 ;
      std::multimap<double, double> MAP_rho_chi2;
      
        double chi2buffer = 0.;
        double iter = 0.5 ;
        int Nloop = 0.;
        do{
            if(!Do_kernel) break;
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
            if(IsVerbose){
                std::cout<< " smoothing spectrum  by smooth from TH1"<<std::endl;          
                std::cout<<" chi2 /ndof "<<calculatechi2(Histo_from_tree,hist_kernel,0.08,1.5)  <<std::endl;
            }
            leg->AddEntry(hist_kernel,Form("KDE rho=%1.2f #chi^{2} / ndof=%2.2f ",iter,calculatechi2(Histo_from_tree,hist_kernel,0.5,10.)),"l");
            MAP_rho_chi2.insert(std::make_pair(calculatechi2(Histo_from_tree,hist_kernel,0.5,10.),iter)); 
            if(IsVerbose || Nloop % 10 == 0)std::cout<<" chi  "<<calculatechi2(Histo_from_tree,hist_kernel,0.5,10.) <<" rho "<<  iter <<std::endl;
            chi2buffer= calculatechi2(Histo_from_tree,hist_kernel,0.5,10.);
            iter+=0.1;
            Nloop++;
        
        }while(chi2buffer < 1.09);
        
        if(Do_kernel){
            leg->Draw();
            c->SaveAs(Form("Smoothed_kernel_%s.pdf",Detector.c_str()));
        }
        TH1D* hist_buff_kernel;
        double opti_rho     ;
        double final_chi2   ;
        unsigned int N_data ;
        if(Do_kernel){
            double opti_rho     =  MAP_rho_chi2.lower_bound(1.)->second;
            double final_chi2   = MAP_rho_chi2.lower_bound(1.)->first;
            unsigned int N_data = event_tree -> GetEntries();
            kde_opti = new TKDE(N_data, &data[0], 0.01,15., "", opti_rho);
            TH1D* hist_buff_kernel =(TH1D*) kde_opti->GetFunction(1000)->GetHistogram();
            hist_buff_kernel->Scale(1./hist_buff_kernel->Integral("WIDTH"));
        }
        if(Do_kernel)hist_buff_kernel->Scale(Histo_from_tree->Integral("WIDTH"));

//------------------------------------- Function used for modelisation here double exponential + cte + 3 fixed gaussian for Ge K L M ---------------------------------------------
        TF1 * f_bkground = new TF1("analytical_bkg","[0]*exp([1]*x)+[2]*exp([3]*x) + [4] + [5]*exp(-0.5*pow(((x-[6])/[7]),2)) + [8]*exp(-0.5*pow((x-[9])/[10],2)) + [11]*exp(-0.5*pow((x-[12])/[13],2))",6.9e-1,300.);
       
        f_bkground->SetParameter(0,2.65558e+06);
        f_bkground->SetParameter(1,-6.74769);
        f_bkground->SetParameter(2,579.215);
        f_bkground->SetParameter(3,-0.676673);
        f_bkground->SetParameter(4,5.17614);
        f_bkground->SetParameter(5,52.2695);
        f_bkground->SetParameter(6,3.66);
        f_bkground->SetParameter(7,2.87e-01);
        f_bkground->SetParameter(8, 105.193);
        f_bkground->SetParameter(9,2.97e+1);
        f_bkground->SetParameter(10,9.44e-01);
        f_bkground->SetParameter(11, 226.983);
        f_bkground->SetParameter(12,2.38e+2);
        f_bkground->SetParameter(13,3.5);
        f_bkground->FixParameter(6,3.66);
        f_bkground->FixParameter(7,2.87e-01);
        f_bkground->FixParameter(9,2.97e+1);
        f_bkground->FixParameter(10,9.44e-01);
        f_bkground->FixParameter(12,2.38e+2);
        f_bkground->FixParameter(13,3.5);
        Histo_from_tree->Fit(f_bkground,"rWL");
        
        f_bkground->SetParameter(0,f_bkground->GetParameter(0));
        f_bkground->SetParameter(1,f_bkground->GetParameter(1));
        f_bkground->SetParameter(2,f_bkground->GetParameter(2));
        f_bkground->SetParameter(3,f_bkground->GetParameter(3));
        f_bkground->SetParameter(4,f_bkground->GetParameter(4));
        f_bkground->SetParameter(5,f_bkground->GetParameter(5));
        f_bkground->SetParameter(6,f_bkground->GetParameter(6));
        f_bkground->SetParameter(7,f_bkground->GetParameter(7));
        f_bkground->SetParameter(8,f_bkground->GetParameter(8));
        f_bkground->SetParameter(9,f_bkground->GetParameter(9));
        f_bkground->SetParameter(10,f_bkground->GetParameter(10));
        f_bkground->SetParameter(11,f_bkground->GetParameter(11));
        f_bkground->SetParameter(12,f_bkground->GetParameter(12));
        f_bkground->SetParameter(13,f_bkground->GetParameter(13));
        
        f_bkground->FixParameter(6,3.63);
        f_bkground->FixParameter(7,2.8e-01);
        f_bkground->FixParameter(9,2.98e+1);
        f_bkground->FixParameter(10,7.78e-01);
        f_bkground->FixParameter(12,2.38e+2);
        f_bkground->FixParameter(13,4.08);
        Histo_from_tree->Fit(f_bkground,"rWL");
        Int_t fit_status = Histo_from_tree->Fit(f_bkground,"rWL");
        int NFit = 0 ;
//------------------------------------- fit until fit status = converged ! 4 times max -----------------------------------------------------------------------
        do{

            f_bkground->SetParameter(0,f_bkground->GetParameter(0));
            f_bkground->SetParameter(1,f_bkground->GetParameter(1));
            f_bkground->SetParameter(2,f_bkground->GetParameter(2));
            f_bkground->SetParameter(3,f_bkground->GetParameter(3));
            f_bkground->SetParameter(4,f_bkground->GetParameter(4));
            f_bkground->SetParameter(5,f_bkground->GetParameter(5));
            f_bkground->SetParameter(6,f_bkground->GetParameter(6));
            f_bkground->SetParameter(7,f_bkground->GetParameter(7));
            f_bkground->SetParameter(8,f_bkground->GetParameter(8));
            f_bkground->SetParameter(9,f_bkground->GetParameter(9));
            f_bkground->SetParameter(10,f_bkground->GetParameter(10));
            f_bkground->SetParameter(11,f_bkground->GetParameter(11));
            f_bkground->SetParameter(12,f_bkground->GetParameter(12));
            f_bkground->SetParameter(13,f_bkground->GetParameter(13));
            f_bkground->FixParameter(6,3.63);
            f_bkground->FixParameter(7,2.8e-01);
            f_bkground->FixParameter(9,2.98e+1);
            f_bkground->FixParameter(10,7.78e-01);
            f_bkground->FixParameter(12,2.38e+2);
            f_bkground->FixParameter(13,4.08);
            Histo_from_tree->Fit(f_bkground,"rwL");
            NFit++;
        }while(fit_status == 4 && NFit < 5);
        if(Do_kernel){hist_buff_kernel->SetName("Smooth_spectrum");
        hist_buff_kernel->SetFillStyle(0);}              
        TFile * ouput_smooth_spectrum = new TFile(Form("Smooth_%s_%s.root",Detector.c_str(),RUN.c_str()),"RECREATE");
        if(Do_kernel)hist_buff_kernel->Write();
        f_bkground->SetNpx(1500000) ;
        f_bkground->Write();
        TParameter<double>* Nhit = new TParameter<double>("Ncount",Histo_from_tree->Integral(Histo_from_tree->FindBin(0.690),Histo_from_tree->FindBin(300.),"Width"));
        Nhit->Write();
        TH1D * tempt_Hist = (TH1D*)Input_file->Get(Form("Ephonon_pos%2.0f_tot_noweight",Voltage_user));
        tempt_Hist->SetLineColor(kRed);
        f_bkground->SetLineColor(kRed + 3);
        f_bkground->SetTitle("");
        f_bkground->Draw("Same");
        Histo_from_tree->Draw("HISTSAME");
        if(Do_kernel)hist_buff_kernel->Draw("HISTSAME");        
        Histo_from_tree->Write("no_weight_spectrum");
        if(true){
            TParameter<double>* exposition = new TParameter<double>("Renorm_to_DRU",Time_exposition *Detector_Gemass);
            TParameter<double>* Integral_HO_LE = new TParameter<double>("Integral_HO_LE",tempt_Hist->Integral(tempt_Hist->FindBin(1.),tempt_Hist->FindBin(2.),"Width"));
            TParameter<double>* Integral_HO_HE = new TParameter<double>("Integral_HO_HE",tempt_Hist->Integral(tempt_Hist->FindBin(2.),tempt_Hist->FindBin(5.),"Width"));
            exposition->Write();
            Integral_HO_LE->Write();
            Integral_HO_HE->Write();
        }
        leg = new TLegend(0.55,0.6,0.9,0.9);
        leg->AddEntry(Histo_from_tree,"Data","l");
        leg->AddEntry((TObject*)0,Form("#frac{#chi^{2}}{N_{dof}} = %2.2f",f_bkground->GetChisquare() / f_bkground->GetNDF()),"");
        if(Do_kernel)leg->AddEntry(hist_buff_kernel,Form("KDE rho=%1.2f #chi^{2} / ndof=%2.2f ",opti_rho,final_chi2),"l");
        leg->AddEntry(f_bkground,"Analytical Fit [0]*exp([1]*x)+[2]*exp([3]*x)+[4]","l");
        leg->AddEntry((TObject*)0,Form("[0] = %1.2e +- %1.2e ",f_bkground->GetParameter(0),f_bkground->GetParError(0)),"");
        leg->AddEntry((TObject*)0,Form("[1] = %1.2e +- %1.2e ",f_bkground->GetParameter(1),f_bkground->GetParError(1)),"");
        leg->AddEntry((TObject*)0,Form("[2] = %1.2e +- %1.2e ",f_bkground->GetParameter(2),f_bkground->GetParError(2)),"");
        leg->AddEntry((TObject*)0,Form("[3] = %1.2e +- %1.2e ",f_bkground->GetParameter(3),f_bkground->GetParError(3)),"");
        leg->AddEntry((TObject*)0,Form("[4] = %1.2e +- %1.2e ",f_bkground->GetParameter(4),f_bkground->GetParError(4)),"");
        leg->Draw();
        
        ofstream streamoutparam("parametre_vs_time.txt",ios::app);
        streamoutparam<<RUN<<" "<<f_bkground->GetParameter(0)<<" "<<f_bkground->GetParError(0)<<" "<<f_bkground->GetParameter(1)<<" "<<f_bkground->GetParError(1)<<" "<<f_bkground->GetParameter(2)<<" "<<f_bkground->GetParError(2)<<" "<<f_bkground->GetParameter(3)<<" "<<f_bkground->GetParError(3)<<" "<<f_bkground->GetParameter(4)<<" "<<f_bkground->GetParError(4)<<std::endl;
        streamoutparam.close();
        
        ofstream streamoutpeak("Peak_vs_time.txt",ios::app);
        Double_t error_int1 = 0;
        Double_t error_int2 = 0;
        Double_t error_int3 = 0;
        Double_t error_int4 = 0;
        Double_t int1 = Histo_from_tree->IntegralAndError(Histo_from_tree->FindBin(3),Histo_from_tree->FindBin(4.2),error_int1);
        Double_t int2 = Histo_from_tree->IntegralAndError(Histo_from_tree->FindBin(26),Histo_from_tree->FindBin(33),error_int2);
        Double_t int3 = Histo_from_tree->IntegralAndError(Histo_from_tree->FindBin(10),Histo_from_tree->FindBin(15),error_int3);
        Double_t int4 = Histo_from_tree->IntegralAndError(Histo_from_tree->FindBin(228),Histo_from_tree->FindBin(248),error_int4);
//------------------------------------- 1.3 keV / 10.36 keV ratio should be around 11 % -----------------------------------------------------------------------
        std::cout<<"10.37 keV over 1.3 keV line ratio "<< int2 / int4 <<std::endl;
        streamoutpeak<<RUN<<" "<<int1<<"  "<<error_int1<<" "<<int2<<"  "<<error_int2<<" "<<int3<<"  "<<error_int3<<" "<<int4<<" "<<error_int4<<std::endl;
        streamoutpeak.close();
        if(Do_DRU){
            c->SaveAs(Form("Smoothed_picked_kernel_DRU_%s_%s.png",Detector.c_str(),RUN.c_str()));
        }else{
            c->SaveAs(Form("Smoothed_picked_kernel_%s_%s.png",Detector.c_str(),RUN.c_str()));
        }
        ouput_smooth_spectrum->Close();
}
//------------------------------------- Method to display numerous histo on the same pad-----------------------------------------------------------------------
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

//------------------------------------- Method to estimate the quality of gaussian kernel modelisation-----------------------------------------------------------------------
double Histo_reader_from_rootfile::calculatechi2(TH1D *hdata,TH1* hkde,double emin,double emax)
{
     int binstart=hdata->FindBin(emin);
     int binend=hdata->FindBin(emax);
     int ntotbinshist=hdata->GetNbinsX();
     if(binend>ntotbinshist){binend=ntotbinshist;}
     if(IsVerbose){
        cout<<"emin ="<<emin<<" bin = "<<binstart<<endl;
        cout<<"emax ="<<emax<<"    bin = "<<binend<<endl;
        cout<<"Nbinsmax ="<<hdata->GetNbinsX()<<endl;
     }
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
     if(IsVerbose) cout<<chi2<<"    "<<ndf<<" "<<chi2/ndf<<endl;
     //for(int bin=)
     return chi2/ndf;
}
int main(int argc, char** argv){  
   TCLAP::CmdLine cmd("Fitting bkg in time", ' ', "0.1");
   TCLAP::ValueArg<std::string> inputfileList("", "inputfile", " input file", true, "", "string",cmd);
   TCLAP::ValueArg<std::string> inputplotList("", "inputplot-list", "Text file containing input files", true, "", "string",cmd);
   TCLAP::ValueArg<std::string> Detector_name("d", "detector", "Which detector", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Output_name  ("o", "output-name", "Output name", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> run_arg      ("r", "Run-name", "Run name", true, "", "string", cmd);
   TCLAP::ValueArg<double> Voltage           ("", "Voltage", "Voltage Voltage", true, 0, "double", cmd);
   TCLAP::SwitchArg verboseArg               ("", "v", "verbose?", cmd);
   TCLAP::SwitchArg treeArg                  ("", "tree", "tree/hist?", cmd);
   TCLAP::SwitchArg DoKernel                 ("", "kernel", "kernel/Analytical?", cmd);
   TCLAP::SwitchArg DRU                      ("", "dru", "dru/count?", cmd);
   cmd.parse(argc, argv);
   Histo_reader_from_rootfile * reader_and_fitter_instance = new Histo_reader_from_rootfile(inputfileList.getValue(),inputplotList.getValue(),Detector_name.getValue() ,Output_name.getValue(),verboseArg.getValue(),treeArg.getValue(),Voltage.getValue(),run_arg.getValue(),DoKernel.getValue(),DRU.getValue());
   return 0;
}
