#include "Plot_HEnergy_Voltage.h"

using namespace std;

Plot_HEnergy_Voltage::Plot_HEnergy_Voltage( const std::string & list_name_in, Double_t HEAT, bool IsRun, bool On_processed , const std::string & outputdir, bool local_list_){
	

	TH1::SetDefaultSumw2();	
	IS_PROCESSED = On_processed;
	OutputDir = outputdir ;
	system(("Create_outputdir.sh "+OutputDir).c_str());
	list_name = list_name_in;
	local_list = local_list_ ;
	if (IsRun){
		RunOnly(HEAT);
	}
	else{
		std::cout<<"enter number of Run in the list"<<std::endl;
		//cin>> Nlist;
		allRUN = 0 ;
		Parse_List();	

	}
	
	Clean();
}

void Plot_HEnergy_Voltage::RunOnly(Double_t HEAT ){
	TH1::SetDefaultSumw2();

	SetTemp(HEAT);
	SetRunname(list_name);
	std::cout<<"Create Ploting instance "<<std::endl;
	Init();
	

}

void Plot_HEnergy_Voltage::RunList(Double_t Heat, const std::string & list){
	TH1::SetDefaultSumw2();

	SetTemp(Heat);
	SetRunname(list);
	std::cout<<"Create Ploting instance "<<std::endl;
	Init();
	

}
void Plot_HEnergy_Voltage::Clean(){

	//std::cout<<" memory release "<<std::endl;
	delete chain_voltage ;
	delete chain_index ;		
	delete chain_HeatEnergy ;
	delete chain_chi2A;
	delete chain_event_processed ;
	delete chain_voltage_pro ;
	delete chain_event_processed_fast ;
	delete H_Eh;
	delete H_Eh_lowres;
	delete H_Ehee;
	delete H2_Eh_chi2;
	delete Time_per_voltage;
	delete Ionration_vs_Ei;
	delete Dchi2_vs_Ep_pass;
	delete Dchi2_vs_Ep_fail;		
	delete chi2_cut_vs_Ep_pass;
	delete chi2_cut_vs_Ep_fail;		
	delete reso_vs_time;		
	delete PSD_plot;
	delete PSD_plot_reso;		
	delete PSD_spectrum; 
	delete PSD_spectrum_summed;	
	delete Input_Files;
	delete Output_Files;
}


void Plot_HEnergy_Voltage::Help(){
	std::cout<<" Please use the following option in the command line : list_name  "<<std::endl;
}
void Plot_HEnergy_Voltage::Parse_List(){

	double Heat_cat[100] ;	
	std::string inputList = list_name;
	std::cout << " opening list : "<<inputList<<std::endl;
	ifstream Listfile(inputList.c_str());	
	char pNamemc[500];
	std::string list_name_temp="";
	double Heat_perRun; 
	int ilist = 0;
	int count_line_mc = 0;
	std::string line_list_IN ="";
	if(!Listfile.fail()){
		while(  std::getline(Listfile, line_list_IN) )
		{
			count_line_mc = 0;        	     
        	        unsigned int Index_space = line_list_IN.find(" ");
        	        std::string temp_str_heat = line_list_IN.substr(Index_space , 4);
        	        std::string temp_str_list = line_list_IN.substr(0, Index_space );
        	        Double_t temp_heat = std::stod(temp_str_heat);                
			Heat_cat[ilist]	= temp_heat; 					
		        std::cout<<" adding : /sps/edelweis/rootDataRun317/streams/prodg/lists/"<<temp_str_list<< ".list Temp "<< Heat_cat[ilist]<< std::endl;
		        string file = temp_str_list ;		       
		        std::string prefix_list = "/pbs/home/h/hlattaud/PLOTEH/"  ;		        
		        if(local_list == 0) prefix_list = "/sps/edelweis/rootDataRun317/streams/prodg/lists/" ;		        
		        ifstream efficiencies((prefix_list+file+".list").c_str());
		        std::string ListRun_name[100] = {""};
		        if(!efficiencies.fail()){
		       		 while(std::getline( efficiencies,ListRun_name[count_line_mc] )){		
					std::cout<<" List content "<< ListRun_name[count_line_mc]<<std::endl;
					system(("Create_list.sh "+ListRun_name[count_line_mc]).c_str());					
					std::cout<<" Run name ? " << ListRun_name[count_line_mc]<< " Heat  "<< Heat_cat[ilist] <<std::endl;						RunList( Heat_cat[ilist] , ListRun_name[count_line_mc] );					
		       		 	count_line_mc += 1; 
		       	 	}		        
			}
			ilist++;		
		}
	}	
}

void Plot_HEnergy_Voltage::Open_file( const std::string & file_name ){

		//Everest output

		chain_voltage    = new TChain("heatCalibData") ;
		chain_index      = new TChain("boloHeader");		
		chain_HeatEnergy = new TChain("Energies_Trig_Filt_Decor") ;
		chain_chi2A      = new TChain("Amplitudes_Trig_Filt_Decor");	
        	std::string inputListMC_ = file_name;
		ifstream ismc(("List/"+Run_name).c_str());
		count_line = 0;
		char pNamemc[500];
		while( ismc.getline(pNamemc, 500, '\n') )
		{
        		if (pNamemc[0] == '#') continue;
        		if (pNamemc[0] == '\n') continue;
        		count_line += 1; 
		        std::cout<<" adding : "<<pNamemc<<std::endl;
		        chain_voltage     ->Add(pNamemc);
		        chain_index       ->Add(pNamemc);
		        chain_HeatEnergy  ->Add(pNamemc);
			chain_chi2A       ->Add(pNamemc);
		}
		N_partition = count_line ;		
		cout<<N_partition << " partition in the Run " << endl;
		cout << "[+] Linking variable...                                 " << endl;
		chain_voltage      -> SetBranchAddress ("Voltage",&Voltage);
		chain_index        -> SetBranchAddress ("indexHeatCalibData",&Index_Calib);
		chain_HeatEnergy   -> SetBranchAddress ("Eh",&Eh);
		chain_chi2A        -> SetBranchAddress ("heatChi2A",&chi2_A);
		chain_chi2A        -> SetBranchAddress ("heatResoA",&Reso_cat);
		chain_HeatEnergy   -> SetBranchAddress ("Ei",&Ei);	
		Nb_voltage    = chain_voltage    -> GetEntries();
		Nb_index      = chain_index      -> GetEntries();
		Nb_HeatEnergy = chain_HeatEnergy -> GetEntries();
		Nb_chi2       = chain_chi2A      -> GetEntries();
	
		//Processed stuff (Nepal output) 
		
		chain_voltage_pro             = new TChain("RunTree_Normal") ;
		chain_event_processed         = new TChain("EventTree_trig_Normal_filt_decor");
		chain_event_processed_fast    = new TChain("EventTree_trig_Fast_filt");
	
        	inputListMC_ = file_name;
		ifstream ismc_2(("List/"+Run_name+"_processed").c_str());
		count_line = 0;
		 pNamemc[500];
		while( ismc_2.getline(pNamemc, 500, '\n') )
		{
        		if (pNamemc[0] == '#') continue;
        		if (pNamemc[0] == '\n') continue;
        		count_line += 1; 
		        std::cout<<" adding : "<<pNamemc<<std::endl;
		        chain_voltage_pro           ->Add(pNamemc);
		        chain_event_processed       ->Add(pNamemc);
		        chain_event_processed_fast  ->Add(pNamemc);
		}
		
		chain_voltage_pro   	   -> SetBranchAddress ("PSD_Filt",&PDS_noise);
		chain_voltage_pro   	   -> SetBranchAddress ("Polar_Ion",&polarion);
		chain_voltage_pro  	   -> SetBranchAddress ("cutoff_freq",&cutofffreq);
		chain_voltage_pro  	   -> SetBranchAddress ("filter_order",&filter_order);
		chain_voltage_pro  	   -> SetBranchAddress ("Chan_Gain",&nVtoADU);
		chain_voltage_pro   	   -> SetBranchAddress ("PSD_Freq",&PDS_freq);
		chain_voltage_pro  	   -> SetBranchAddress ("f_max_heat",&f_max_heat);		 
		chain_event_processed      -> SetBranchAddress ("MicroStp",&micro_step);
		chain_event_processed      -> SetBranchAddress ("Time_crate",&Time_Crate);		
		chain_event_processed      -> SetBranchAddress ("chi2_OF_h",&chi2_norm);
		chain_event_processed_fast -> SetBranchAddress ("chi2_OF_h",&chi2_fast);	
		cout << "[+] Linking variable... done                           " << endl;
}

void Plot_HEnergy_Voltage::Write_histo_tofile(float temp, int voltage, const std::string & run_name ){

	std::string reso_CAT="";
	std::string Processed = "";	
	//categorization based on ADU heat reso
	if(Reso_cat_buffer < 1.3) reso_CAT = "highres";
	if(Reso_cat_buffer < 2. && Reso_cat_buffer >= 1.3) reso_CAT = "mediumres";
	if( Reso_cat_buffer >= 2.) reso_CAT = "lowres";	
	//new categorisation 
	bool PSD_cat = false ;
	if(PSD_cat){
		if(Reso_cat_buffer < 2.) reso_CAT = "highres";
		if(Reso_cat_buffer < 5. && Reso_cat_buffer >= 2.) reso_CAT = "mediumres";
		if( Reso_cat_buffer >= 5.) reso_CAT = "lowres";
	}	
	if( IS_PROCESSED==1) Processed = "processed";	
	if(allRUN == 1) Output_Files = new TFile((OutputDir+"/Eh_allruns_"+to_string(temp)+"mk.root").c_str(),"UPDATE");
	if(allRUN == 0) Output_Files = new TFile((OutputDir+"/Eh_perrun_"+to_string(temp)+"mk_"+run_name.c_str()+"_"+to_string(voltage)+"_"+reso_CAT+"_"+Processed+".root").c_str(),"RECREATE");	
	TParameter<Double_t> * resoHEAT = new TParameter<Double_t>("Resolution_heat", Reso_cat_buffer);
	if( run_name != "" && IS_PROCESSED==0) {
		H_Ehee       ->Write();    
		H_Eh         ->Write();
		H2_Eh_chi2   ->Write();
		H_Eh_lowres  ->Write();
		resoHEAT     ->Write();
		Time_per_voltage->Write();
		Ionration_vs_Ei->Write();		
		Dchi2_vs_Ep_pass->Write();
		Dchi2_vs_Ep_fail->Write();		
		chi2_cut_vs_Ep_pass->Write();
		chi2_cut_vs_Ep_fail->Write();
		
	}else{
	
		if( run_name != "") PSD_plot->Write("PSD_spectrum");	
	}	
	delete resoHEAT;
	Output_Files->Close();
}

Double_t Plot_HEnergy_Voltage::EpBinIndex(Double_t pt,std::vector<Double_t> binning) {

	for (unsigned int i(0); i <  binning.size() ;  i++) {
		if (pt < binning.at(i) && pt > binning.at(1)) return binning.at(i) - binning.at(i-1) ;
		return binning.at(1) ;
	}
	return binning_vec.size();
}

Double_t Plot_HEnergy_Voltage::Kevee_weight(Double_t Eh){

	if(Eh < 1) return 0.01 ;
	if(Eh < 5) return 0.04 ;
	return 0.2;
}


void Plot_HEnergy_Voltage::Loop_over_Chain(){

// Calculate time spent in whole run
	Double_t temp_freq_heat_max = 0 ;
	chain_voltage_pro->GetEntry(3);
	temp_freq_heat_max = f_max_heat ;	
	Double_t Ellapsed_time = 0;	
	chain_event_processed->GetEntry(0);	
	Double_t time_1 = micro_step / temp_freq_heat_max ;	
	chain_event_processed->GetEntry(chain_event_processed->GetEntries() - 1 );	
	Double_t time_2 = (micro_step / temp_freq_heat_max) +  3600 * (N_partition - 1)  ;	
	Ellapsed_time = ( time_2 - time_1 );	
// Test time calculus alternatives 
	chain_event_processed->GetEntry(0);
	Double_t altime1 = 3600. - micro_step / temp_freq_heat_max ;
	chain_event_processed->GetEntry(chain_event_processed->GetEntries() - 1 );
	Double_t altime2 = altime1 + ((micro_step / temp_freq_heat_max)  + (3600. *( N_partition - 2))) ;
// End time spent in whole run	calculus.
	chain_index   ->GetEntry(0) ;
	chain_voltage ->GetEntry(Index_Calib);	
	//std::cout<<" Time 1 "<<   time_1 << " Time 2 "<< time_2 << "  alternative calculus "<< altime2 / 3600. << " h " <<std::endl;	
	std::cout<<" Voltage and Temp for run : "<< Run_name <<" "<<Voltage<<" V "<<" "<< heat<<" mK run lasted for "<< (Ellapsed_time) / 3600.<< " h" <<std::endl; 
	std::string voltname = "pos"+to_string(fabs( Voltage));
	//if(Voltage < 0) voltname = "neg"+to_string(fabs( Voltage));
	std::string histname         = "";
	std::string histname_ee      = "";
	std::string histname_calib   = "";
	if(allRUN == 0){ 	
		histname = "Ephonon_"+voltname+"_"+to_string(heat)+"mk" ;
		histname_ee = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_keVee" ;
		histname_calib = "Ioniratio_vs_Ei_"+voltname+"_"+to_string(heat)+"mk" ;		
	}else{	
		histname = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_"+Run_name ;
		histname_ee = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_"+Run_name+"_keVee" ;
		histname_calib = "Ioniratio_vs_Ei_"+voltname+"_"+to_string(heat)+"mk" ;	
	}
	std::cout<<" Creating  Histo "<<histname<<std::endl;	
	Int_t Nbin_sub100   = 0;
	Int_t Nbin          = 0;
	//Double_t Binning_keVee[251] = {0.};	
	Double_t Binning_chi2[10000] = {0.} ;	
	Double_t itbin = 0.;
	Double_t itbinchi2 = 0.;
	for(int bin = 0 ; bin < 10000; bin ++ ){
			itbinchi2 += 0.1;
			Binning_chi2[bin] = itbinchi2;		 
	}	
	Double_t Binning_Dchi2[100] = {0.} ;	
	itbinchi2 = -0.05;
	for(int bin = 0 ; bin < 100; bin ++ ){
			itbinchi2 += 0.001;
			Binning_Dchi2[bin] = itbinchi2;		 
	}		
	binning_vec_kevee.clear();	
	binning_vec_kevee.push_back(0.);
	Double_t iterator_bin_ee = 0 ;
	do{
		if(iterator_bin_ee == 0){		
			 iterator_bin_ee += 0.005;		
		}else{
			Double_t sigma =std::sqrt( std::pow(0.005,2)  + std::pow(0.02*iterator_bin_ee,2));
			iterator_bin_ee += sigma ;
			binning_vec_kevee.push_back(iterator_bin_ee);
		}
	}while(iterator_bin_ee <= 100);	
	Double_t Binning_keVee[binning_vec_kevee.size()];	
	for(int it = 0 ; it < binning_vec_kevee.size() ; it ++){	
		Binning_keVee [it] = binning_vec_kevee.at(it);	
	}	
	binning_vec.clear();	
	binning_vec.push_back(0.);
	Double_t iterator_bin = 0 ;
	do{
		if(iterator_bin == 0){		
			 iterator_bin += 0.04;		
		}else{
			Double_t sigma =std::sqrt( std::pow(0.02,2)  + std::pow(0.02*iterator_bin,2));
			iterator_bin += sigma ;
			binning_vec.push_back(iterator_bin);
		}
	}while(iterator_bin <= 500);	
	Double_t Binning_keV[binning_vec.size()];	
	for(int it = 0 ; it < binning_vec.size() ; it ++){	
		Binning_keV [it] = binning_vec.at(it);	
	}	
	binning_vec_low_res.clear();	
	binning_vec_low_res.push_back(0.);
	iterator_bin = 0 ;
	do{
		if(iterator_bin == 0){		
			 iterator_bin += 0.04;		
		}else{
			Double_t sigma =std::sqrt( std::pow(0.04,2)  + std::pow(0.05*iterator_bin,2));
			iterator_bin += sigma ;
			binning_vec_low_res.push_back(iterator_bin);
		}

	}while(iterator_bin < 500);	
	Double_t Binning_keVlow[binning_vec.size()];	
	for(unsigned int it = 0 ; it < binning_vec_low_res.size() ; it ++){	
		Binning_keVlow [it] = binning_vec_low_res.at(it);	
	}	
	H_Eh          = new TH1D((histname).c_str(), (histname).c_str(),750., 0.,300.);
	H_Eh          ->SetBins((int)binning_vec.size()-1, Binning_keV);	
	H_Eh_lowres   = new TH1D((histname+"_lowres").c_str(), (histname+"_lowres").c_str(),750., 0.,300.);
	H_Eh_lowres   ->SetBins((int)binning_vec_low_res.size()-1, Binning_keVlow);	
	H_Ehee        = new TH1D((histname_ee).c_str(), (histname_ee).c_str(),375, 0.,15.);	
	H_Ehee        ->SetBins((int) binning_vec_kevee.size()-1, Binning_keVee) ; 	
	H2_Eh_chi2    = new TH2D((histname+"_vs_chi2").c_str(), (histname+"_vs_chi2").c_str(),750, 0.,300.,100,0.,1500.); 
	H2_Eh_chi2    -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);	
	Double_t ratio_binninh [100] = {0};	
	iterator_bin = 0 ;
	for( int i = 0 ; i < 100 ; i++){	
	 	ratio_binninh [i] = iterator_bin ;	 	
		iterator_bin += 0.02 ;
	}	
	Ionration_vs_Ei     = new TH2D((histname_calib).c_str(), (histname_calib).c_str(),750, 0.,300.,100,0.,1500.);
	Ionration_vs_Ei     -> SetBins((int) binning_vec_kevee.size()-1,Binning_keVee,99 , ratio_binninh );	
	Dchi2_vs_Ep_pass    = new TH2D((histname+"_vs_Dchi2_pass").c_str(), (histname+"_vs_Dchi2_pass").c_str(),750, 0.,300.,100,0.,1500.);
	Dchi2_vs_Ep_pass    -> SetBins((int)binning_vec.size() -1, Binning_keV , 99, Binning_Dchi2);	
	Dchi2_vs_Ep_fail    = new TH2D((histname+"_vs_Dchi2_fail").c_str(), (histname+"_vs_Dchi2_fail").c_str(),750, 0.,300.,100,0.,1500.);
	Dchi2_vs_Ep_fail    -> SetBins((int)binning_vec.size() -1, Binning_keV , 99, Binning_Dchi2);		
	chi2_cut_vs_Ep_pass = new TH2D((histname+"_vs_chi2_pass").c_str(), (histname+"_vs_chi2_pass").c_str(),750, 0.,300.,100,0.,1500.);
	chi2_cut_vs_Ep_pass -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);
	chi2_cut_vs_Ep_fail = new TH2D((histname+"_vs_chi2_fail").c_str(), (histname+"_vs_chi2_fail").c_str(),750, 0.,300.,100,0.,1500.);
	chi2_cut_vs_Ep_fail -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);	
	Reso_cat_buffer = 0 ; 	
	int rejected_dchi2 = 0;
// Loop on both processed and calib , apply quality cuts
	for(int it = 0; it < Nb_HeatEnergy; it++, point_time_reso++ ){		
		chain_HeatEnergy           ->GetEntry(it);
		chain_chi2A                ->GetEntry(it);
		chain_event_processed      ->GetEntry(it);
		chain_event_processed_fast ->GetEntry(it);		
		Double_t Ep = Eh * (1 + (fabs(Voltage)/3.));
		H2_Eh_chi2->Fill(Ep,(chi2_A/1024.), 1./EpBinIndex(Ep, binning_vec));		
	       // if ( it%1000 == 0 )std::cout<< "absolute time  "<< Time_Crate/100000.  << " resolution "<< Reso_cat <<" entry "<<point_time_reso<<std::endl;
	        
		reso_vs_time->SetPoint(point_time_reso,Time_Crate / 100000., Reso_cat);		
		//chi2 cut
		if((chi2_A/1024.) > (1.15 + 100 * TMath::Power(fabs(Ep)/300. , 3.)) ){			 
			 chi2_cut_vs_Ep_fail-> Fill(Ep, (chi2_A/1024.), 1./EpBinIndex(Ep, binning_vec));
			 continue ;			 
		}
		chi2_cut_vs_Ep_pass -> Fill(Ep, (chi2_A/1024.), 1./EpBinIndex(Ep, binning_vec));		
		//delta chi2 cut
		if((chi2_norm[0]/1024.) - (chi2_fast[0]/1024.) > 0.001){			
			Dchi2_vs_Ep_fail -> Fill(Ep, (chi2_norm[0]/1024.) - (chi2_fast[0]/1024.), 1./EpBinIndex(Ep, binning_vec));
			rejected_dchi2++;
			continue;
		}
		Reso_cat_buffer += Reso_cat; 
		Dchi2_vs_Ep_pass -> Fill(Ep, (chi2_norm[0]/1024.) - (chi2_fast[0]/1024.), 1./EpBinIndex(Ep, binning_vec));
		H_Ehee           -> Fill(Eh, 1./EpBinIndex(Eh, binning_vec_kevee));
		H_Eh             -> Fill(Ep, 1./EpBinIndex(Ep, binning_vec) );
		H_Eh_lowres      -> Fill(Ep, 1./EpBinIndex(Ep, binning_vec_low_res) );
		Ionration_vs_Ei  -> Fill(fabs(Ei) , Eh / fabs(Ei) , 1./EpBinIndex(fabs(Ei), binning_vec_kevee));	
	} 	
	std::cout<<" delta chi2 rejected "<<rejected_dchi2<<std::endl;
	Reso_cat_buffer = Reso_cat_buffer/Nb_HeatEnergy;
	Double_t psd_freq [15] = {0.};
	Double_t psd_filt [15] = {0.};
	Double_t temp_max_filt = 0 ;
	for(int it = 0; it <  chain_voltage_pro->GetEntries() ; it++){	
		chain_voltage_pro->GetEntry(it);
		for(int it2 = 0; it2 <  15 ; it2++){
			psd_filt[it2] += std::pow (nVtoADU[0]* 1./(sqrt (1+ std::pow(cutofffreq/PDS_noise[it2][0],2*filter_order))),2);
			psd_freq[it2] = PDS_freq [it2] ;						
			if(it == chain_voltage_pro->GetEntries() - 1 ) {
				psd_filt[it2] = std::sqrt(psd_filt[it2]);				
				if(psd_filt[it2] > temp_max_filt )temp_max_filt = psd_filt[it2] ;
			}
		}		
	}	
	//PSD_plot_reso = new TGraphErrors(15,psd_freq ,psd_filt );
	//Reso_cat_buffer = temp_max_filt;
	std::cout<< " Max PSD noise : "<< Reso_cat_buffer <<std::endl;		
	Time_per_voltage = new TH1D ((histname+"_ellapsed_time").c_str(), (histname+"_ellapsed_time").c_str(),1,0.,1. );		
	Time_per_voltage->SetBinContent(1, Ellapsed_time );	
	std::cout<<" Integral for renormalization "<< Time_per_voltage -> Integral() <<std::endl;	
	Write_histo_tofile(heat, Voltage, Run_name);	
}

void Plot_HEnergy_Voltage::Loop_over_Chain_processed(){	
	chain_voltage_pro->GetEntry(0);
	std::cout<<" Voltage and Temp for run : "<< Run_name <<" "<<polarion[0] -polarion[2]<<" V "<<" "<< heat<<" mK "<<std::endl; 
	std::string voltname = "pos"+to_string(fabs(polarion[0] -polarion[2]));
	if(polarion[0] -polarion[2] < 0) voltname = "neg"+to_string(fabs( polarion[0] -polarion[2]));	
	std::string histname ="" ;
	Int_t Nbin_freq =500 ;
	if(allRUN == 0){ 	
		histname = "PSD_"+voltname+"_"+to_string(heat)+"mk" ;		
	}else{	
		histname = "PSD_"+voltname+"_"+to_string(heat)+"mk_"+Run_name ;	
	}	
	Double_t psd_filt [1024] = {0.};
	Double_t psd_freq [1024] = {0.};
	Double_t temp_freq_heat_max = 0 ;	
	for(int it = 0; it <  chain_voltage_pro->GetEntries() ; it ++){		
		chain_voltage_pro->GetEntry(it);	
		for(int it2 = 0; it2 <  1024 ; it2++){
			psd_filt[it2] += std::pow (nVtoADU[0]* 1./(sqrt (1+ std::pow(cutofffreq/PDS_noise[it2][0],2*filter_order))),2);
			psd_freq[it2] = PDS_freq [it2] ;			
			if(it = chain_voltage_pro->GetEntries() - 1 ) psd_filt[it2] = std::sqrt(psd_filt[it2]);
		}	
		temp_freq_heat_max = f_max_heat ;		
	}	
	PSD_plot = new TGraphErrors(1024,psd_freq, psd_filt );		
	Write_histo_tofile(heat, polarion[0] -polarion[2], Run_name);
}

void Plot_HEnergy_Voltage::SetTemp(){
	std::cout<<" type heat  : "<<std::endl; 
	std::cin>> heat ;	
}

void Plot_HEnergy_Voltage::SetRunname(){
	std::cout<<" type run name  : "<<std::endl; 
	std::cin>> Run_name ;
}

void Plot_HEnergy_Voltage::SetTemp(Double_t heat_){
	std::cout<<"  heat   = "<<heat_<<std::endl; 
	 heat = heat_;	
}

void Plot_HEnergy_Voltage::SetRunname(const std::string & runName_ ){
	std::cout<<" run name  : "<<(runName_).c_str()<<std::endl; 
	 Run_name = runName_;
}
void Plot_HEnergy_Voltage::Init(){
	//SetTemp();
	//SetRunname();	
	//Cryo_Run Run_317_reso_buffer ("reso_vs_time");
	point_time_reso = 0 ;
	reso_vs_time = new TGraph();
	Open_file(Run_name.c_str());	
	if(N_partition != 0) {	
	 if(IS_PROCESSED==0) {	 
	 	Loop_over_Chain();
	}else{
		Loop_over_Chain_processed();
	}	
	Write_timed_reso(("Resolution_timed_"+to_string(heat)+".root").c_str());	
	}else{	
		std::cout<<"No partition in the run"<<std::endl;
	} 
}

void Plot_HEnergy_Voltage::Write_timed_reso(const std::string & name_output_reso){

	TFile* output = new TFile(name_output_reso.c_str(), "RECREATE");
	reso_vs_time->Write("reso_vs_time");	
	output->Close();
	delete output ;
}

void Plot_HEnergy_Voltage::Estimate_Run_ellapsed_time(){

}

void Plot_HEnergy_Voltage::loop_over_generic_chain(TChain* chain){

}


void Plot_HEnergy_Voltage::cleaning(){


}
