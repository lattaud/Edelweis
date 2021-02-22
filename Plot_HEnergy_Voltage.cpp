#include "Plot_HEnergy_Voltage.h"

using namespace std;

//------------------------------------------------------Constructor of the Class-----------------------------------------------------------------------------------
Plot_HEnergy_Voltage::Plot_HEnergy_Voltage(  std::string const & list_name_in, Double_t const & HEAT, bool const & IsRun, bool const & On_processed , std::string const & outputdir, bool const & local_list_, std::string const & detector_,unsigned int const & runonpairpart, Double_t const & ionCut, bool const & cut_ion_rej, std::string const &prod){
	
	TH1::SetDefaultSumw2();	
	IS_PROCESSED = On_processed;
	Detector = detector_;
	Neg_cutIon = cut_ion_rej ;
	OutputDir = Form("%s_%3.2feV_%s",outputdir.c_str(),ionCut,Detector.c_str()) ;
	Pair_partition = runonpairpart;
	IonCut = ionCut ;
	Prod = prod;
	std::cout<<" Running for detector "<< detector_<<" Partition "<<Pair_partition<<std::endl;
	if(detector_ == "RED30"){
	    weight_detector = 0.034; 
	    ndof_chi2 = 1024 ; 
	}else if(detector_ == "FID848"){
	    weight_detector = 0.87; 
	    ndof_chi2 = 1024 ; 
	}else if(detector_ == "NbSi209"){ 
	   weight_detector = 0.2;
	   ndof_chi2 = 512 ; 
	}
	system(("./Create_outputdir.sh "+OutputDir).c_str());
	list_name = list_name_in;
	local_list = local_list_ ;
	if (IsRun){
		RunOnly(HEAT);
	}
	else{
		std::cout<<"enter number of Run in the list"<<std::endl;
		allRUN = 0 ;
		Parse_List();	
	}
}

void Plot_HEnergy_Voltage::RunOnly(Double_t const & HEAT ){
	TH1::SetDefaultSumw2();
	SetTemp(HEAT);
	SetRunname(list_name);
	std::cout<<"Create Skimming instance "<<std::endl;
	Init();
}

void Plot_HEnergy_Voltage::RunList(Double_t const & Heat, std::string const & list){
	TH1::SetDefaultSumw2();
	SetTemp(Heat);
	SetRunname(list);
	std::cout<<"Create Skimming instance "<<std::endl;
	Init();	
}
void Plot_HEnergy_Voltage::Clean(){
//---------------- memory release----------------
	delete chain_voltage ;
	delete chain_index ;		
	delete chain_HeatEnergy ;
	delete chain_chi2A;
	delete chain_event_processed ;
	delete chain_voltage_pro ;
	//delete chain_event_processed_fast ;
	delete chain_event_Reso_processed ;
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
	delete time_lengh;
	delete H_Eh_noweight;
	delete H_EhB_noweight;
	delete H_Ehtot_noweight;
	delete H_Ehee_noweight;
	delete H_EhBee_noweight;
	delete H_Ehtotee_noweight;
	
}


void Plot_HEnergy_Voltage::Help(){
	std::cout<<" Please use the following option in the command line : list_name  "<<std::endl;
}
//--------------------------------------Chi 2 cut currently hard codded soon will be loaded from Jules G Calib files ----------------------------------------------
bool Plot_HEnergy_Voltage::Pass_chi2_cut(std::string const & detector, Double_t const & Ep, Double_t const & chi2A, Double_t const & chi2B){

   bool output = true;
//--------------------------------------Cut on heat channel A chi 2------------------------------------------------------------------------------------------------------
   double cutA = 1000000, cutB =1000000; 
   if(detector == "RED30") cutA = (1.15 + 100 * TMath::Power(fabs(Ep)/300. , 3.)) ;
   if(detector == "NbSi209") cutA = (1.2 + 100 * TMath::Power(fabs(Ep)/200. , 3.)) ;
   if(detector == "FID848") cutA = (1.15 + 100 * TMath::Power(fabs(Ep)/700. , 2.2)) ;
    
   output &=  chi2A <= cutA ;
//--------------------------------------Cut on heat channel B chi 2 if existing------------------------------------------------------------------------------------------   
   if(detector == "NbSi209") cutB = (1.2 + 100 * TMath::Power(fabs(Ep)/100. , 3.)) ;
   if(detector == "FID848") cutB = (1.15 + 100 * TMath::Power(fabs(Ep)/600. , 2.2)) ;
   
   if(detector != "RED30") output &= chi2B <= cutB ;
   return output;
}
//-------------------------------------Delta Chi 2 cut currently hard codded soon will be loaded from Jules G Calib files ----------------------------------------------
bool Plot_HEnergy_Voltage::Pass_Deltachi2_cut(std::string const & detector, Double_t const & chi2_1, Double_t  const & chi2_2){

   bool output = true;
   double cut = chi2_1 - chi2_2;
//-------------------------------------Delta chi2 cut only applied on RED30 so far----------------------------------------------------------------------------------
   if(detector == "RED30")   output &= cut < 0.001;
   if(detector == "NbSi209") return true; 
   if(detector == "FID848")  return true;
   
   
   return output;
}
//-------------------------------------List parser method ---------------------------------------------------------------------------------------------------------
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
	point_time_reso = 0 ;
	reso_vs_time = new TGraph();
	std::string prod = Prod;
	std::string script_parity_partition ;
	if(!Listfile.fail()){
		while(  std::getline(Listfile, line_list_IN) )
		{
		        TString temp_str_1 =line_list_IN ;
		        if(temp_str_1.Contains("#", TString::kIgnoreCase)) continue ;	
			    count_line_mc = 0;        	     
        	    unsigned int Index_space = line_list_IN.find(" ");
        	    std::string temp_str_heat = line_list_IN.substr(Index_space , 4);
        	    std::string temp_str_list = line_list_IN.substr(0, Index_space );
        	    Double_t temp_heat = std::stod(temp_str_heat);                
			    Heat_cat[ilist]	= temp_heat; 					
		        std::cout<<" adding : /sps/edelweis/rootDataRun317/streams/"+prod+"/lists/"<<temp_str_list<< ".list Temp "<< Heat_cat[ilist]<< std::endl;
		        string file = temp_str_list ;		       
		        std::string prefix_list = ""  ;		        
		        if(local_list == 0) prefix_list = "/sps/edelweis/rootDataRun317/streams/"+prod+"/lists/" ;		
		        std::cout<< prefix_list+file+".list"<<std::endl;       
		        ifstream efficiencies((prefix_list+file+".list").c_str());
		        std::string ListRun_name[100] = {""};
		        if(!efficiencies.fail()){
		       		 while(std::getline( efficiencies,ListRun_name[count_line_mc] )){
		       		        TString temp_str = ListRun_name[count_line_mc];
		       		        if(temp_str.Contains("#", TString::kIgnoreCase)) continue ;	
					           std::cout<<" List content "<< ListRun_name[count_line_mc]<<" Partition parity "<<Pair_partition<<std::endl;
					           if(Pair_partition == 2){
					            script_parity_partition ="Create_list.sh";
					           }else if(Pair_partition == 1){
					            script_parity_partition="Create_list_Odd.sh";
					           }else{
					            script_parity_partition="Create_listAll.sh";
					           }
					           std::cout<<" Create list for loop "<<script_parity_partition+" "+ListRun_name[count_line_mc]+" "+prod+" "+Detector<<std::endl;				
					           system(("./"+script_parity_partition+" "+ListRun_name[count_line_mc]+" "+prod+" "+Detector).c_str());					
					           std::cout<<" Run name ? " << ListRun_name[count_line_mc]<< " Heat  "<< Heat_cat[ilist] <<std::endl;						RunList( Heat_cat[ilist] , ListRun_name[count_line_mc] );				      
					           count_line_mc += 1; 
		       	 	}		        
			}
			    ilist++;		
		}
	}
	Write_timed_reso(("Resolution_timed_"+to_string(heat)+".root").c_str());	
}

void Plot_HEnergy_Voltage::Open_file( std::string const & file_name ){

//-------------------------------------Loading tree from Everest root file---------------------------------------------------------------------------------------------------

		chain_voltage              = new TChain("heatCalibData") ;
		chain_index                = new TChain("boloHeader");		
		chain_HeatEnergy           = new TChain("Energies_Trig_Filt_Decor") ;
		chain_chi2A                = new TChain("Amplitudes_Trig_Filt_Decor");
		chain_event_Reso_processed = new TChain("resolution_Noise_Filt_Decor");	
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
		      chain_voltage                   ->Add(pNamemc);
		      chain_index                     ->Add(pNamemc);
		      chain_HeatEnergy                ->Add(pNamemc);
			   chain_chi2A                     ->Add(pNamemc);
			   chain_event_Reso_processed      ->Add(pNamemc);
		}
//-------------------------------------Setting address of needed variable--------------------------------------------------------------------------------------------------------
		N_partition = count_line ;		
		cout<<N_partition << " partition in the Run " << endl;
		cout << "[+] Linking variable...                                 " << endl;
		chain_voltage              -> SetBranchAddress ("Voltage",&Voltage);
		chain_index                -> SetBranchAddress ("indexHeatCalibData",&Index_Calib);
		chain_HeatEnergy           -> SetBranchAddress ("Eh",&Eh);
		chain_HeatEnergy           -> SetBranchAddress ("EhA",&EhA);
		chain_HeatEnergy           -> SetBranchAddress ("EhB",&EhB);
		chain_chi2A                -> SetBranchAddress ("heatChi2A",&chi2_A);
		chain_chi2A                -> SetBranchAddress ("heatResoA",&Reso_cat);
		chain_HeatEnergy           -> SetBranchAddress ("EiB",&EiB);	
		chain_HeatEnergy           -> SetBranchAddress ("EiA",&EiA);
		chain_event_Reso_processed -> SetBranchAddress ("resoHeatA", &reso_eV); 
		Nb_voltage    = chain_voltage    -> GetEntries();
		Nb_index      = chain_index      -> GetEntries();
		Nb_HeatEnergy = chain_HeatEnergy -> GetEntries();
		Nb_chi2       = chain_chi2A      -> GetEntries();
	               
//-------------------------------------Loading tree from Processed root file (Nepal output-------------------------------------------------------------------------------------------------------- 
		
		chain_voltage_pro             = new TChain("RunTree_Normal") ;
		chain_event_processed         = new TChain("EventTree_trig_Normal_filt_decor");
		if(Detector == "RED30") chain_event_processed_fast    = new TChain("EventTree_trig_Fast_filt");
		if(Detector == "FID848" || Detector == "NbSi209"){
		   chain_event_processed_Slow    = new TChain("EventTree_trig_Slow_filt");
		   chain_event_processed_NTD    = new TChain("EventTree_trig_NTD_filt");
	}
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
		      if(Detector == "RED30")  chain_event_processed_fast  ->Add(pNamemc);
		      if(Detector == "FID848" || Detector == "NbSi209") chain_event_processed_Slow  ->Add(pNamemc);
		}
//-------------------------------------Setting address of needed variable--------------------------------------------------------------------------------------------------------		
		chain_voltage_pro   	   -> SetBranchAddress ("PSD_Filt",&PDS_noise);
		chain_voltage_pro   	   -> SetBranchAddress ("Polar_Ion",&polarion);
		chain_voltage_pro  	       -> SetBranchAddress ("cutoff_freq",&cutofffreq);
		chain_voltage_pro  	       -> SetBranchAddress ("filter_order",&filter_order);
		chain_voltage_pro  	       -> SetBranchAddress ("Chan_Gain",&nVtoADU);
		chain_voltage_pro   	   -> SetBranchAddress ("PSD_Freq",&PDS_freq);
		chain_voltage_pro  	       -> SetBranchAddress ("f_max_heat",&f_max_heat);		 
		chain_event_processed      -> SetBranchAddress ("MicroStp",&micro_step);
		chain_event_processed      -> SetBranchAddress ("Time_unix",&Time_Crate);
		chain_event_processed      -> SetBranchAddress ("Energy_OF",&Energy_OF);		
		chain_event_processed      -> SetBranchAddress ("chi2_OF_h",&chi2_norm);
		chain_event_processed      -> SetBranchAddress ("chi2_OF_i",&chi2_i);
		chain_event_processed      -> SetBranchAddress ("chi2_OF_half",&chi2_half);
		chain_event_processed      -> SetBranchAddress ("MegaStp",&Mega_stp);
		chain_event_processed      -> SetBranchAddress ("NumPart",&N_partitiontree);
	  if(Detector == "RED30") chain_event_processed_fast -> SetBranchAddress ("chi2_OF_h",&chi2_fast);
      if(Detector == "FID848" || Detector == "NbSi209"){
		   chain_event_processed_Slow -> SetBranchAddress ("chi2_OF_h",&chi2_fast);
	   }
		cout << "[+] Linking variable... done                           " << endl;
}
//-------------------------------------Method to write to output root file----------------------------------------------------------------------------------------------------------
void Plot_HEnergy_Voltage::Write_histo_tofile(float const & temp, int const & voltage, std::string const & run_name ){

	std::string reso_CAT="";
	std::string Processed = "";	
//-------------------------------------categorization based on ADU heat reso---------------------------------------------------------------------------------------------------------
	bool eV_reso = true ;
	double luke_boost = (1.+ fabs(voltage)/3.);
	luke_boost = 1. ;
	std::cout<<" LUKE BOOST "<<luke_boost<<std::endl;
	if(eV_reso ){
	    if(Detector == "NbSi209"){
	        if(Reso_cat_buffer < 0.09*luke_boost) reso_CAT = "highres";
	        if(Reso_cat_buffer >= 0.09*luke_boost && Reso_cat_buffer <= 0.12*luke_boost) reso_CAT = "mediumres";
	        if(Reso_cat_buffer >= 0.12*luke_boost) reso_CAT = "lowres";
	    }else if(Detector == "RED30"){
	        if(Reso_cat_buffer < 0.05*luke_boost) reso_CAT = "highres";
	        if(Reso_cat_buffer >= 0.05*luke_boost) reso_CAT = "mediumres";
	    }
	 }
//-------------------------------------categorization based on ADU heat reso--------------------------------------------------------------------------------------------------------------
	bool ADU_reso = false ;
	if(ADU_reso){
	        if(Reso_cat_buffer < 1.3) reso_CAT = "highres";
	        if(Reso_cat_buffer < 2. && Reso_cat_buffer >= 1.3) reso_CAT = "mediumres";
	        if( Reso_cat_buffer >= 2.) reso_CAT = "lowres";	
	 }
//-------------------------------------new categorisation -----------------------------------------------------------------------------------------------------------------------------------
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
//-------------------------------------Writting to output root file ------------------------------------------------------------------------------------------------------------------------------
	if( run_name != "" && IS_PROCESSED==0) {
		H_Ehee       ->Write();    
		H_Eh         ->Write();
		H_EhBee       ->Write();    
		H_EhB         ->Write();
		H_Ehtotee       ->Write();    
		H_Ehtot         ->Write();
		H2_Eh_chi2   ->Write();
		H_Eh_lowres  ->Write();
		resoHEAT     ->Write();
		Time_per_voltage->Write();
		Ionration_vs_Ei->Write();		
		Dchi2_vs_Ep_pass->Write();
		Dchi2_vs_Ep_fail->Write();	
		Dchi2Slow_vs_Ep_pass->Write();
		Dchi2Slow_vs_Ep_fail->Write();	
		Dchi2NTD_vs_Ep_pass->Write();
		Dchi2NTD_vs_Ep_fail->Write();		
		chi2_cut_vs_Ep_pass->Write();
		chi2_cut_vs_Ep_fail->Write();
		H_Eh_noweight->Write();
	    H_EhB_noweight->Write();
	    H_Ehtot_noweight->Write();
	    H_Ehee_noweight->Write();
	    H_EhBee_noweight->Write();
	    H_Ehtotee_noweight->Write();
		outTree_->Write();
		time_lengh->Write();
		G2_Eh_chi2->SetLineWidth(0);
		G2_Eh_chi2->SetMarkerStyle(8);
		G2_Eh_chi2->SetMarkerSize(0.5);
		H_EiFid          -> Write();
		H_EiA            -> Write();
		H_EiB            -> Write();
		EiFid_vs_Eh_passcut->SetLineWidth(0);
		EiFid_vs_Eh_rejected->SetLineWidth(0);
		EiFid_vs_Eh_rejected->SetMarkerColor(kRed);
		
		EiFid_vs_chi2_passcut->SetLineWidth(0);
		EiFid_vs_chi2_passcut->SetMarkerStyle(8);
		EiFid_vs_chi2_rejected->SetLineWidth(0);
		EiFid_vs_chi2_rejected->SetMarkerColor(kRed);
		EiFid_vs_chi2_passcut  -> Write();
		EiFid_vs_chi2_rejected -> Write();
        EiFid_vs_Eh_passcut    -> Write();
        EiFid_vs_Eh_rejected   -> Write();
		std::string signTension = "pos";
		if (voltage < 0.) signTension = "neg";
		G2_Eh_chi2         ->Write(("Graphchi2vsEP_"+signTension+to_string(fabs(voltage))).c_str());
		
	}else{
	
		if( run_name != "") PSD_plot->Write("PSD_spectrum");	
	}	
	delete resoHEAT;
	Output_Files->Close();
}
//------------------------------------- Method to get the size of the desired Energy bin---------------------------------------------------------------------------------------------
Double_t Plot_HEnergy_Voltage::EpBinIndex(Double_t const & pt,std::vector<Double_t> const & binning) {

	for (unsigned int i(1); i <  binning.size() ;  i++) {
		if (pt < binning.at(i) && pt > binning.at(i-1)) return binning.at(i) - binning.at(i-1) ;
		
	}
	return binning_vec.size();
}
//-------------------------------------Obsolete method not used anymore---------------------------------------------------------------------------------------------------------------
Double_t Plot_HEnergy_Voltage::Kevee_weight(Double_t const &  Eh){

	if(Eh < 1) return 0.01 ;
	if(Eh < 5) return 0.04 ;
	return 0.2;
}


void Plot_HEnergy_Voltage::Loop_over_Chain(){

//------------------------------------- Calculate time spent in whole run--------------------------------------------------------------------------------------------------------------
	Double_t temp_freq_heat_max = 0 ;
	chain_voltage_pro->GetEntry(3);
	temp_freq_heat_max = f_max_heat ;	
	Double_t Ellapsed_time = 0;	
	chain_event_processed->GetEntry(0);	
	Double_t time_1 = micro_step / temp_freq_heat_max ;	
	chain_event_processed->GetEntry(chain_event_processed->GetEntries() - 1 );	
	Double_t time_2 = (micro_step / temp_freq_heat_max) +  3600 * (N_partition - 1)  ;	
	Ellapsed_time = ( time_2 - time_1 );	
// //-------------------------------------Test time calculus alternatives ----------------
	chain_event_processed->GetEntry(0);
	Double_t altime1 = 3600. - micro_step / temp_freq_heat_max ;
	chain_event_processed->GetEntry(chain_event_processed->GetEntries() - 1 );
	Double_t altime2 = altime1 + ((micro_step / temp_freq_heat_max)  + (3600. *( N_partition - 2))) ;
//------------------------------------- End time spent in whole run	calculus.----------------
	chain_index   ->GetEntry(0) ;
	chain_voltage ->GetEntry(Index_Calib);	
	std::cout<<" Voltage and Temp for run : "<< Run_name <<" "<<Voltage<<" V "<<" "<< heat<<" mK run lasted for "<< (Ellapsed_time) / 3600.<< " h" <<std::endl; 
	std::string voltname = "pos"+to_string(int(fabs( Voltage)));
	std::string histname         = "";
	std::string histname_ee      = "";
	std::string histname_calib   = "";
	std::string histname_ion     = "";
	std::string graphname        = "";
	if(allRUN == 0){ 	
		histname       = "Ephonon_"+voltname ;
		histname_ee    = "Ephonon_"+voltname+"_keVee" ;
		histname_calib = "Ioniratio_vs_Ei_"+voltname ;
		histname_ion   = "Eion_"+voltname ;	
		graphname      = "EiFid_vs_Eh_"+voltname;	
	}else{	
		histname       = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_"+Run_name ;
		histname_ee    = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_"+Run_name+"_keVee" ;
		histname_calib = "Ioniratio_vs_Ei_"+voltname+"_"+to_string(heat)+"mk" ;	
	}
	std::cout<<" Creating  Histo "<<histname<<std::endl;	
	
//-------------------------------------Creating a binning for Chi2 studies------------------------------------------------------------------------------
	Int_t Nbin_sub100   = 0;
	Int_t Nbin          = 0;
	Double_t Binning_chi2[10000] = {0.} ;	
	Double_t itbin = 0.;
	Double_t itbinchi2 = 0.;
	for(int bin = 0 ; bin < 10000; bin ++ ){
			itbinchi2 += 0.1;
			Binning_chi2[bin] = itbinchi2;		
			Nbin++;
	}	
	Double_t Binning_Dchi2[100] = {0.} ;	
	itbinchi2 = -0.05;
	for(int bin = 0 ; bin < 100; bin ++ ){
			itbinchi2 += 0.001;
			Binning_Dchi2[bin] = itbinchi2;		 
	}		
	
//-------------------------------------Creating binning for heat spectrum in ee studies-----------------------------------------------------------------
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
//-------------------------------------Creating binning for heat spectrum in kev studies------------------------------------------------------------------------
	binning_vec.clear();	
	binning_vec.push_back(0.);
	Double_t iterator_bin = 0 ;
	int sumbin = 0 ; 
	do{
		if(iterator_bin == 0){		
			 iterator_bin += 0.04;	
			 sumbin++	;
		}else{
			Double_t sigma =std::sqrt( std::pow(0.02,2)  + std::pow(0.02*iterator_bin,2));
			iterator_bin += sigma ;
			sumbin++;
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
//-------------------------------------Histo declaration--------------------------------------------------------------------------------------------------------------		
	H_Eh          = new TH1D((histname).c_str(), (histname).c_str(),750., 0.,300.);
	H_Eh          ->SetBins((int)binning_vec.size()-1, Binning_keV);	
	H_Eh_lowres   = new TH1D((histname+"_lowres").c_str(), (histname+"_lowres").c_str(),750., 0.,300.);
	H_Eh_lowres   ->SetBins((int)binning_vec_low_res.size()-1, Binning_keVlow);	
	H_Ehee        = new TH1D((histname_ee).c_str(), (histname_ee).c_str(),375, 0.,15.);	
	H_Ehee        ->SetBins((int) binning_vec_kevee.size()-1, Binning_keVee) ; 	
	H2_Eh_chi2    = new TH2D((histname+"_vs_chi2").c_str(), (histname+"_vs_chi2").c_str(),750, 0.,300.,100,0.,1500.); 
	H2_Eh_chi2    -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);
	G2_Eh_chi2    = new TGraph();
	
	H2_Ei_chi2    = new TH2D((histname+"i_vs_chi2").c_str(), (histname+"_vs_chi2").c_str(),750, 0.,300.,100,0.,1500.); 
	H2_Ei_chi2    -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);
	
	H_EhB  = new TH1D((histname+"_B").c_str(), (histname).c_str(),750., 0.,300.);
	H_EhB  -> SetBins((int)binning_vec.size()-1, Binning_keV);	
    H_Ehtot = new TH1D((histname+"_tot").c_str(), (histname).c_str(),750., 0.,300.);
	H_Ehtot ->SetBins((int)binning_vec.size()-1, Binning_keV);	
	H_EhBee  = new TH1D((histname_ee+"_B").c_str(), (histname).c_str(),750., 0.,300.);
	H_EhBee  ->SetBins((int) binning_vec_kevee.size()-1, Binning_keVee) ; 	
    H_Ehtotee = new TH1D((histname_ee+"_tot").c_str(), (histname).c_str(),750., 0.,300.);
	H_Ehtotee ->SetBins((int) binning_vec_kevee.size()-1, Binning_keVee) ; 		
	H_Eh_noweight      = new TH1D((histname+"_noweight").c_str(), (histname).c_str(),5000., 0.,500.);
	H_EhB_noweight     = new TH1D((histname+"_B_noweight").c_str(), (histname).c_str(),5000., 0.,500.);
	H_Ehtot_noweight   = new TH1D((histname+"_tot_noweight").c_str(), (histname).c_str(),5000., 0.,500.);
	H_Ehee_noweight    = new TH1D((histname_ee+"_noweight").c_str(), (histname_ee).c_str(),500, 0.,150.);	
	H_EhBee_noweight   = new TH1D((histname_ee+"_B_noweight").c_str(), (histname).c_str(),500., 0.,150.);
	H_Ehtotee_noweight = new TH1D((histname_ee+"_tot_noweight").c_str(), (histname).c_str(),500., 0.,150.);
	
	H_EiA              = new TH1D((histname_ion+"_A").c_str(), (histname_ion+"_A").c_str(),1500., 0.,150.);
    H_EiB              = new TH1D((histname_ion+"_B").c_str(), (histname_ion+"_B").c_str(),1500., 0.,150.);
    H_EiFid            = new TH1D((histname_ion+"_Fid").c_str(), (histname_ion+"_Fid").c_str(),15000., 0.,150.);
    EiFid_vs_Eh_passcut  = new TGraph();
    EiFid_vs_Eh_passcut->SetName(graphname.c_str());
	EiFid_vs_Eh_rejected = new TGraph();	
    EiFid_vs_Eh_rejected->SetName(Form("%s_rejected",graphname.c_str()));
    
    EiFid_vs_chi2_passcut  = new TGraph();
    EiFid_vs_chi2_rejected = new TGraph();
    
    EiFid_vs_chi2_passcut  ->SetName(Form("%s_vschi2noChal",graphname.c_str()));
    EiFid_vs_chi2_rejected ->SetName(Form("%s_vschi2_all",graphname.c_str()));
	
	int Ngraph_point = 0 ;	
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
	Dchi2Slow_vs_Ep_pass    = new TH2D((histname+"_vs_Dchi2Slow_pass").c_str(), (histname+"_vs_Dchi2Slow_pass").c_str(),750, 0.,300.,100,0.,1500.);
	Dchi2Slow_vs_Ep_pass    -> SetBins((int)binning_vec.size() -1, Binning_keV , 99, Binning_Dchi2);	
	Dchi2Slow_vs_Ep_fail    = new TH2D((histname+"_vs_Dchi2Slow_fail").c_str(), (histname+"_vs_Dchi2Slow_fail").c_str(),750, 0.,300.,100,0.,1500.);
	Dchi2Slow_vs_Ep_fail    -> SetBins((int)binning_vec.size() -1, Binning_keV , 99, Binning_Dchi2);	
	Dchi2NTD_vs_Ep_pass    = new TH2D((histname+"_vs_Dchi2NTD_pass").c_str(), (histname+"_vs_Dchi2NTD_pass").c_str(),750, 0.,300.,100,0.,1500.);
	Dchi2NTD_vs_Ep_pass    -> SetBins((int)binning_vec.size() -1, Binning_keV , 99, Binning_Dchi2);	
	Dchi2NTD_vs_Ep_fail    = new TH2D((histname+"_vs_Dchi2NTD_fail").c_str(), (histname+"_vs_Dchi2NTD_fail").c_str(),750, 0.,300.,100,0.,1500.);
	Dchi2NTD_vs_Ep_fail    -> SetBins((int)binning_vec.size() -1, Binning_keV , 99, Binning_Dchi2);		
	chi2_cut_vs_Ep_pass = new TH2D((histname+"_vs_chi2_pass").c_str(), (histname+"_vs_chi2_pass").c_str(),750, 0.,300.,100,0.,1500.);
	chi2_cut_vs_Ep_pass -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);
	chi2_cut_vs_Ep_fail = new TH2D((histname+"_vs_chi2_fail").c_str(), (histname+"_vs_chi2_fail").c_str(),750, 0.,300.,100,0.,1500.);
	chi2_cut_vs_Ep_fail -> SetBins((int)binning_vec.size() -1, Binning_keV , 9999, Binning_chi2);		
	Reso_cat_buffer           = 0 ; 	
	int rejected_dchi2        = 0;
	int rejected_dchi2slow    = 0;
	int rejected_dchi2NTD     = 0;
	int buffer_time           = 0;
	int Added_time            = 0 ;
//-------------------------------------Loop on both processed and calib , apply quality cuts----------------------------------------------------------------------------
    int total_entries     = Nb_HeatEnergy;
    int cutchiion_entries = 0;
    int cutEion_entries   = 0;
    int cutchiA_entries   = 0;
    int cutFid_entries    = 0;
	for(int it = 0; it < Nb_HeatEnergy; it++  ){		
		chain_HeatEnergy           ->GetEntry(it);
		chain_chi2A                ->GetEntry(it);
		chain_event_processed      ->GetEntry(it);
		if(Detector == "RED30") chain_event_processed_fast ->GetEntry(it);
		if(Detector == "FID848" || Detector == "NbSi209") chain_event_processed_Slow ->GetEntry(it);		
		Double_t chi2A, chi2B , chi2ionA, chi2ionB ;
		if( Detector == "NbSi209"){
		   chi2A = chi2_norm[0];
		   chi2B = chi2_half[1];
		}else{
		   chi2A = chi2_norm[0];
		   chi2B = chi2_norm[1];
		}
		double Ei = (EiA + EiB) /2.;
//-------------------------------------TIME COMPUTATION--------------------------		
//-------------------------------------this method miss the last partition -> it is added after the loop onto events.	Not the final stored value -------------------------
		bool verbose = false;	
		if( Mega_stp <= 10 && buffer_time > 3500. ){		
		    Added_time += buffer_time;
		     if(verbose == true){
		        std::cout<<" increasing time "<<	Added_time <<std::endl;	
		        std::cout<<" buffer "<<	buffer_time <<" mega step "<<Mega_stp <<std::endl;	
		        std::cout<<" Partition "<<N_partitiontree<<std::endl;
		      }
		}
		buffer_time = Mega_stp ;
//-------------------------------------END TIME COMPUTATION------------------------------------------------------------------------------------------------------------------
		chi2ionA = chi2_norm[2];
		chi2ionB = chi2_norm[3];	
		Double_t EpA = EhA * (1 + (fabs(Voltage)/3.));
		Double_t EpB = EhB * (1 + (fabs(Voltage)/3.));
		Double_t Eptot = Eh * (1 + (fabs(Voltage)/3.));
		H2_Eh_chi2->Fill(Eptot,(chi2A/ndof_chi2), 1./EpBinIndex(Eptot, binning_vec));		
        G2_Eh_chi2->SetPoint(Ngraph_point, Eptot,(chi2_A/ndof_chi2));
        Ngraph_point++;
		
//-------------------------------------Cut on chi2 ionization-------------------------------------------------------------------------------------------------------------------------
		if((chi2ionA / 1024. ) > 1.8  && (chi2ionB / 1024. ) > 1.8){
		    cutchiion_entries++;
		    EiFid_vs_Eh_rejected->SetPoint(it,Eh, Ei);
		    continue;
		}
		
//-------------------------------------Cut on Ion Energy channel ---------------------------------------------------------------------------------------------------------------------------
		if(Ei < IonCut && IonCut >= 0 && !Neg_cutIon) {
		    cutEion_entries++;
		    EiFid_vs_Eh_rejected->SetPoint(it,Eh, Ei);
		    continue ;
		}else if(Ei >= IonCut && Neg_cutIon){
		    cutEion_entries++;
		    continue;
		}
//-------------------------------------Cut on NbSi fiducial volume Heat channel only----------------------------------------------------------------------------------------------------------------
		if(Detector == "NbSi209"  && pow(EpA - EpB,2) >= pow(0.6,2) + pow(0.045*Eptot,2) ) {
		    cutFid_entries++;
		    EiFid_vs_Eh_rejected->SetPoint(it,Eh, Ei);
		    continue ;
		}
//-------------------------------------Cut on chi2 heat ---------------------------------------------------------------------------------------------------------------------------------------
		if(!Pass_chi2_cut(Detector, Eptot, chi2A/1024., chi2B/ndof_chi2) ){			 
			 chi2_cut_vs_Ep_fail-> Fill(Eptot, (chi2A/ndof_chi2), 1./EpBinIndex(Eptot, binning_vec));
			 cutchiA_entries++;
			 EiFid_vs_Eh_rejected->SetPoint(it,Eh, Ei);
			 continue ;			 
		}
		chi2_cut_vs_Ep_pass -> Fill(Eptot, (chi2_A/ndof_chi2), 1./EpBinIndex(Eptot, binning_vec));		
		
//-------------------------------------Cut on delta chi2 -----------------------------------------------------------------------------------------------------------------------------------------
		if(!Pass_Deltachi2_cut(Detector,chi2A/ndof_chi2, (chi2_fast[0]/ndof_chi2) ) ){			
			Dchi2_vs_Ep_fail -> Fill(Eptot, (chi2A/ndof_chi2) - (chi2_fast[0]/ndof_chi2), 1./EpBinIndex(Eptot, binning_vec));
			rejected_dchi2++;
			continue;
		}
		
		Reso_cat_buffer += Reso_cat; 
		Dchi2_vs_Ep_pass -> Fill(EpA, (chi2_A/ndof_chi2) - (chi2_fast[0]/ndof_chi2), 1./EpBinIndex(EpA, binning_vec));
		H_Ehee           -> Fill(EhA, 1./(EpBinIndex(EhA, binning_vec_kevee)));
		H_Eh             -> Fill(EpA, 1./(EpBinIndex(EpA, binning_vec)));
		H_EhBee          -> Fill(EhB, 1./(EpBinIndex(EhB, binning_vec_kevee)));
		H_EhB            -> Fill(EpB, 1./(EpBinIndex(EpB, binning_vec)));
		H_Ehtotee        -> Fill(Eh, 1./(EpBinIndex(Eh, binning_vec_kevee)));
		H_Ehtot          -> Fill(Eptot, 1./(EpBinIndex(Eptot, binning_vec)));
		H_Eh_lowres      -> Fill(EpA, 1./(EpBinIndex(EpA, binning_vec_low_res)));
		H_EiFid          -> Fill(Ei);
		H_EiA            -> Fill(EiA);
		H_EiB            -> Fill(EiB);
		
	   if( Eh < 0.04)  EiFid_vs_chi2_passcut  ->SetPoint(it,Ei, chi2A/1024.);
       if( Eh > 0.04)  EiFid_vs_chi2_rejected ->SetPoint(it,Ei, chi2A/1024.);
		
	   H_Eh_noweight      -> Fill(EpA)  ;
	   H_EhB_noweight     -> Fill(EpB)  ;
	   H_Ehtot_noweight   -> Fill(Eptot,1./(0.1));
	   H_Ehee_noweight    -> Fill(EhA)  ; 
	   H_EhBee_noweight   -> Fill(EhB)  ; 
	   H_Ehtotee_noweight -> Fill(Eh);		
	   EiFid_vs_Eh_passcut-> SetPoint(it,Eh, Ei);
	   Ionration_vs_Ei    -> Fill(fabs(Ei) , Eh / fabs(Ei) , 1./(EpBinIndex(fabs(Ei), binning_vec_kevee)));
	   E_h_buf = Eh ;
	   E_p_buf = Eptot ;
	   weight  = (Ellapsed_time/(3600.*24.));
	   voltage = Voltage ;

		outTree_->Fill();
	}
//-------------------------------------END EVENT LOOP---------------------------------------------------------------------------------------------------------------------------------------
	chain_event_processed      ->GetEntry(chain_event_processed      ->GetEntries()-1);
	Added_time  += Mega_stp ;
	time_lengh = new TParameter<double> (("timed_lenght_"+voltname+"V").c_str(),Added_time/(3600.*24.));	 
	std::cout<<"*********************************TEST TIME ***********************************"	<<std::endl;
	std::cout<<"*********************OLD METHODS "<<Ellapsed_time/(3600.)<<std::endl;
	std::cout<<"*********************NEW METHODS "<<Added_time/(3600.)<<std::endl;
	std::cout<<"*********************************CUTFLOW***********************************"	<<std::endl;
	std::cout<<"*********************CHI2 ION    "<<(float(cutchiion_entries) / float(total_entries))*100.<<std::endl;	
	std::cout<<"*********************E ION       "<<(float(cutEion_entries)   / float(total_entries))*100.<<std::endl;		
	std::cout<<"*********************CHI2 A      "<<(float(cutchiA_entries)   / float(total_entries))*100.<<std::endl;
	std::cout<<"*********************Fid         "<<(float(cutFid_entries)    / float(total_entries))*100.<<std::endl;		
		
	Reso_cat_buffer = Reso_cat_buffer/Nb_HeatEnergy;
	Double_t psd_freq [15] = {0.};
	Double_t psd_filt [15] = {0.};
	Double_t temp_max_filt = 0 ;
//------------------------------------- PSD related compuation----------------------------------------------------------------------------------------------------------------------------
	for(int it = 0; it <  chain_voltage_pro->GetEntries() ; it++){	
		chain_voltage_pro->GetEntry(it);
		chain_event_Reso_processed ->GetEntry(it);
		
		Reso_cat_buffer = pow(reso_eV,2) ;
		
		for(int it2 = 0; it2 <  15 ; it2++){
			psd_filt[it2] += std::pow (nVtoADU[0]* 1./(sqrt (1+ std::pow(cutofffreq/PDS_noise[it2][0],2*filter_order))),2);
			psd_freq[it2] = PDS_freq [it2] ;						
			if(it == chain_voltage_pro->GetEntries() - 1 ) {
				psd_filt[it2] = std::sqrt(psd_filt[it2]);				
				if(psd_filt[it2] > temp_max_filt )temp_max_filt = psd_filt[it2] ;
			}
		}		
	}	
	Reso_cat_buffer = sqrt(Reso_cat_buffer);
	point_time_reso++;
	std::cout<< " Max PSD noise : "<< Reso_cat_buffer <<std::endl;		
	Time_per_voltage = new TH1D ((histname+"_ellapsed_time").c_str(), (histname+"_ellapsed_time").c_str(),1,0.,1. );		
	Time_per_voltage->SetBinContent(1, Added_time);	
	std::cout<<" Integral for renormalization "<< Time_per_voltage -> Integral() <<std::endl;	
	Write_histo_tofile(heat, Voltage, Run_name);	
}
//-------------------------------------Method to loop only on processed events-------------------------------------------------------------------------------------------------------------
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

void Plot_HEnergy_Voltage::SetTemp(Double_t const & heat_){
	std::cout<<"  heat   = "<<heat_<<std::endl; 
	 heat = heat_;	
}

void Plot_HEnergy_Voltage::SetRunname(std::string const & runName_ ){
	std::cout<<" run name  : "<<(runName_).c_str()<<std::endl; 
	 Run_name = runName_;
}
void Plot_HEnergy_Voltage::Init(){	
	Open_file(Run_name.c_str());	
	if(N_partition != 0) {	
	    if(IS_PROCESSED==0) {
	      outTree_ = new TTree("selected_events","selected_events");
	      outTree_->Branch("Eh",&E_h_buf,"Eh/D");
	      outTree_->Branch("Ep",&E_p_buf,"Ep/D");
	      outTree_->Branch("weight", &weight, "weight/D");
	      outTree_->Branch("voltage", &voltage, "voltage/D");		 
	      Loop_over_Chain();
	    }else{
		  Loop_over_Chain_processed();
	    }			
	}else{	
		std::cout<<"No partition in the run"<<std::endl;
	} 
}

void Plot_HEnergy_Voltage::Write_timed_reso( std::string const & name_output_reso){
	
	TFile* output = new TFile(name_output_reso.c_str(), "RECREATE");
	reso_vs_time->SetMarkerStyle(2);
	reso_vs_time->SetLineWidth(0);
	reso_vs_time->Write("reso_vs_time");	
	output->Close();
	delete output ;
}

//-----------------------UNDER DEV-----------------------
void Plot_HEnergy_Voltage::Estimate_Run_ellapsed_time(){

}

void Plot_HEnergy_Voltage::loop_over_generic_chain(TChain* chain){

}


void Plot_HEnergy_Voltage::cleaning(){


}
