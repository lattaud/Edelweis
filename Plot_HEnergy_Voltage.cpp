#include "Plot_HEnergy_Voltage.h"

using namespace std;

Plot_HEnergy_Voltage::Plot_HEnergy_Voltage(const std::string list_name_in , Double_t HEAT, bool IsRun, bool On_processed , std::string outputdir){
	

	TH1::SetDefaultSumw2();	
	IS_PROCESSED = On_processed;
	OutputDir = outputdir ;
	system(("Create_outputdir.sh "+OutputDir).c_str());
	list_name = list_name_in;
	if (IsRun){
		RunOnly(HEAT);
	}
	else{
		std::cout<<"enter number of Run in the list"<<std::endl;
		//cin>> Nlist;
		allRUN = 0 ;
		Parse_List();		

	}
	

}

void Plot_HEnergy_Voltage::RunOnly(Double_t HEAT ){
	TH1::SetDefaultSumw2();

	SetTemp(HEAT);
	SetRunname(list_name);
	std::cout<<"Create Ploting instance "<<std::endl;
	Init();
	
	

}

void Plot_HEnergy_Voltage::RunList(Double_t Heat, std::string list){
	TH1::SetDefaultSumw2();

	SetTemp(Heat);
	SetRunname(list);
	std::cout<<"Create Ploting instance "<<std::endl;
	Init();
	
	

}
Plot_HEnergy_Voltage::~Plot_HEnergy_Voltage() {

	delete chain_voltage;
	delete chain_index ;	
	delete chain_HeatEnergy ;
	delete chain_chi2A;
	delete H_Eh;

}

void Plot_HEnergy_Voltage::Help(){
	std::cout<<" Please use the following option in the comand line : list_name  "<<std::endl;
}
void Plot_HEnergy_Voltage::Parse_List(){

	double Heat_cat[100] ;
	
	
	std::string inputList = list_name;
	std::cout << " opening list : "<<inputList<<std::endl;
	ifstream Listfile(inputList.c_str(),ios::in);
	
	char pNamemc[500];
	std::string list_name_temp;
	double Heat_perRun; 
	int ilist = 0;
	int count_line_mc = 0;
	std::string temp_namelistIN = "";
	while(  !Listfile.eof() )
	{
		count_line_mc = 0;
               	Listfile>>list_name_temp>>Heat_perRun;
		Heat_cat[ilist]	= Heat_perRun;
	        std::cout<<" adding : /sps/edelweis/rootDataRun317/streams/prodg/lists/"<<list_name_temp<< ".list Temp "<< Heat_cat[ilist]<< std::endl;
	        string file = list_name_temp ;
	        if(file == temp_namelistIN) continue;
	        ///ifstream efficiencies(("/sps/edelweis/rootDataRun317/streams/prodg/lists/"+file+".list").c_str(),ios::in);
	        ifstream efficiencies((file+".list").c_str(),ios::in);
	        std::string ListRun_name[100];
	        while(!efficiencies.eof() ){

	       	 	efficiencies>>ListRun_name[count_line_mc];
			std::cout<<" List content "<< ListRun_name[count_line_mc]<<std::endl;
			system(("Create_list.sh "+ListRun_name[count_line_mc]).c_str());
			if (ListRun_name[count_line_mc] != " ") {

				std::cout<<" Run name ? " << ListRun_name[count_line_mc]<< " Heat  "<< Heat_cat[ilist] <<std::endl;						
				RunList( Heat_cat[ilist] , ListRun_name[count_line_mc] );
				
			}
	       	 	count_line_mc += 1; 
	        }
	        

		ilist++;
		temp_namelistIN = file ;
		
	}
	
}

void Plot_HEnergy_Voltage::Open_file( std::string file_name){


	//if(IS_PROCESSED==0){
	
	// everest output
	
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
	
		Nb_voltage    = chain_voltage    -> GetEntries();
		Nb_index      = chain_index      -> GetEntries();
		Nb_HeatEnergy = chain_HeatEnergy -> GetEntries();
		Nb_chi2       = chain_chi2A      -> GetEntries();
	//}else{
	
		//Processed stuff (Nepal output) 
		
		chain_voltage_pro            = new TChain("RunTree_Normal") ;
		chain_event_processed    = new TChain("EventTree_trig_Normal_filt_decor");
	
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
		        chain_voltage_pro     ->Add(pNamemc);
		        chain_event_processed ->Add(pNamemc);
		        
		}
		
		chain_voltage_pro   	  -> SetBranchAddress ("PSD_Filt",&PDS_noise);
		chain_voltage_pro   	  -> SetBranchAddress ("Polar_Ion",&polarion);
		chain_voltage_pro  	  -> SetBranchAddress ("cutoff_freq",&cutofffreq);
		chain_voltage_pro  	  -> SetBranchAddress ("filter_order",&filter_order);
		chain_voltage_pro  	   -> SetBranchAddress ("Chan_Gain",&nVtoADU);
		chain_voltage_pro   	  -> SetBranchAddress ("PSD_Freq",&PDS_freq);
		chain_voltage_pro  	  -> SetBranchAddress ("f_max_heat",&f_max_heat);
		
		chain_event_processed     -> SetBranchAddress ("MicroStp",&micro_step);
		
		Nb_voltage    = chain_voltage_pro    -> GetEntries();
		
	
	//}
	
	//if(Nb_chi2 != Nb_HeatEnergy) throw "Trees have different sized : CHECK REQUESTED";
	cout << "[+] Linking variable... done                           " << endl;
}
void Plot_HEnergy_Voltage::Write_histo_tofile(float temp, int voltage, std::string run_name){

	std::string reso_CAT;
	std::string Processed = "";
	if(Reso_cat_buffer < 1.) reso_CAT = "highres";
	if(Reso_cat_buffer < 2. && Reso_cat_buffer >= 1.) reso_CAT = "mediumres";
	if( Reso_cat_buffer >= 2.) reso_CAT = "lowres";
	if( IS_PROCESSED==1) Processed = "processed";
	
	if(allRUN == 1) Output_Files = new TFile((OutputDir+"/Eh_allruns_"+to_string(temp)+"mk.root").c_str(),"UPDATE");
	if(allRUN == 0) Output_Files = new TFile((OutputDir+"/Eh_perrun_"+to_string(temp)+"mk_"+run_name.c_str()+"_"+to_string(voltage)+"_"+reso_CAT+"_"+Processed+".root").c_str(),"RECREATE");
	
	TParameter<Double_t> * resoHEAT = new TParameter<Double_t>("Resolution_heat", Reso_cat_buffer);
	if( run_name != "" && IS_PROCESSED==0) {
		H_Ehee      ->Write();    
		H_Eh        ->Write();
		H2_Eh_chi2  ->Write();
		H_Eh_lowres  ->Write();
		resoHEAT     ->Write();
		
	}else{
	
		if( run_name != "") PSD_plot->Write("PSD_spectrum");	
	}
	
	
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
	
	
	
	Ellapsed_time = 1. / ( time_2 - time_1 );
	
	
	chain_index   ->GetEntry(0) ;
	chain_voltage ->GetEntry(Index_Calib);
	
	std::cout<<" Voltage and Temp for run : "<< Run_name <<" "<<Voltage<<" V "<<" "<< heat<<" mK run lasted for "<< Ellapsed_time / 3600.<< " h" <<std::endl; 
	std::string voltname = "pos"+to_string(fabs( Voltage));
	//if(Voltage < 0) voltname = "neg"+to_string(fabs( Voltage));
	std::string histname  ;
	std::string histname_ee  ;


	if(allRUN == 0){ 
	
		histname = "Ephonon_"+voltname+"_"+to_string(heat)+"mk" ;
		histname_ee = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_keVee" ;
		
	}else{
	
		histname = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_"+Run_name ;
		histname_ee = "Ephonon_"+voltname+"_"+to_string(heat)+"mk_"+Run_name+"_keVee" ;
	
	}
	std::cout<<" Creating  Histo "<<histname<<std::endl;
	
	
	Int_t Nbin_sub100   = 0;
	Int_t Nbin          = 0;
	



	//Double_t Binning_keVee[251] = {0.};
	
	Double_t Binning_chi2[1000] = {0.} ;
	
	Double_t itbin = 0.;
	Double_t itbinchi2 = 0.;
	for(int bin = 0 ; bin < 1000; bin ++ ){

			itbinchi2 += 0.1;
			Binning_chi2[bin] = itbinchi2;
		 

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

	}while(iterator_bin_ee <= 15);
	
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
	
	
	H_Eh   = new TH1D((histname).c_str(), (histname).c_str(),750., 0.,300.);
	H_Eh->SetBins((int)binning_vec.size()-1, Binning_keV);
	
	H_Eh_lowres   = new TH1D((histname+"_lowres").c_str(), (histname+"_lowres").c_str(),750., 0.,300.);
	H_Eh_lowres->SetBins((int)binning_vec_low_res.size()-1, Binning_keVlow);
	
	H_Ehee = new TH1D((histname_ee).c_str(), (histname_ee).c_str(),375, 0.,15.);	
	H_Ehee->SetBins((int) binning_vec_kevee.size()-1, Binning_keVee) ; 
	
	H2_Eh_chi2 = new TH2D((histname+"_vs_chi2").c_str(), (histname+"_vs_chi2").c_str(),750, 0.,300.,100,0.,1500.); 
	H2_Eh_chi2 -> SetBins((int)binning_vec.size() -1, Binning_keV , 999, Binning_chi2);
	
	
	Reso_cat_buffer = 0 ; 
	for(int it = 0; it < Nb_HeatEnergy; it++ ){
		
		chain_HeatEnergy   ->GetEntry(it);
		chain_chi2A        ->GetEntry(it);
		
		Double_t Ep = Eh * (1 + (fabs(Voltage)/3.));
		H2_Eh_chi2->Fill(Ep,(chi2_A/1024.));
		if((chi2_A/1024.) > (1.15 + 100 * TMath::Power(fabs(Ep)/300. , 3.)) )	 continue ;

		Reso_cat_buffer += Reso_cat; 

		H_Ehee       -> Fill(Eh, Ellapsed_time*EpBinIndex(Eh, binning_vec_kevee));
		H_Eh         -> Fill(Ep, Ellapsed_time*EpBinIndex(Ep, binning_vec) );
		H_Eh_lowres  ->Fill(Ep, Ellapsed_time*EpBinIndex(Ep, binning_vec_low_res) );
		//std::cout<<" testing hist weight "<< Kevee_weight(Eh) << "  " << EpBinIndex(Ep, binning_vec) <<std::endl;
		
	} 
	
	Reso_cat_buffer = Reso_cat_buffer/Nb_HeatEnergy;
	
	
	Write_histo_tofile(heat, Voltage, Run_name);
	 

}

void Plot_HEnergy_Voltage::Loop_over_Chain_processed(){
	
	chain_voltage_pro->GetEntry(0);
	std::cout<<" Voltage and Temp for run : "<< Run_name <<" "<<polarion[0] -polarion[2]<<" V "<<" "<< heat<<" mK "<<std::endl; 
	std::string voltname = "pos"+to_string(fabs(polarion[0] -polarion[2]));
	if(polarion[0] -polarion[2] < 0) voltname = "neg"+to_string(fabs( polarion[0] -polarion[2]));

	
	
	std::string histname  ;
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

void Plot_HEnergy_Voltage::SetRunname(std::string runName_){

	std::cout<<" run name  : "<<(runName_).c_str()<<std::endl; 
	 Run_name = runName_;


}
void Plot_HEnergy_Voltage::Init(){

	//SetTemp();
	//SetRunname();
	Open_file(Run_name.c_str());	
	 if(IS_PROCESSED==0) {	 
	 	Loop_over_Chain();
	}else{
		Loop_over_Chain_processed();
	}	
 
}
void Plot_HEnergy_Voltage::Estimate_Run_ellapsed_time(){





}

void Plot_HEnergy_Voltage::loop_over_generic_chain(TChain* chain){


}


void Plot_HEnergy_Voltage::cleaning(){


	for(int itHist = 0; itHist < count_line ; itHist ++ ){
		//delete H_Eh [itHist];
	}
}
