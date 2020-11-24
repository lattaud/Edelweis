#include "Smearing_spectrum.h"
using namespace std;

//GLOBAL Graph for TF1 estimation -> Ugly ugly fix
TGraph* Global_spectrum    = new TGraph();
TGraph* Global_efficiency  = new TGraph();
TGraph* Global_signal_Wimp = new TGraph();

double ProbCDMS(unsigned int k,double mu,double F, double Verbosity);
double Efficiency_from_graph(double *x,double *par);
double Spectrum_from_graph(double *x,double *par);
double Wimp_from_graph(double *x,double *par);
double Quantified_spectrum(double *x,double *par);
double Quantified_spectrum_slow(double *x,double *par);
double EffCorr_Quantified_spectrum(double *x,double *par);
double EffCorr_WimpDMN(double *x,double *par);
double Smeared_spectrumDMN(double *x,double *par);
double Smeared_spectrum(double *x,double *par);
void Smearing_distribution::Get_graph(){
   cout<<"Which mass to smear : "<<endl;
   cout<<Form("%sMeV_spectrum_%2.0fV",Mass_str.c_str(),Tension)<<endl;
   Spectrum_graph    = (TGraph*) Input_file->Get(Form("%sMeV_spectrum_51V",Mass_str.c_str(),Tension));
   
}
void Smearing_distribution::Get_wimpgr(){
    
    Wimp_graph      = (TGraph*) Wimp_file->Get("GeV5_spectrum");
    std::cout<<" Npoint "<<Wimp_graph->Integral()<<std::endl;
}

void Smearing_distribution::Get_efficiencies(){
   Efficiency_graph  = (TGraph*) Eff_file->Get(Form("Efficienciy_cutChihB_%s",Eff_NAME.c_str()));   
}
void Smearing_distribution::Apply_efficiencies(){
   
   
}      

void Smearing_distribution::Write_to_file(){

   if(Verbosity >= 1 )cout<<"Storing to file "<<Output_name<<endl;
   Efficiency_function = new TF1("eff_function",Efficiency_from_graph,0.01,1.5,1);
   Efficiency_function->SetParameter(0,Tension);
   Spectrum_function  = new TF1("Raw_spectrum_function",Spectrum_from_graph,0.,10.,0);
   TFile* Ouput_file = new TFile(Form("%s_%s_%s",Mass_str.c_str(),Eff_NAME.c_str(),Output_name.c_str()),"RECREATE");
   Smeared_spectrum_fct->SetNpx(5000);
   Smeared_spectrum_fct->Write();   
   Spectrum_function->Write();
   Quantified_spectrum_fct->Write();
//---------TEST ON WIMPS -----------------//   
   /*Smeared_Wimp_fct->SetNpx(5000);
   Smeared_Wimp_fct->Write();   
   Efficiency_function->Write();
   Effcorr_Wimp_fct->SetNpx(5000);
   Effcorr_Wimp_fct->Write();
   TF1*   WIMPfct =  new TF1("Raw_wimp",Wimp_from_graph,0.0,10.,0);
   WIMPfct->SetNpx(5000);
   WIMPfct->Write();*/
//---------END TEST ON WIMPS -----------------//   

   //Effcorr_Quantified_spectrum_fct->Write();
   /*TGraph* Rawspectrum             = new TGraph();
   TGraph* Qspectrum               = new TGraph();
   TGraph* Qspectrumeffcorr        = new TGraph();
   TGraph* QspectrumeffcorrSmeared = new TGraph();*/
   Ouput_file->Close();
}
Smearing_distribution::Smearing_distribution(const std::string input_file, const std::string &list_histo, const std::string output_name, const std::string efficiency_name, const std::string detector_name, double tension, double percent_running, double fano, double resolution, double verbosity, std::string mass_str,std::string const& EIcut){

   Output_name      = output_name;
   Detector         = detector_name;
   Input_file       = new TFile(input_file.c_str());  
   Eff_file         = new TFile(Form("Efficiencies_%s_%s_cutFid.root",EIcut.c_str(),efficiency_name.c_str()));
   Wimp_file        = new TFile("Spectrum_DMNR_phonon5_Gev.root");
   List_mass        = list_histo;
   Tension          = tension;
   Percent_running  = percent_running;
   Fano             = fano;
   Verbosity        = verbosity;
   Resolution       = resolution;
   Mass_str         = mass_str ;
   Eff_NAME         = efficiency_name ;
   Eicut            = EIcut;
}

void Smearing_distribution::Quantify_spectrum(){

  Quantified_spectrum_fct =  new TF1("spectrum_phonon_quantified",Quantified_spectrum_slow,0.03,(1+Tension/3.),2);
  Quantified_spectrum_fct->SetParameter(0,Tension);
  Quantified_spectrum_fct->SetParameter(1,Fano);
  
  Effcorr_Quantified_spectrum_fct = new TF1("spectrum_phonon_quantified_EffCorr", EffCorr_Quantified_spectrum ,0.03,(1+Tension/3.),2); 
  Effcorr_Quantified_spectrum_fct->SetParameter(0,Tension);
  Effcorr_Quantified_spectrum_fct->SetParameter(1,Fano);
  
  
  Effcorr_Wimp_fct = new TF1("spectrum_WIMP_EffCorr", EffCorr_WimpDMN ,0.03,(1+Tension/3.),1); 
  Effcorr_Wimp_fct->SetParameter(0,Tension);
}


double EffCorr_WimpDMN(double *x,double *par){
    double Ep = x[0];
    double V  = par[0];
    TF1*   spectrum =  new TF1("spectrum",Wimp_from_graph,0.0,3.,0);
	TF1*   efficiency =  new TF1("eff",Efficiency_from_graph,0.01,3.,1);
	efficiency->SetParameter(0,V);
	
	double Efficiency = efficiency->Eval(Ep/(1+V/3.)); 
	double dRDER = spectrum->Eval(Ep);
	if(Efficiency < 0) return 0.;
    return dRDER*Efficiency;
}


double EffCorr_Quantified_spectrum(double *x,double *par){
    double Ep = x[0];
    double V  = par[0];
	double F  = par[1];
	double verbosity = 0;
	TF1*   spectrum =  new TF1("spectrum",Spectrum_from_graph,0.0,1.,0);
	TF1*   efficiency =  new TF1("eff",Efficiency_from_graph,0.01,3.,1);
	efficiency->SetParameter(0,V);
	double trueErmax=200.;//spectrum->GetXmax();
	double sum=0;
	double Globalgapenergy     = 0.67;
	double Globalepsilongamma  = 3;
	double Efficiency = efficiency->Eval(Ep/(1+V/3.)); 
	if(Efficiency < 0) return 0.;
	if(verbosity > 0 ) cout<<"Ep "<<Ep*1000.<<" Eee "<<Ep*1000./(1+V/3.)<<endl;
	if(verbosity > 0 )cout<<" efficiency "<< Efficiency<<endl;
	for(int npair = 0 ; npair < 150; npair++){
	   if(Ep*1000. - (npair * V)  < 0 ){
	      //cout<<" Paire "<<npair<<" Er rejected "<<Ep*1000. - (npair * V)<<endl;
	      continue ;
	   }
	   double dRDER = spectrum->Eval(Ep - (npair * V)/1000. );
	   if(Ep - (npair * V)/1000. > 1.)  dRDER = 0 ;
	   sum += dRDER * ProbCDMS(npair,(Ep*1000 - npair * V)/3.,F,0);
	   if(verbosity > 1 ){     
	      cout<<"dR/dER "<<dRDER<<" E_{r} "<<Ep*1000. - (npair * V)<<" n "<<npair<<endl;
	      cout<<" proba "<<ProbCDMS(npair,Ep - (npair * V)/1000.,F,0)<<endl;
	   }
	}
	
   return sum*Efficiency;
   

}


double Quantified_spectrum_slow(double *x,double *par){
   double Ep = x[0];
   double V  = par[0];
	double F  = par[1];
	double verbosity = 0;
	TF1*   spectrum =  new TF1("spectrum",Spectrum_from_graph,0.0,1.,0);
	double trueErmax=200.;//spectrum->GetXmax();
	double sum=0;
	double Globalgapenergy     = 0.67;
	double Globalepsilongamma  = 3;
	if(verbosity > 0 ) cout<<"Ep "<<Ep*1000.<<endl;
	for(int npair = 0 ; npair < 150; npair++){
	   if(Ep*1000. - (npair * V)  < 0 ){
	      //cout<<" Paire "<<npair<<" Er rejected "<<Ep*1000. - (npair * V)<<endl;
	      continue ;
	   }
	   double dRDER = spectrum->Eval(Ep - (npair * V)/1000. );
	   if(Ep - (npair * V)/1000. > 1.)  dRDER = 0 ;
	   sum += dRDER * ProbCDMS(npair,(Ep*1000 - npair * V)/3.,F,0);
	   if(verbosity > 0 ){     
	      cout<<"dR/dER "<<dRDER<<" E_{r} "<<Ep*1000. - (npair * V)<<" n "<<npair<<endl;
	      cout<<" proba "<<ProbCDMS(npair,Ep - (npair * V)/1000.,F,0)<<endl;
	   }
	}
	
   return sum;
}

double ProbCDMS(unsigned int k,double mu,double F, double Verbosity)
{
	if(F>1 || F<0 || mu<0 ){cout<<"Can't do magic (For F>1 I could but I am lazy)"<<endl; exit(1);}	
	if(F==1){ 
		if(Verbosity >=1 ) cout<<"Poisson distribution"<<endl;		
		return ROOT::Math::poisson_pdf(k,mu);
	}
	if(mu==0){
		if(Verbosity >=1 ) cout<<" if mu=0, then necessarily the probability to get 0 is 1" << endl;
		if(k==0){return 1.;}else{return 0;}	
	}	
	if(mu>100){
		if(Verbosity >=1 ) cout<<" swtich to a Gaussian when mu>100"<<endl;
		return TMath::Gaus(k,mu,TMath::Sqrt(F*mu),kTRUE);
	}
	unsigned int nl=floor(mu/(1.-F));
   unsigned int nh=ceil(mu/(1.-F));		
   double Fl=1.-mu/nl;
	double Fh=1.-mu/nh;		
	if(Fl<=0 && Fh>0){		
		unsigned int Bernouilli_Nmax=ceil(mu);
		unsigned int Bernouilli_Nmin=floor(mu);
		double p=Bernouilli_Nmax-mu;
		double Proba_Bernouilli;
		double Fb=p*(1.-p)/mu;		
// Weigthed Binomial with Prob=1 of getting mu
		if(Bernouilli_Nmax==Bernouilli_Nmin){
			if(k==Bernouilli_Nmin){
				 Proba_Bernouilli=1;
				 Fb=0;
			}else{
				Proba_Bernouilli=0.;
			}			
// Weigthed Binomial and Bernouilli distribution	
		}else{
			if(k==Bernouilli_Nmin){
				 Proba_Bernouilli=p;
			}else if(k==Bernouilli_Nmax){
				Proba_Bernouilli=1.-p;
			}else{
				Proba_Bernouilli=0.;
			}		
		}				
		if(F<Fb){
			if(Verbosity >=1 ) cout<<"Bernouilli distribution"<<endl;
			return Proba_Bernouilli;		
		}else{
			if(Verbosity >=1 ) cout<<"Weighted Bernouilli and Binomial distributions"<<endl;
			double DeltaF=(F-Fb)/(Fh-Fb);
			double wb=(1.-DeltaF);
  			double wh=DeltaF;
			double Binomial_h=ROOT::Math::binomial_pdf(k,1.-Fh,nh);
			double Weighted=wb*Proba_Bernouilli+wh*Binomial_h; 			
			return Weighted;
		}		
	}	    	    	
    	if(nl==nh){ 
    		if(Verbosity >=1 ) cout<<"Binomial distribution"<<endl;
    		return ROOT::Math::binomial_pdf(k,1.-Fl,nl);
    	}else{
    		
    		if(Verbosity >=1 ) cout<<"Weighted double Binomial distribution"<<endl;
    		double DeltaF=(F-Fl)/(Fh-Fl);
    		double wl=(1.-DeltaF);
  		double wh=DeltaF;
	    	double Binomial_l=ROOT::Math::binomial_pdf(k,1.-Fl,nl);
	   	double Binomial_h=ROOT::Math::binomial_pdf(k,1.-Fh,nh);
	   	double Weighted=wl*Binomial_l+wh*Binomial_h; 	   	
	   	return Weighted;	
    	}
}
double maximum(double a,double b)
{
	if(a<=b){
	   return b;
	}else{   
	   return a;
	}
}
void Smearing_distribution::Smear(){
   Smeared_spectrum_fct = new TF1("smeared_fct",Smeared_spectrum,0.0,(1+Tension/3.),7);
   Smeared_spectrum_fct->SetParameter(0,Resolution);
   Smeared_spectrum_fct->SetParameter(1,Percent_running);
   Smeared_spectrum_fct->SetParameter(2,5);//Nsigma smearing gaussian
   Smeared_spectrum_fct->SetParameter(3,0.001);//Bin width for Riemman integration eV by eV
   Smeared_spectrum_fct->SetParameter(4,Verbosity); 
   Smeared_spectrum_fct->SetParameter(5,Fano);
   Smeared_spectrum_fct->SetParameter(6,Tension); 
   
   
   Smeared_Wimp_fct = new TF1("smeared_WIMP_fct",Smeared_spectrumDMN,0.0,(1+Tension/3.),6);
   Smeared_Wimp_fct->SetParameter(0,Resolution);
   Smeared_Wimp_fct->SetParameter(1,Percent_running);
   Smeared_Wimp_fct->SetParameter(2,5);//Nsigma smearing gaussian
   Smeared_Wimp_fct->SetParameter(3,0.001);//Bin width for Riemman integration eV by eV
   Smeared_Wimp_fct->SetParameter(4,Verbosity); 
   Smeared_Wimp_fct->SetParameter(5,Tension);
   
   
   
   if(Verbosity >=1 ){
      cout<< "DEBBUGGING "<<endl;
      cout<< "Spectrum evaluated in : 1.keV "<<Smeared_spectrum_fct->Eval(1.)<<endl;
   }
}
///////////////////////////////////NEEDED Function /////////////////////////////////////////////
double Smeared_spectrum(double *x,double *par){
    double Ep          = x[0] ;
	double sigmaphonon = par[0] ;
	double percent     = par[1]/100. ;
	double N           = par[2];
	double binwidth    = par[3];
	double verbosity   = par[4];
	double Fano        = par[5];
	double V           = par[6];
	TF1*   spectrum = new TF1("spectrum",EffCorr_Quantified_spectrum,0.0,3.,2);
	spectrum->SetParameter(0,V);
    spectrum->SetParameter(1,Fano);
	double sum = 0;
	double Riemman_iterator = Ep - N*sqrt(  pow(sigmaphonon,2) + pow(percent*Ep ,2) ) ;
	double Riemman_stop = Ep + N*sqrt(pow(sigmaphonon,2) + pow(percent*Ep ,2)) ;
	while(Riemman_iterator <= Riemman_stop){
	if(verbosity > 0 ){
	   cout<<"iterator "<< Riemman_iterator <<" stopping point "<<Riemman_stop <<endl;
	   cout<<"dR/dEpQ "<< spectrum->Eval(Riemman_iterator)<<" sigma tot "<<sqrt(pow(sigmaphonon,2)+pow(percent*Ep ,2))<<endl;
	}
	   sum += spectrum->Eval(Riemman_iterator) * (1./sqrt(2*M_PI*(pow(sigmaphonon,2) + pow(percent*Ep ,2))))*exp(-(pow(Ep-Riemman_iterator,2)/(2*pow(sigmaphonon,2) + pow(percent*Ep ,2))));
	   Riemman_iterator+=	binwidth ;	
	}
	
	
	return sum*binwidth  ;
}

double Smeared_spectrumDMN(double *x,double *par){

    double Ep          = x[0] ;
	double sigmaphonon = par[0] ;
	double percent     = par[1]/100. ;
	double N           = par[2];
	double binwidth    = par[3];
	double verbosity   = par[4];
	double V           = par[5];
	TF1*   spectrum = new TF1("spectrum",EffCorr_WimpDMN,0.0,10.,1);
    spectrum->SetParameter(0,V);
    double sum = 0;
	double Riemman_iterator = Ep - N*sqrt(  pow(sigmaphonon,2) + pow(percent*Ep ,2) ) ;
	double Riemman_stop = Ep + N*sqrt(pow(sigmaphonon,2) + pow(percent*Ep ,2)) ;
	while(Riemman_iterator <= Riemman_stop){
	if(verbosity > 0 ){
	   cout<<"iterator "<< Riemman_iterator <<" stopping point "<<Riemman_stop <<endl;
	   cout<<"dR/dEpQ "<< spectrum->Eval(Riemman_iterator)<<" sigma tot "<<sqrt(pow(sigmaphonon,2)+pow(percent*Ep ,2))<<endl;
	}
	   sum += spectrum->Eval(Riemman_iterator) * (1./sqrt(2*M_PI*(pow(sigmaphonon,2) + pow(percent*Ep ,2))))*exp(-(pow(Ep-Riemman_iterator,2)/(2*pow(sigmaphonon,2) + pow(percent*Ep ,2))));
	   Riemman_iterator+=	binwidth ;	
	}
	
	
	return sum*binwidth  ;

}

double Efficiency_from_graph(double *x,double *par)
{
    Double_t Voltage = par[0];
    Double_t Ehee = x[0];
	Double_t f =  Global_efficiency->Eval(Ehee *1000 ) ;
	return f;
}

double Spectrum_from_graph(double *x,double *par)
{  
   Double_t Ehee = x[0] ;
	Double_t f =  Global_spectrum->Eval(Ehee) ;
	return f;   
}

double Wimp_from_graph(double *x,double *par)
{  
    Double_t Ehee = x[0] ;
	Double_t f =  Global_signal_Wimp->Eval(Ehee*1000.) ;
	return f;   
}
//Get function :
TGraph* Smearing_distribution::Get_graph_eff(){
   return Efficiency_graph;
}
TGraph* Smearing_distribution::Get_graph_spectrum(){
   return Spectrum_graph;
}
TGraph* Smearing_distribution::Get_graph_WIMP(){
    return Wimp_graph;
}
int main(int argc, char** argv){  
   TCLAP::CmdLine cmd("Smearing signal", ' ', "0.1");
   TCLAP::ValueArg<std::string> inputfileList("", "inputfile", "input files", true, "", "string",cmd);
   TCLAP::ValueArg<std::string> inputplotList("", "inputplot-list", "Text file containing input graph", true, "", "string",cmd);
   TCLAP::ValueArg<std::string> Detector_name("d", "detector", "Which detector", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Output_name("o", "output-name", "Output name", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> Efficiency_name("", "efficiency-name", "Efficiency name", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> MASS_name("", "Mass", "Mass signal to smear", true, "", "string", cmd);
   TCLAP::ValueArg<std::string> EICUT("", "Eicut", "EICUT", true, "", "string", cmd);
   TCLAP::ValueArg<double> Voltage("", "Voltage", "Voltage Voltage", true, 0, "double", cmd);
   TCLAP::ValueArg<double> percent_resolution("", "Percent", "Percent resolution running", true, 0, "double", cmd);
   TCLAP::ValueArg<double> fano_arg("", "Fano", "Fano", true, 0, "double", cmd);
   TCLAP::ValueArg<double> reso_arg("", "Resolution_phonon", "Reso", true, 0, "double", cmd);
   TCLAP::ValueArg<double> verboseArg("", "v", "Verbosity level ? 0 mute ", true, 0, "double", cmd);

   cmd.parse(argc, argv);
   Smearing_distribution * Smearing_instance = new Smearing_distribution(inputfileList.getValue(),inputplotList.getValue(), Output_name.getValue(),Efficiency_name.getValue(),Detector_name.getValue(),Voltage.getValue(),percent_resolution.getValue(), fano_arg.getValue(), reso_arg.getValue() ,verboseArg.getValue(),MASS_name.getValue(),EICUT.getValue());
   Smearing_instance->Get_efficiencies();
   Smearing_instance->Get_graph();
   Smearing_instance->Get_wimpgr();
   Global_spectrum    = Smearing_instance->Get_graph_spectrum();
   Global_efficiency  = Smearing_instance->Get_graph_eff();
   Global_signal_Wimp = Smearing_instance->Get_graph_WIMP();
   Smearing_instance->Quantify_spectrum();
   Smearing_instance->Smear();
   Smearing_instance->Write_to_file();
   return 0;
}
