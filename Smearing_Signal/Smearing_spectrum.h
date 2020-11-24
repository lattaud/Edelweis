#ifndef __SMEARING_SPECTRUM__
#define __SMEARING_SPECTRUM__

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
#include "TMath.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TParameter.h"
#include "../Smoothing_procedure/tclap/CmdLine.h"
#include "TKDE.h"
#include "TParameter.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Math.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/DistFunc.h"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <math.h>
using namespace std;

class Smearing_distribution {

   public:
      void Get_graph               ();
      void Get_wimpgr              ();
      void Get_efficiencies        ();
      void Apply_efficiencies      ();      
      void Write_to_file           ();
      
      Smearing_distribution        (const std::string input_file, const std::string &list_histo, const std::string Output_name, const std::string efficiency_name, const std::string detector_name, double tension, double percent_running, double fano, double resolution, double verbosity, std::string mass_str,std::string const& EIcut);
      ~Smearing_distribution()     = default;
      void    Smear                 ();
      void    Quantify_spectrum     ();
      TGraph* Get_graph_eff         (); 
      TGraph* Get_graph_spectrum    ();
      TGraph* Get_graph_WIMP        ();
   private: 
      TFile *                  Input_file;
      TFile *                  Eff_file;
      TFile *                  Wimp_file;
      TF1   *                  Smeared_spectrum_fct;
      TF1   *                  Quantified_spectrum_fct;
      TF1   *                  Effcorr_Quantified_spectrum_fct;
      TF1   *                  Spectrum_function;
      TF1   *                  Efficiency_function;
      TF1   *                  Effcorr_Wimp_fct;
      TF1   *                  Smeared_Wimp_fct;
      TGraph*                  Efficiency_graph;
      TGraph*                  Spectrum_graph;
      TGraph*                  Wimp_graph;
      std::string              Detector;
	  std::string              Output_name;
	  std::string              List_mass;	
	  std::string              Eff_NAME;	 
      double                   Tension;
      double                   Percent_running;
      double                   Fano;
      double                   Verbosity;
      double                   Mass=1000.;
      std::string              Mass_str;
      bool                     Use_probCDMS = true ;
      double                   Resolution;  
      std::string              Eicut;
};
#endif
