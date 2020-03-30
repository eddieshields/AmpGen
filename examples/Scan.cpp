#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"
#include "AmpGen/ErrorPropagator.h"
#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>

using namespace AmpGen;

template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList& data, EventList& mc, MinuitParameterSet& MPS, int nBins );

std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps);
std::map<std::string, double > getPulls(std::map<std::string, std::vector<double> > fits, std::map<std::string, std::vector<double> > inits);
void writePulls(std::string fileName, std::map<std::string, double> pulls);



int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log", "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root", "Name of the output plot file");
  
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();

  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
  


  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
  size_t      nEvents  = NamedParameter<size_t>     ("nEvents"   , 10000       , "Number of events to fill in") ;
  int nBins = NamedParameter<int> ("nBins", 100, "number of bins for projection");




   

  if( dataFile == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("LogFile: " << logFile << "; Plots: " << plotFile );
  
#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  /* A MinuitParameterSet is (unsurprisingly) a set of fit parameters, and can be loaded from 
     the parsed options. For historical reasons, this is referred to as loading it from a "Stream" */
  MinuitParameterSet MPS;
  MPS.loadFromStream();
  

  std::map<std::string, std::vector<double> > inits = getParams(MPS);
  /* An EventType specifies the initial and final state particles as a vector that will be described by the fit. 
     It is typically loaded from the interface parameter EventType. */
  EventType evtType(pNames);
  /* A CoherentSum is the typical amplitude to be used, that is some sum over quasi two-body contributions 
     weighted by an appropriate complex amplitude. The CoherentSum is generated from the couplings described 
     by a set of parameters (in a MinuitParameterSet), and an EventType, which matches these parameters 
     to a given final state and a set of data. A common set of rules can be matched to multiple final states, 
     i.e. to facilitate the analysis of coupled channels. 
     The CoherentSum is only appropriate for decays involving only (pseudo)scalars in the inital / final state, 
     otherwise the sum must also be over initial / final spin states. In this case, as PolarisedSum should be used. 
     See FitterWithPolarisation for an example of this use case.    
  */
  CoherentSum sig(evtType, MPS);
  
  /* Events are read in from ROOT files. If only the filename and the event type are specified, 
     the file is assumed to be in the specific format that is defined by the event type, 
     unless the branches to load are specified in the user options */
  EventList events(dataFile, evtType, Branches(bNames), GetGenPdf(false) );
  
  /* Generate events to normalise the PDF with. This can also be loaded from a file, 
     which will be the case when efficiency variations are included. Default number of normalisation events 
     is 5 million. */
  EventList eventsMC = intFile == "" ? Generator<>(evtType, &rndm).generate(2e6) : EventList(intFile, evtType, GetGenPdf(true));
  
  sig.setMC( eventsMC );
  sig.prepare();

  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
  

   


  /* Do the fit and return the fit results, which can be written to the log and contains the 
     covariance matrix, fit parameters, and other observables such as fit fractions */
  FitResult* fr = doFit(make_likelihood(events, sig), events, eventsMC, MPS , nBins);
  /* Calculate the `fit fractions` using the signal model and the error propagator (i.e. 
     fit results + covariance matrix) of the fit result, and write them to a file. 
   */
  auto fitFractions = sig.fitFractions( fr->getErrorPropagator() ); 
  
  fr->addFractions( fitFractions );
  fr->writeToFile( logFile );
       std::map<std::string, std::vector<double> > fits = getParams(MPS);
    std::map<std::string, double> pulls = getPulls(fits, inits);
      for(std::map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
         INFO("Pull = "<<it->first<<" "<<it->second);
       }
    writePulls(logFile, pulls);


  output->cd();
  
  /* Write out the data plots. This also shows the first example of the named arguments 
     to functions, emulating python's behaviour in this area */

  auto plots = events.makeDefaultProjections(Prefix("Data"), Bins(nBins));
  for ( auto& plot : plots ) plot->Write();

  output->Close();
}

template <typename likelihoodType>
FitResult* doFit( likelihoodType&& likelihood, EventList& data, EventList& mc, MinuitParameterSet& MPS, int nBins )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();
  /* Minimiser is a general interface to Minuit1/Minuit2, 
     that is constructed from an object that defines an operator() that returns a double 
     (i.e. the likielihood, and a set of MinuitParameters. */
  Minimiser mini( likelihood, &MPS );


      std::ofstream scanfile;
      auto pName = "D0{K*(892)bar-{K0S0,pi-},pi+}_Re";
        scanfile.open("Scan.txt", std::ios_base::app);
        auto param = MPS[pName];
        double minimum=param->minInit();
        double maximum=param->maxInit();
        double val = minimum;
        double stepSize = param->stepInit();

        MPS[pName]->setCurrentFitVal(val);   


        while (val < maximum){
    

     MPS[pName]->setCurrentFitVal(val);    

        
       
    INFO(pName<<" = "<<val<<", Function = "<<mini.FCN() );
        scanfile<<val<<"\t"<<mini.FCN()<<"\n";      
        val += stepSize;
        
        }
       
    scanfile.close();


  auto covar = mini.covMatrix();
  INFO("Printing Covariant matrix");
  covar.Print();
  auto covarFull = mini.covMatrixFull();
  //INFO("Printing full Covariant matrix");
  //covarFull.Print();
  mini.gradientTest();
  mini.doFit();

  FitResult* fr = new FitResult(mini);
  
  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );

  /* Make the plots for the different components in the PDF, i.e. the signal and backgrounds. 
     The structure assumed the PDF is some SumPDF<eventListType, pdfType1, pdfType2,... >. */
  unsigned int counter = 1;
  for_each(likelihood.pdfs(), [&](auto& pdf){
    auto pfx = Prefix("Model_cat"+std::to_string(counter));
    auto mc_plot3 = mc.makeDefaultProjections(WeightFunction(pdf), Bins(nBins), pfx);
    for( auto& plot : mc_plot3 )
    {
      plot->Scale( ( data.integral() * pdf.getWeight() ) / plot->Integral() );
      plot->Write();
    }
    counter++;
  });

  /* Estimate the chi2 using an adaptive / decision tree based binning, 
     down to a minimum bin population of 15, and add it to the output. */
  Chi2Estimator chi2( data, mc, likelihood, 15 );
  chi2.writeBinningToFile("chi2_binning.txt");
  fr->addChi2( chi2.chi2(), chi2.nBins() );
  
  fr->print();
  return fr;
}
std::map<std::string, std::vector<double> > getParams(MinuitParameterSet & mps){
  std::map<std::string, std::vector<double> > out;
  for (auto& param : mps){
    std::vector<double> vect = {param->mean(), param->err()};
    out[param->name()] = vect;
  }
  return out;
}




std::map<std::string, double > getPulls(std::map<std::string, std::vector<double> > fits, std::map<std::string, std::vector<double> > inits)  {
  std::string output = "";
 std::map<std::string, double> out;
 for(std::map<std::string, std::vector<double> >::iterator init = inits.begin(); init != inits.end(); ++init) {
  for(std::map<std::string, std::vector<double> >::iterator fit = fits.begin(); fit != fits.end(); ++fit) {
    std::string initName = init->first;
    std::vector<double> initParams = init->second;
    std::string fitName = fit->first;
    std::vector<double> fitParams = fit->second;
    if (initName == fitName){
        double pull = fitParams[0] - initParams[0];
        if (fitParams[1] != 0){
          pull /= fitParams[1];
        }
        out[fitName] = pull;
    }
  }
 }
  return out;

}


void writePulls(std::string fileName, std::map<std::string, double> pulls){
  std::ofstream outfile;
  outfile.open(fileName, std::ios_base::app);
      for(std::map<std::string, double >::iterator it = pulls.begin(); it != pulls.end(); ++it) {
        outfile<<"Pull "<<it->first<<" "<<it->second<<"\n";
       }

}


