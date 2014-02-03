#ifndef HGCSSCaloProperties_h
#define HGCSSCaloProperties_h

#include <map>
#include <string>

#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpline.h"

#include "Math/WrappedTF1.h"
#include "Math/BrentMinimizer1D.h"


#include "HGCSSParameters.hh"

typedef std::pair<Float_t,Float_t> Measurement_t;

/**
   @class GraphToTF1
   @short a wrapper to use a graph as a TF1 through the TSplines
 */
class GraphToTF1{
public:
  GraphToTF1(TString name, TGraph *g);
  ~GraphToTF1(){};

  double operator()(double *x,double *p);
  TSpline *sp_;
};

/**
   @class ShowerProfile
   @short store histograms and graphs regarding a given shower profile
 */
class ShowerProfile
{
public:
  ShowerProfile();
  ~ShowerProfile(){};

  void writeTo(TDirectory *dir);
  bool buildShowerProfile(Float_t eElec, TString version);

  TH1F         *h_rawEn,          *h_en,    *h_enFit, *h_showerMax;
  TH2F         *h_enVsOverburden, *h_enVsDistToShowerMax;
  TGraphErrors *gr_raw,            *gr_centered;
  TGraphErrors *gr_unc,            *gr_relUnc;
};

/**
   @class ShowerProfile
   @short store histograms and graphs regarding a given shower profile
 */
class DigiShowerProfile
{
public:
  DigiShowerProfile();
  ~DigiShowerProfile(){};

  void writeTo(TDirectory *dir);
  bool buildShowerProfile(Float_t eElec, TString version);

  TH1F         *h_energy;
};

/**
   @class CaloProperties
   @short stores all shower properties built and the fits to the calibration and performance
 */

class CaloProperties
{
public:
  CaloProperties(TString tag);
  void setEnergiesToScan(std::vector<Float_t> &enList) { genEn_=enList; }
  ~CaloProperties(){};

  void writeTo(TDirectoryFile *dir);

  void characterizeCalo();

  TString tag_;
  TGraph *gr_showerMax, *gr_centeredShowerMax;
  std::vector<TGraphErrors *> calibCurve_,resCurve_;
  std::vector<Measurement_t> stochTerms_, constTerms_;
  std::map<Float_t,ShowerProfile> showerProfiles_;
  std::vector<Float_t> genEn_;
};

void drawHeader();
void setStyle();

#endif
