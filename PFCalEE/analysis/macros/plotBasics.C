#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"

int plotBasics(){//main


  TString plotDir = "../PLOTS/version_20/scenario_0/";

  //const double Emip = 0.0559;//in MeV
  //const unsigned nLayers = 30;

  TFile *inputFile = TFile::Open(plotDir+"/CalibHistos.root");
  if (!inputFile) {
    std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TCanvas *myc = new TCanvas("myc","myc",1);//1500,1000);

  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //parameters to change
  TString genEn = "500";
  //list histos
  TH1F *p_nSimHits = (TH1F*)gDirectory->Get("p_nSimHits_"+genEn);
  TH1F *p_nRecHits = (TH1F*)gDirectory->Get("p_nRecHits_"+genEn);
  
  std::ostringstream saveName;
  saveName.str("");
  saveName << plotDir << "/nHits_" << genEn;
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  

  gStyle->SetOptStat(1110);
  myc->Divide(2,1);
  myc->cd(1);
  p_nSimHits->Rebin(200);
  float xMin = p_nSimHits->GetMean()-5*p_nSimHits->GetRMS();
  float xMax = p_nSimHits->GetMean()+5*p_nSimHits->GetRMS();
  float yMax = p_nSimHits->GetMaximum();
  p_nSimHits->GetXaxis()->SetRangeUser(xMin,xMax);
  p_nSimHits->Draw("PE");

  TLatex lat;
  lat.DrawLatex(xMin,yMax*1.2,"e-, "+genEn+" GeV, 2.5#times2.5 mm^{2} cells");

  myc->cd(2);
  p_nRecHits->Rebin(20);
  xMin = p_nRecHits->GetMean()-5*p_nRecHits->GetRMS();
  xMax = p_nRecHits->GetMean()+5*p_nRecHits->GetRMS();
  yMax = p_nRecHits->GetMaximum();
  p_nRecHits->GetXaxis()->SetRangeUser(xMin,xMax);
  p_nRecHits->Draw("PE");
  lat.DrawLatex(xMin,yMax*1.2,"e-, "+genEn+" GeV, 1#times1 cm^{2} cells");

  myc->Update();
  myc->Print("Basics/nHits_"+genEn+".png");
  myc->Print("Basics/nHits_"+genEn+".pdf");




  return 0;

}//main
