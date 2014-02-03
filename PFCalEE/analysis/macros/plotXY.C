#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

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


int plotXY(){//main  

  TString plotDir = "../PLOTS/version_0/";

  TFile *inputFile = TFile::Open(plotDir+"CalibHistos.root");
  if (!inputFile) {
    std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  const unsigned nLayers = 30;

  unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
  //unsigned genEn[]={10};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
  TH2F *p_xy[nGenEn][nLayers];
  double Emax[nGenEn];

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    Emax[iE] = 0;
    
    for (unsigned iL(0); iL<nLayers; ++iL){
      std::ostringstream lName;
      lName << "p_xy_" << genEn[iE] << "_" << iL;
      p_xy[iE][iL] = (TH2F*)gDirectory->Get(lName.str().c_str());
      if (!p_xy[iE][iL]) {
	std::cout << " -- ERROR, pointer for histogram is null for layer: " << iL << ". Exiting..." << std::endl;
	return 1;
      }
      double Etot = p_xy[iE][iL]->GetMaximum();
      if (Etot > Emax[iE]) Emax[iE] = Etot;
    }

    std::cout << " -- max energy " << Emax[iE] << std::endl;
    
    TCanvas *mycAll = new TCanvas("mycAll","mycAll",1500,1000);
    TCanvas *myc = new TCanvas("myc","myc",1);
    mycAll->Divide(6,5);
    
    gStyle->SetOptStat(0);
    
    std::ostringstream saveName;
    for (unsigned iL(0); iL<nLayers; ++iL){
      mycAll->cd(iL+1);
      p_xy[iE][iL]->SetMaximum(Emax[iE]);
      p_xy[iE][iL]->Draw("colz");
      myc->cd();
      p_xy[iE][iL]->Draw("colz");
      myc->Update();
      saveName.str("");
      saveName << plotDir << "/xySimHits_layer" << iL << "_" << genEn[iE] << "GeV";
      myc->Print((saveName.str()+".png").c_str());
      myc->Print((saveName.str()+".pdf").c_str());
    }
    saveName.str("");
    saveName << plotDir << "/xySimHits_" << genEn[iE] << "GeV";
    mycAll->Update();
    mycAll->Print((saveName.str()+".png").c_str());
    mycAll->Print((saveName.str()+".pdf").c_str());
    
  }//loop on energies

  return 0;


}//main
