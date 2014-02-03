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

int plotTime(){//main


  TString plotDir = "../PLOTS/version_20/";

  TFile *inputFile = TFile::Open(plotDir+"CalibHistos.root");
  if (!inputFile) {
    std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  const unsigned nLayers = 30;

  //unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
  unsigned genEn[]={100};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
  TH1F *p_time[nGenEn][nLayers];
  TH1F *p_timeTot[nGenEn];

  TCanvas *myc = new TCanvas("myc","myc",1);
    
  gStyle->SetOptStat(0);
  
  for (unsigned iE(0); iE<nGenEn; ++iE){
    
    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    myc->cd();    
    std::ostringstream lName;
    double maxY = 0;
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName.str("");
      lName << "p_time_" << genEn[iE] << "_" << iL;
      p_time[iE][iL] = (TH1F*)gDirectory->Get(lName.str().c_str());
      if (!p_time[iE][iL]) {
	std::cout << " -- ERROR, pointer for histogram Etot is null for layer: " << iL << ". Exiting..." << std::endl;
	return 1;
      }
      if (p_time[iE][iL]->GetMaximum()>maxY) maxY = p_time[iE][iL]->GetMaximum();
      //if (iL==0) p_timeTot[iE] = (TH1F*)p_time[iE][iL]->Clone();
      //else p_timeTot[iE]->Add(p_time[iE][iL]);
      std::cout << " Layer " << iL << " mean = " <<p_time[iE][iL]->GetMean() << " RMS = " << p_time[iE][iL]->GetRMS() << std::endl;
    }
    for (unsigned iL(0); iL<nLayers; ++iL){
      p_time[iE][iL]->SetLineColor(iL%9+1);
      p_time[iE][iL]->SetMaximum(maxY);
      if (iL==0){
	p_time[iE][iL]->Draw();
      }
      else {
	//p_time[iE][iL]->Add(p_time[iE][iL-1]);
	p_time[iE][iL]->Draw("same");
      }
    }
    
    std::ostringstream saveName;
    saveName.str("");
    saveName << plotDir << "/SimG4Time_" << genEn[iE] << "GeV";
    myc->Update();
    myc->Print((saveName.str()+".png").c_str());
    myc->Print((saveName.str()+".pdf").c_str());
    
  }//loop on energies

  return 0;

  
}//main
