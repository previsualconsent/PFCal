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

#include "HGCSSSimHit.hh"
#include "HGCSSParameters.hh"

double getWeight(unsigned layer,std::string aVersion){
  if (layer<10) return 1;
  if (aVersion.find("20")!= aVersion.npos || aVersion.find("21") != aVersion.npos){
    if (layer < 20) return 0.8/0.5;
    else return 1.2/0.5;
  }
  else {
    if (layer < 20) return 0.8/0.4;
    else return 1.2/0.4;
  }
}


int main(int argc, char** argv){//main  

  if (argc < 3) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
              << " <input_dir>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  const std::string inDir = argv[2];
  bool debug = false;
  if (argc >3) debug = atoi(argv[3]);

  TFile *outputFile = TFile::Open("DigiHistos.root","RECREATE");

  const unsigned nLayers = N_LAYERS;

  //unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
  unsigned genEn[]={10};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
  TH2F *p_xy[nGenEn][nLayers];
  TH1F *p_Etot[nGenEn][nLayers];

  double Emax[nGenEn];

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    Emax[iE] = 0;
    
    double Etot[nLayers];

    for (unsigned iL(0); iL<nLayers; ++iL){
      std::ostringstream lName;
      lName << "p_xy_" << genEn[iE] << "_" << iL;
      p_xy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100);
      lName.str("");
      lName << "p_Etot_" << genEn[iE] << "_" << iL;
      p_Etot[iE][iL] = new TH1F(lName.str().c_str(),";Etot (GeV)",1000,0,100);
      lName.str("");
      Etot[iL] = 0;
    }
    
    std::ostringstream input;
    input << inDir << genEn[iE] << "_0.250" << "/PFcal.root";
    TFile *inputFile = TFile::Open(input.str().c_str());

    if (!inputFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
      return 1;
    }

    TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
    if (!lTree){
      std::cout << " -- Error, ntuple CaloStack  cannot be opened. Exiting..." << std::endl;
      return 1;
    }
    
    float event;
    float volNb;
    float volX0;
    float volX0trans;
    float den;
    float denWeight;
    float denAbs;
    float denTotal;
    float gFrac;
    float eFrac;
    float muFrac;
    float hadFrac;
    float avgTime;
    float nSiHits;
    std::vector<HGCSSSimHit> * hitvec = 0;

    lTree->SetBranchAddress("event",&event);
    lTree->SetBranchAddress("volNb",&volNb);
    lTree->SetBranchAddress("volX0",&volX0);
    lTree->SetBranchAddress("volX0trans",&volX0trans);
    lTree->SetBranchAddress("den",&den);
    lTree->SetBranchAddress("denWeight",&denWeight);
    lTree->SetBranchAddress("denAbs",&denAbs);
    lTree->SetBranchAddress("denTotal",&denTotal);
    lTree->SetBranchAddress("gFrac",&gFrac);
    lTree->SetBranchAddress("eFrac",&eFrac);
    lTree->SetBranchAddress("muFrac",&muFrac);
    lTree->SetBranchAddress("hadFrac",&hadFrac);
    lTree->SetBranchAddress("avgTime",&avgTime);
    lTree->SetBranchAddress("nSiHits",&nSiHits);
    lTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
    
    const unsigned nEvts = (pNevts > lTree->GetEntries()/30. || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/30.) : pNevts;

    std::cout << "- Processing = " << nEvts  << " events out of " << lTree->GetEntries()/30. << std::endl;

    for (unsigned ievt(0); ievt<nEvts*30; ++ievt){//loop on entries
      if (ievt%3000 == 0) std::cout << "... Processing event: " << ievt/30 << std::endl;
 
      lTree->GetEntry(ievt);

      unsigned layer = volNb;
      if (debug) std::cout << "... Processing layer " << layer << " with " << (*hitvec).size() << " simhits." << std::endl;
      Etot[layer] = 0;
      
      for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*hitvec)[iH];
	lHit.layer(layer);
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	if (debug) {
	  std::cout << " --  Hit " << iH << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}
	double weightedE = lHit.energy()*getWeight(layer,pVersion);
	if (weightedE > Emax[iE]) Emax[iE] = weightedE ;
	p_xy[iE][layer]->Fill(posx,posy,weightedE);
	Etot[layer] += weightedE;
      }//loop on hits
      p_Etot[iE][layer]->Fill(Etot[layer]);
    }//loop on entries

    std::cout << " -- max energy " << Emax[iE] << std::endl;
    
  }//loop on energies

  TCanvas *mycAll = new TCanvas("mycAll","mycAll",1);
  TCanvas *myc = new TCanvas("myc","myc",1);
  mycAll->Divide(5,6);
  
  gStyle->SetOptStat(0);
  
  std::ostringstream saveName;
  for (unsigned iE(0); iE<nGenEn; ++iE){
    for (unsigned iL(0); iL<nLayers; ++iL){
      mycAll->cd(iL+1);
      p_xy[iE][iL]->SetMaximum(Emax[iE]);
      p_xy[iE][iL]->Draw("colz");
      myc->cd();
      p_xy[iE][iL]->Draw("colz");
      myc->Update();
      saveName.str("");
      saveName << "PLOTS/xySimHits_layer" << iL << "_" << genEn[iE] << "GeV";
      myc->Print((saveName.str()+".png").c_str());
      myc->Print((saveName.str()+".pdf").c_str());
      outputFile->cd();
      p_xy[iE][iL]->Write();
      p_Etot[iE][iL]->Write();
    }
    saveName.str("");
    saveName << "PLOTS/xySimHits_" << genEn[iE] << "GeV";
    mycAll->Update();
    mycAll->Print((saveName.str()+".png").c_str());
    mycAll->Print((saveName.str()+".pdf").c_str());
  }

  outputFile->Write();
  return 0;


}//main
