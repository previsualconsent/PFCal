#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

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

  TFile *outputFile = TFile::Open("vectorfit.root","RECREATE");

  const unsigned nLayers = N_LAYERS;

  //unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
  unsigned genEn[]={10};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  TH2F *p_xy[nGenEn][nLayers];
  TH3F *p_xyz[nGenEn];
  TGraph * xgraph[nGenEn];
  TGraph * ygraph[nGenEn];

  double Emax[nGenEn];

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    Emax[iE] = 0;


    std::ostringstream lName;
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName << "p_xy_" << genEn[iE] << "_" << iL;
      p_xy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100);
      lName.str("");
    }
    lName << "p_xyz_" << genEn[iE];
    p_xyz[iE] = new TH3F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100,nLayers,0,nLayers);
    lName.str("");

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

    const unsigned nEvts = (pNevts > lTree->GetEntries()/nLayers || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/nLayers) : pNevts;

    std::cout << "- Processing = " << nEvts  << " events out of " << lTree->GetEntries()/nLayers << std::endl;

    unsigned nPoints = 0;
    float xavg[nLayers];
    float yavg[nLayers];
    float xerr[nLayers];
    float yerr[nLayers];
    float layers[nLayers];
    for (unsigned ievt(0); ievt<nEvts*nLayers; ++ievt){//loop on entries
      if (ievt%(nLayers*100) == 0) std::cout << "... Processing event: " << ievt/nLayers << std::endl;

      lTree->GetEntry(ievt);

      unsigned layer = volNb;
      if (debug) std::cout << "... Processing layer " << layer << " with " << (*hitvec).size() << " simhits." << std::endl;

      bool hits = false;
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
	double weightedE = lHit.energy()*getWeight(layer,inDir);
        if (weightedE > Emax[iE]) Emax[iE] = weightedE ;
        p_xy[iE][layer]->Fill(posx,posy,weightedE);
        p_xyz[iE]->Fill(posx,posy,layer,weightedE);
        if(weightedE) hits = true;
      }//loop on hits
      if(hits)
      {
        layers[nPoints]=float(layer);
        xavg[nPoints]=p_xy[iE][layer]->GetMean(1);
        yavg[nPoints]=p_xy[iE][layer]->GetMean(2);
        //xerr[nPoints]=p_xy[iE][layer]->GetRMS(1);
        //yerr[nPoints]=p_xy[iE][layer]->GetRMS(2);
        float clayer = 0.25f;
        float cint = 0.0f;
        float crms = 1.0f;
        xerr[nPoints]=TMath::Sqrt(TMath::Power(layer*clayer,2)+
            TMath::Power(p_xy[iE][layer]->GetRMS(1)*crms,2) +
            TMath::Power(p_xy[iE][layer]->Integral()*cint,2));
        yerr[nPoints]=TMath::Sqrt(TMath::Power(layer*clayer,2)+
            TMath::Power(p_xy[iE][layer]->GetRMS(2)*crms,2) +
            TMath::Power(p_xy[iE][layer]->Integral()*cint,2));

        std::cout << "nPoints: " << nPoints;
        std::cout << " xavg[nPoints]: " <<  xavg[nPoints];
        std::cout << " yavg[nPoints]: " << yavg[nPoints];
        std::cout << " xerr[nPoints]: " << xerr[nPoints];
        std::cout << " yerr[nPoints]: " << yerr[nPoints];
        std::cout << std::endl;

        nPoints++;
      }

    }//loop on entries

    float * layererr = 0;
    xgraph[iE] = new TGraphErrors(nPoints,layers,xavg,layererr,xerr);
    ygraph[iE] = new TGraphErrors(nPoints,layers,yavg,layererr,yerr);

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

      std::ostringstream title;
      std::string pad="";
      if(iL < 10) pad = "0";
      title << "Layer " << pad << iL;
      p_xy[iE][iL]->SetTitle(title.str().c_str());
      p_xy[iE][iL]->Draw("colz");
      myc->Update();
      saveName.str("");
      saveName << "PLOTS/xySimHits_layer" << pad << iL << "_" << genEn[iE] << "GeV";
      myc->Print((saveName.str()+".png").c_str());
      myc->Print((saveName.str()+".pdf").c_str());
      outputFile->cd();
      p_xy[iE][iL]->Write();
    }
    p_xyz[iE]->Write();
    saveName.str("");
    saveName << "PLOTS/xySimHits_" << genEn[iE] << "GeV";
    mycAll->Update();
    mycAll->Print((saveName.str()+".png").c_str());
    mycAll->Print((saveName.str()+".pdf").c_str());

    TCanvas * c = new TCanvas("xgraph","xgraph");
    xgraph[iE]->Draw("AP");
    xgraph[iE]->Fit("pol2");

    TF1 * fit = xgraph[iE]->GetFunction("pol2");
    std::ostringstream title;
    title << "x slope: " << fit->Derivative(0);
    xgraph[iE]->SetTitle(title.str().c_str());

    saveName.str("");
    saveName << "PLOTS/xgraph_" << genEn[iE] << "GeV";
    c->Print((saveName.str()+".png").c_str());
    delete c;

    c = new TCanvas("ygraph","ygraph");
    ygraph[iE]->SetTitle("Y graph");
    ygraph[iE]->Draw("AP");
    ygraph[iE]->Fit("pol2");

    fit = ygraph[iE]->GetFunction("pol2");
    title.str("");
    title << "y slope: " << fit->Derivative(0);
    ygraph[iE]->SetTitle(title.str().c_str());

    saveName.str("");
    saveName << "PLOTS/ygraph_" << genEn[iE] << "GeV";
    c->Print((saveName.str()+".png").c_str());
    delete c;
  }

  outputFile->Write();
  return 0;


}//main
