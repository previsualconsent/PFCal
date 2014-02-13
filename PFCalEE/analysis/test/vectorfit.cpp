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
#include "TFitResult.h"

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

  if (argc < 4) {
    std::cout << " Usage: " 
      << argv[0] << " <nEvts to process (0=all)>"
      << " <input_dir>" 
      << " <angle>" 
      << " <optional: debug (default=0)>"
      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  const std::string inDir = argv[2];
  bool debug = false;
  std::string angle(argv[3]);
  if (argc >4) debug = atoi(argv[4]);

  TFile *outputFile = TFile::Open("vectorfit.root","RECREATE");

  const unsigned nLayers = N_LAYERS;

  unsigned genEn[]={5,10,25,50,75,100,150};
  //unsigned genEn[]={10};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  TH2F *p_angle[nGenEn];


  double Emax[nGenEn];

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    Emax[iE] = 0;


    std::ostringstream lName;

    lName << "p_angle_" << genEn[iE];
    p_angle[iE] = new TH2F(lName.str().c_str(),";x angle(rad);y angle(rad)",100,-.5,.5,100,-.5,.5);
    lName.str("");

    std::ostringstream input;
    input << inDir << genEn[iE] << "_" << angle << "/PFcal.root";
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


    lTree->GetEntry(0);
    float curEvent(event);

    TH2F *p_xy[nEvts][nLayers];
    float sumw2[nEvts][nLayers];
    for (unsigned i(0); i<nEvts; ++i){
        for (unsigned iL(0); iL<nLayers; ++iL){
            lName << "p_xy_" << i << "_" << iL;
            p_xy[i][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100);
            lName.str("");

            sumw2[i][iL] = 0;
        }
    }

    for (unsigned ievt(0); ievt < nEvts*nLayers; ++ievt){//loop on entries

        lTree->GetEntry(ievt);


      if (ievt%(nLayers*100) == 0) std::cout << "... Processing event: " << ievt/nLayers << std::endl;


      unsigned layer = volNb;
      if (debug) std::cout << "... Processing layer " << layer << " with " << (*hitvec).size() << " simhits." << std::endl;

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
        if(weightedE) 
        {
            p_xy[int(event)][layer]->Fill(posx,posy,weightedE);
            sumw2[int(event)][layer]+=TMath::Power(weightedE,2);
        }
      }//loop on hits

    }//loop on entries

    TGraph * xgraph;
    TGraph * ygraph;
    for (unsigned i(0); i<nEvts; ++i)
    {
        if (i%100 == 0) std::cout << "... Fitting event: " << i << std::endl;
        int nPoints = -1;
        int nFitPoints = 5;
        float xavg[nLayers];
        float yavg[nLayers];
        float xerr[nLayers];
        float yerr[nLayers];
        float layers[nLayers];
        for (int iL(0); iL<nLayers; ++iL)
        {
            if(p_xy[i][iL]->GetEntries() > 0)
            {
                nPoints++;
                //std::cout << "supposedly valid: " << layers[iL] << ":" << xavg[iL] << ":" << xerr[iL] << std::endl;;
            }
            else continue;

            layers[nPoints]=float(iL);
            xavg[nPoints]=p_xy[i][iL]->GetMean(1);
            yavg[nPoints]=p_xy[i][iL]->GetMean(2);

            //std::cout << sumw2[i][iL]  << ":"
            //          << p_xy[i][iL]->GetSumOfWeights() << ":"
            //          << p_xy[i][iL]->Integral() << ":"
            //          << p_xy[i][iL]->GetEntries() << ":"
            //          << std::endl;
            xerr[nPoints]=200.0f/(80.0f*TMath::Sqrt(12))*TMath::Sqrt(sumw2[i][iL]/TMath::Power(p_xy[i][iL]->Integral(),2));
            yerr[nPoints]=xerr[nPoints];

            if(nPoints <= nFitPoints && xerr[nPoints]<=0.0 )
            {
                std::cout << "where are the zeros: " << layers[nPoints] << ":" << xavg[nPoints] << ":" << xerr[nPoints] << std::endl;
            }


            //xerr[nPoints]=200.0f/(80.0f*TMath::Sqrt(12))*TMath::Sqrt(;
            //yerr[nPoints]=200.0f/(80.0f*TMath::Sqrt(12));



        }
        float * layererr = 0;

        if(nPoints >1)
        {
            TCanvas * c = new TCanvas("asalkfj","asfka");

            xgraph = new TGraphErrors(nPoints,layers,xavg,layererr,xerr);
            xgraph->Draw("AP");
            xgraph->Fit("pol1","Q","",0,nFitPoints);

            TF1 * fit = xgraph->GetFunction("pol1");
            float xangle =  fit->Derivative(0)/6.80f;

            std::ostringstream title;
            title << "angle: " << xangle;
            xgraph->SetTitle(title.str().c_str());

            if(event<10)
            {
                std::ostringstream saveName;
                saveName << "PLOTS/Fit_" << i << "_" << genEn[iE] << "GeV" << ".png";
                c->SaveAs(saveName.str().c_str());
            }
            delete c;



            ygraph = new TGraphErrors(nPoints,layers,yavg,layererr,yerr);
            ygraph->Fit("pol1","Q","",0,nFitPoints);
            fit = ygraph->GetFunction("pol1");
            float yangle =  fit->Derivative(0)/6.80f;

            p_angle[iE]->Fill(xangle,yangle);
        }
        if(i == 9) std::cout << "alskaf" <<  p_xy[i][0]->GetEntries() << std::endl;
    }

  }//loop on energies

  TCanvas *myc = new TCanvas("myc","myc",1);

  //gStyle->SetOptStat(0);

  std::ostringstream saveName;
  for (unsigned iE(0); iE<nGenEn; ++iE){
      myc->cd();
      std::ostringstream title;
      title << "Energy  " << genEn[iE] ;
      p_angle[iE]->SetTitle(title.str().c_str());
      p_angle[iE]->Draw("colz");
      myc->Update();
      saveName.str("");
      saveName << "PLOTS/angle_hist_" << genEn[iE] << "GeV";
      myc->Print((saveName.str()+".png").c_str());
      myc->Print((saveName.str()+".pdf").c_str());
      outputFile->cd();
      p_angle[iE]->Write();
  }

  outputFile->Write();
  return 0;


}//main
