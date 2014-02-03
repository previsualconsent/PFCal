#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
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
	      << " <full path to input file>"
	      << " <file name (PFcal.root, or DigiPFcal.root)>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string fileName = argv[3];//"PFcal.root";
  unsigned debug = 0;
  if (argc >4) debug = atoi(argv[4]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Input file name: " << fileName << std::endl
	    << " -- Processing " << pNevts << " events." << std::endl;

  TString plotDir = "PLOTS/";
  std::vector<std::string> pathVec;
  boost::split( pathVec, filePath, boost::is_any_of("/"));
  std::string pVersion;
  std::string pScenario;
  std::string tmpV, tmpS;
  if (pathVec[pathVec.size()-1]=="")  {
    tmpS = pathVec[pathVec.size()-2];
    tmpV = pathVec[pathVec.size()-3];
  }
  else {
    tmpS = pathVec[pathVec.size()-1];
    tmpV = pathVec[pathVec.size()-2];
  }
  if (tmpS.find("version") != tmpS.npos){
    pVersion = tmpS;
    pScenario = "";
  }
  else {
    pVersion = tmpV;
    pScenario = tmpS;
  }

  plotDir += pVersion;
  plotDir += "/";
  plotDir += pScenario;
  plotDir += "/";

  std::cout << " -- Output file directory is : " << plotDir << std::endl;

  TFile *outputFile = TFile::Open(plotDir+"/CalibHistos.root","RECREATE");
  if (!outputFile) {
    std::cout << " -- Error, output file " << plotDir << "/CalibHistos.root cannot be opened. Please create output directory : " << plotDir << ". Exiting..." << std::endl;
    return 1;
  }


  const unsigned nLayers = N_LAYERS;

  unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
  //unsigned genEn[]={100};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  
  TH2F *p_xy[nGenEn][nLayers];
  TH1F *p_Etot[nGenEn][nLayers];
  TH1F *p_Efrac[nGenEn][nLayers];
  TH1F *p_time[nGenEn][nLayers];
  TH1F *p_Etotal[nGenEn];
  TH1F *p_Ereco[nGenEn];
  TH1F *p_nSimHits[nGenEn];
  TH1F *p_nRecHits[nGenEn];

  bool isG4Tree = true;

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "- Processing energy : " << genEn[iE] << std::endl;

    std::ostringstream input;
    input << filePath << "/e-/e_" << genEn[iE] << "/" << fileName ;
    TFile *inputFile = TFile::Open(input.str().c_str());

    if (!inputFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
      return 1;
    }

    TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
    if (!lTree){
      std::cout << " -- Error, tree HGCSSTree cannot be opened. Trying RecoTree instead." << std::endl;
      isG4Tree = false;
      lTree = (TTree*)inputFile->Get("RecoTree");
      if (!lTree){
	std::cout << " -- Error, tree RecoTree cannot be opened either. Exiting..." << std::endl;
	return 1;
      }
    }

    //Initialise histos
    //necessary to have overflows ?
    gStyle->SetOptStat(1111111);
    double Etot[nLayers];
    for (unsigned iL(0);iL<nLayers;++iL){
      Etot[iL] = 0;
    }
    double Etotal = 0;
    double Ereco = 0;
    std::ostringstream lName;
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName.str("");
      lName << "p_xy_" << genEn[iE] << "_" << iL;
      p_xy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100);
      lName.str("");
      lName << "p_Etot_" << genEn[iE] << "_" << iL;
      p_Etot[iE][iL] = new TH1F(lName.str().c_str(),";G4 Etot (MeV)",800,0,800);
      lName.str("");
      Etot[iL] = 0;
      lName.str("");
      lName << "p_Efrac_" << genEn[iE] << "_" << iL;
      p_Efrac[iE][iL] = new TH1F(lName.str().c_str(),";integrated E_{layer}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_time_" << genEn[iE] << "_" << iL;
      p_time[iE][iL] = new TH1F(lName.str().c_str(),";G4 time (ns)",200,0,2);
    }
    lName.str("");
    lName << "p_Etotal_" << genEn[iE];
    p_Etotal[iE] = new TH1F(lName.str().c_str(),";G4 Etotal (MeV)",7000,0,7000);
    p_Etotal[iE]->StatOverflows();
    lName.str("");
    lName << "p_nSimHits_" << genEn[iE];
    p_nSimHits[iE] = new TH1F(lName.str().c_str(),"; nSimHits",50000,0,50000);
    if (!isG4Tree){
      lName.str("");
      lName << "p_Ereco_" << genEn[iE];
      p_Ereco[iE] = new TH1F(lName.str().c_str(),";Reco Etotal (MIPs)",10000,0,110000);
      p_Ereco[iE]->StatOverflows();
      lName.str("");
      lName << "p_nRecHits_" << genEn[iE];
      p_nRecHits[iE] = new TH1F(lName.str().c_str(),"; nRecHits",6000,0,6000);
    }

    unsigned event = 0;
    float volNb = 0;
    float volX0 = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    if (isG4Tree){
      //lTree->SetBranchAddress("event",&event);
      lTree->SetBranchAddress("volNb",&volNb);
      lTree->SetBranchAddress("volX0",&volX0);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    }
    else {
      lTree->SetBranchAddress("event",&event);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    }

    const unsigned nEvts = 
      isG4Tree ? 
      ((pNevts > lTree->GetEntries()/30. || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/30.) : pNevts) :
      ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;

    std::cout << "- Processing = " << nEvts  << " events out of " ;
    if (isG4Tree) std::cout << lTree->GetEntries()/30. << std::endl;
    else std::cout << lTree->GetEntries() << std::endl;
    

    for (unsigned ievt(0); ievt<(isG4Tree?nEvts*30:nEvts); ++ievt){//loop on entries
      if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
      else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
      lTree->GetEntry(ievt);

      if (debug){
	if (isG4Tree) {
	  std::cout << "... Processing layer " << volNb << " with " << (*simhitvec).size() << " simhits, weight = " << getWeight(volNb,pVersion) << std::endl;
	  Etot[static_cast<unsigned>(volNb)] = 0;
	}
	else {
	  std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
	}
      }
      
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*simhitvec)[iH];
	unsigned layer = lHit.layer();
	if (isG4Tree && layer==volNb+1) layer = volNb;
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	if (debug>1) {
	  std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}
	double weightedE = lHit.energy()*getWeight(layer,pVersion);
	p_xy[iE][layer]->Fill(posx,posy,weightedE);
	p_time[iE][layer]->Fill(lHit.time());
	Etot[layer] += weightedE;
	Etotal += weightedE;
      }//loop on hits
      if (!isG4Tree){
	Ereco = 0;
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
	  HGCSSRecoHit lHit = (*rechitvec)[iH];
	  if (debug>1) {
	    std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		      << " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	    lHit.Print(std::cout);
	  }
	  Ereco += lHit.energy()*getWeight(lHit.layer(),pVersion);
	}
	p_Ereco[iE]->Fill(Ereco);
	p_nRecHits[iE]->Fill((*rechitvec).size());
	if (debug) std::cout << "... recoE = " << Ereco << std::endl;
      }
      if (!isG4Tree || (isG4Tree && (volNb == nLayers-1))){
	if (debug) std::cout << " -- Filling histograms..." << std::endl;
	double Etmp = 0;
	for (unsigned iL(0);iL<nLayers;++iL){
	  p_Etot[iE][iL]->Fill(Etot[iL]);
	  Etmp += Etot[iL];
	  p_Efrac[iE][iL]->Fill(Etmp/Etotal);
	  Etot[iL] = 0;
	}
	p_Etotal[iE]->Fill(Etotal);
	p_nSimHits[iE]->Fill((*simhitvec).size());
	Etotal = 0;
      }
    }//loop on entries
    for (unsigned iL(0); iL<nLayers; ++iL){
      outputFile->cd();
      p_xy[iE][iL]->Write();
      p_time[iE][iL]->Write();
      p_Etot[iE][iL]->Write();
      p_Efrac[iE][iL]->Write();
    }
    p_Etotal[iE]->Write();
    p_nSimHits[iE]->Write();
    if (!isG4Tree) {
      p_Ereco[iE]->Write();
      p_nRecHits[iE]->Write();
    }
    
    std::cout << " -- Summary of energies: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE]->GetEntries() 
	      << " mean " << p_Etotal[iE]->GetMean() 
	      << " rms " << p_Etotal[iE]->GetRMS() 
	      << " overflows " << p_Etotal[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
	      << std::endl;
    if (!isG4Tree) {
      std::cout << " ---- RecHits: entries " << p_Ereco[iE]->GetEntries() 
		<< " mean " << p_Ereco[iE]->GetMean() 
		<< " rms " << p_Ereco[iE]->GetRMS() 
		<< " overflows " << p_Ereco[iE]->GetBinContent(p_Ereco[iE]->GetNbinsX()+1)
		<< std::endl;
    }
    std::cout << " -- Summary of hits: " << std::endl
	      << " ---- SimHits: entries " << p_nSimHits[iE]->GetEntries() 
	      << " mean " << p_nSimHits[iE]->GetMean() 
	      << " rms " << p_nSimHits[iE]->GetRMS() 
	      << " overflows " << p_nSimHits[iE]->GetBinContent(p_nSimHits[iE]->GetNbinsX()+1)
	      << std::endl;
    if (!isG4Tree) {
      std::cout << " ---- RecHits: entries " << p_nRecHits[iE]->GetEntries() 
		<< " mean " << p_nRecHits[iE]->GetMean() 
		<< " rms " << p_nRecHits[iE]->GetRMS() 
		<< " overflows " << p_nRecHits[iE]->GetBinContent(p_nRecHits[iE]->GetNbinsX()+1)
		<< std::endl;
    }
  }//loop on energies

  outputFile->Write();

  return 0;


}//main
