#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"


template <class T>
void extractParameterFromStr(std::string aStr,T vec){ 
  if (aStr == "") return;
  std::vector<std::string> layVec;
  boost::split( layVec, aStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for parameter given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      exit(1);
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      exit(1);
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      std::istringstream(lPair[1])>>vec[iL];
    }
  }//loop on elements
}


void simToDigi(HGCSSRecoHit & aHit, 
	       const double & aNoise, 
	       const double & aMipE, 
	       const unsigned & aMipToADC){
  //convert to MIP
  double oldE = aHit.energy()/aMipE;
  double newE = oldE+aNoise;
  aHit.noiseFraction(newE > 0 ? 
		     (fabs(aNoise) < oldE ? aNoise/oldE : 1) : -1);
  aHit.energy(newE > 0 ? newE : 0);
  //round to the lower integer, not nearer...
  aHit.adcCounts(static_cast<unsigned>(aHit.energy()*aMipToADC));
}

bool digiToReco(HGCSSRecoHit & aRecHit,
		const unsigned & aMipToADC,
		const unsigned & aThreshold){
  aRecHit.energy(aRecHit.adcCounts()*1.0/aMipToADC);
  if (aRecHit.adcCounts() > aThreshold) 
    return true;
  return false;
}

int main(int argc, char** argv){//main  

  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  const unsigned nReqA = 8;
  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < nReqA) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"<< std::endl
	      << "<full path to input file>"<< std::endl
	      << "<full path to output file>"<< std::endl
	      << "<granularities \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<threshold (in ADC counts) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
	      << "<mipToADC (in ADC counts/MIP)> " << std::endl
	      << std::endl
	      << "<optional: randomSeed (default=0)> "  << std::endl
	      << "<optional: debug (default=0)>" << std::endl
	      << "<optional: save sim hits (default=1)> " << std::endl
	      << "<optional: save digi hits (default=0)> " << std::endl
	      << "<optional: mipEnergy (default=0.0559 MeV)>" << std::endl
	      << "<optional: maxTime (default=2ns)>" << std::endl
	      << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string inFilePath = argv[2];
  std::string outFilePath = argv[3];
  std::string granulStr = argv[4];
  std::string noiseStr = argv[5];
  std::string threshStr = argv[6];
  unsigned pMipToADC = 50;
  std::istringstream(argv[7])>>pMipToADC;

  std::cout << " ----------------------------------------" << std::endl
	    << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Processing " << pNevts << " events." << std::endl
	    << " -- Granularities: " << granulStr << std::endl
	    << " -- noise: " << noiseStr << std::endl
	    << " -- thresholds: " << threshStr << std::endl
	    << " -- Mip to ADC conversion: " << pMipToADC << " ADC counts / MIP." << std::endl;
	    
  double pMipEnergy = 0.0559; //from muon deposits at 25 GeV, in version_20...
  double pMaxTime = 2.0;//ns
  unsigned debug = 0;
  unsigned pSeed = 0;
  bool pSaveDigis = 0;
  bool pSaveSims = 1;
  if (nPar > nReqA) std::istringstream(argv[nReqA+1])>>pSeed;
  if (nPar > nReqA+1) {
    debug = atoi(argv[nReqA+1]);
    std::cout << " -- DEBUG output is set to " << debug << std::endl;
  }
  if (nPar > nReqA+2) std::istringstream(argv[nReqA+2])>>pSaveDigis;
  if (nPar > nReqA+3) std::istringstream(argv[nReqA+3])>>pSaveSims;
  if (nPar > nReqA+4) std::istringstream(argv[nReqA+4])>>pMipEnergy;
  if (nPar > nReqA+5) std::istringstream(argv[nReqA+5])>>pMaxTime;
  
  std::cout << " -- Random seed will be set to : " << pSeed << std::endl;
  if (pSaveDigis) std::cout << " -- DigiHits are saved." << std::endl;
  if (!pSaveSims) std::cout << " -- SimHits are not saved." << std::endl;
  std::cout << " -- Mip energy: " << pMipEnergy << " MeV." << std::endl
	    << " -- Max time integrated: " << pMaxTime << " ns." << std::endl
	    << " ----------------------------------------" << std::endl;
  

  const unsigned nLayers = N_LAYERS;
  unsigned granularity[nLayers];
  double pNoiseInMips[nLayers];
  unsigned pThreshInADC[nLayers];
  for (unsigned iL(0); iL<nLayers; ++iL){
    granularity[iL] = 1;
    pNoiseInMips[iL] = 0.1;
    pThreshInADC[iL] = static_cast<unsigned>(5*pNoiseInMips[iL]*pMipToADC+0.5);
  }

  extractParameterFromStr<unsigned[nLayers]>(granulStr,granularity);
  extractParameterFromStr<double[nLayers]>(noiseStr,pNoiseInMips);
  extractParameterFromStr<unsigned[nLayers]>(threshStr,pThreshInADC);

  unsigned nbCells = 0;

  std::cout << " -- Granularities and noise are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL] << " mips, " << pThreshInADC[iL] << " adc - ";
    if (iL%5==4) std::cout << std::endl;
    
    nbCells += N_CELLS_XY_MAX/(granularity[iL]*granularity[iL]);
  }
        
  std::cout << " -- Total number of cells = " << nbCells << std::endl;

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(pSeed);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  /////////////////////////////////////////////////////////////
  //output
  /////////////////////////////////////////////////////////////

  std::ostringstream outputStr;
  outputStr << outFilePath << "/DigiPFcal" ;
  // outputStr << "_" << granulStr << "_";
  // if (noiseStr.size()>0) outputStr << noiseStr << "_";
  // else outputStr << pNoiseInMips[0] << "_";
  // if (threshStr.size()>0) outputStr << threshStr << "_";
  // else outputStr << pThreshInADC[0] << "_";
  // outputStr << pMipToADC << "ADCPerMip";
  // //outputStr << "_" << pMaxTime << "ns";
  if (pSaveDigis)  outputStr << "_withDigiHits";
  if (!pSaveSims)  outputStr << "_withoutSimHits";
  outputStr << ".root";
  
  TFile *outputFile = TFile::Open(outputStr.str().c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outputStr.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- File will be saved as " << outputStr.str() << std::endl;
  }

  TTree *outputTree = new TTree("RecoTree","HGC Standalone simulation reco tree");
  HGCSSSimHitVec lSimHits;
  HGCSSRecoHitVec lDigiHits;
  HGCSSRecoHitVec lRecoHits;
  unsigned maxSimHits = 0;
  unsigned lEvent = 0;
  outputTree->Branch("event",&lEvent);
  if (pSaveSims) outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&lSimHits);
  if (pSaveDigis) outputTree->Branch("HGCSSDigiHitVec","std::vector<HGCSSRecoHit>",&lDigiHits);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);

  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-2,2);

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::string inputStr = inFilePath + "/PFcal.root";
  TFile *inputFile = TFile::Open(inputStr.c_str());
  
  if (!inputFile) {
    std::cout << " -- Error, input file " << inputStr << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TTree *inputTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!inputTree){
    std::cout << " -- Error, tree HGCSSTree  cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //input tree
  /////////////////////////////////////////////////////////////

  float event;
  float volNb;
  std::vector<HGCSSSimHit> * hitvec = 0;

  inputTree->SetBranchAddress("event",&event);
  inputTree->SetBranchAddress("volNb",&volNb);
  inputTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
    
  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////


  const unsigned nEvts = (pNevts > inputTree->GetEntries()/30. || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()/30.) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries()/30. << std::endl;

  //create map used to assemble hits per event.
  std::map<unsigned,HGCSSRecoHit> lHitMap;
  std::pair<std::map<unsigned,HGCSSRecoHit>::iterator,bool> isInserted;

  for (unsigned ientry(0); ientry<nEvts*30; ++ientry){//loop on entries

    unsigned ievt =  ientry/30;

    inputTree->GetEntry(ientry);
    lEvent = event;
    unsigned layer = volNb;
    
    if (debug>1) {
      std::cout << " **DEBUG** Processing evt " << ievt << " layer " << layer << std::endl;
    }
    else if (debug == 1 && layer==0) std::cout << "... Processing event: " << ievt << std::endl;
    else if (ievt%50 == 0 && layer==0) std::cout << "... Processing event: " << ievt << std::endl;
    

    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*hitvec)[iH];
      lSimHits.push_back(lHit);

      //C-AMM: TO DO ?
      //discard simhits not in the right time window
      //but needs to take into account propagation of the particle from interaction point...
      //for the moment accept everything.
      //In CALICE, shaper time much larger means everything is integrated.
      //if (lHit.time() > pMaxTime) continue;

      HGCSSRecoHit lRecHit(lHit,granularity[layer]);
      //double posx = lHit.get_x();
      //double posy = lHit.get_y();

      isInserted = lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lRecHit.fullcellid(),lRecHit));
      if (!isInserted.second) isInserted.first->second.Add(lHit);

    }//loop on input simhits

    if (debug>1) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << " recohits = " << lHitMap.size() << std::endl;
    }

    //create digihits, everywhere to have also pure noise.
    for (unsigned iX(0); iX<SIZE_X/(CELL_SIZE_X*granularity[layer])/2;++iX){
      for (unsigned iY(0); iY<SIZE_Y/(CELL_SIZE_Y*granularity[layer])/2;++iY){
	HGCSSRecoHit lNoiseHit;
	lNoiseHit.layer(layer);
	lNoiseHit.encodeCellId(true,true,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
	lNoiseHit.encodeCellId(true,false,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
	lNoiseHit.encodeCellId(false,true,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
	lNoiseHit.encodeCellId(false,false,iX,iY,granularity[layer]);
	lHitMap.insert(std::pair<unsigned,HGCSSRecoHit>(lNoiseHit.fullcellid(),lNoiseHit));
      }
    }

    
    if (layer == nLayers-1){
      //save event
      std::map<unsigned,HGCSSRecoHit>::iterator lIter = lHitMap.begin();
      lDigiHits.reserve(lHitMap.size());

      for (; lIter != lHitMap.end(); ++lIter){
	HGCSSRecoHit & lHit = lIter->second;
	//convert to MIP and add noise.
	double lNoise = lRndm->Gaus(0,pNoiseInMips[layer]);
	simToDigi(lHit,lNoise,pMipEnergy,pMipToADC);
	lDigiHits.push_back(lHit);
	p_noise->Fill(lNoise);
      }
      //apply threshold
      lRecoHits.reserve(lDigiHits.size());
      for (unsigned iH(0); iH<lDigiHits.size(); ++iH){
	//copy digihit to modify it
	HGCSSRecoHit lHit = lDigiHits[iH];
	if (digiToReco(lHit,pMipToADC,pThreshInADC[layer]))
	  lRecoHits.push_back(lHit);
      }

      if (debug) {
	std::cout << " **DEBUG** sim-digi-reco hits = " << lSimHits.size() << "-" << lDigiHits.size() << "-" << lRecoHits.size() << std::endl;
      }

      
      outputTree->Fill();
      //reserve necessary space and clear vectors.
      if (lSimHits.size() > maxSimHits) {
	maxSimHits = 2*lSimHits.size();
	std::cout << " -- INFO: event " << ievt << " maxSimHits updated to " << maxSimHits << std::endl;
      }
      lSimHits.clear();
      lDigiHits.clear();
      lRecoHits.clear();
      lHitMap.clear();
      lSimHits.reserve(maxSimHits);
    }
    
  }//loop on entries

  outputFile->cd();
  outputTree->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;

}//main
