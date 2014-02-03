#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TMath.h"

#include "HGCSSSimHit.hh"
#include "HGCSSParameters.hh"

using namespace std;

//
void printHelp()
{
  printf("-h       --> print this\n");
  printf("-i       --> input directory with the simulations\n");
  printf("-v       --> version\n");
  printf("-e       --> csv list with energies\n");
  printf("-d       --> activate printouts\n");
  printf("command line example: plotXY -i /afs/cern.ch/user/a/amagnan/SLHC/PFCal/PFCalEE -d -v 20 -e 5\n");
}

//
int main(int argc, char** argv){//main  

  bool debug(false);
  TString inputDir("/afs/cern.ch/user/a/amagnan/SLHC/PFCal/PFCalEE");
  TString version("version_");
  std::vector<TString> genEn;
  TString gun("e-");

  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("-h")!=string::npos) { printHelp();    return 0; }
      if(arg.find("-d")!=string::npos) { debug=true; }
      if(arg.find("-g")!=string::npos) { gun=argv[i+1]; i++; cout << "Particle gun is : " << gun << endl; }
      if(arg.find("-i")!=string::npos) { inputDir=argv[i+1]; i++; gSystem->ExpandPathName(inputDir); cout << "Input directory: " << inputDir << endl; }     
      if(arg.find("-v")!=string::npos) { version+=argv[i+1]; i++; cout << "Version: " << version << endl; }
      if(arg.find("-e")!=string::npos) { 
	TString csvListStr(argv[i+1]); 
	i++; 
	TObjArray *csvList = csvListStr.Tokenize(",");
	for(Int_t itkn=0; itkn<csvList->GetEntriesFast(); itkn++)
	    genEn.push_back( csvList->At(itkn)->GetName() );
	cout << "Will study the following energies: " << csvListStr << endl;
      } 
    }
  if(genEn.size()==0 || version=="version_") { printHelp(); return 0; }


  const unsigned nGenEn(genEn.size());     
  TH2F *sensorXY        = new TH2F("sensorxy",";x(mm);y(mm);Energy [MeV]",   80,-100,100,80,-100,100); sensorXY->Sumw2();        sensorXY->SetDirectory(0);
  TH1F *sensorEnHits    = new TH1F("sensorhitsen",";Energy [MeV];Events",    100,0, 1);                sensorEnHits->Sumw2();    sensorEnHits->SetDirectory(0);
  TH1F *sensorTotEn     = new TH1F("sensoren",";Energy [MeV];Events",        100,0, 10);               sensorTotEn->Sumw2();     sensorTotEn->SetDirectory(0);
  TH1F *sensorEnHitsvsR = new TH1F("sensorenhitsvsr",";#rho(mm); Energy [MeV]", 100,0, 250);            sensorEnHitsvsR->Sumw2(); sensorEnHitsvsR->SetDirectory(0);
  TH1F *sensorEnHitsvsX = new TH1F("sensorenhitsvsx",";x(mm); Energy [MeV]", 100,0, 250);               sensorEnHitsvsX->Sumw2(); sensorEnHitsvsX->SetDirectory(0);
  TH1F *sensorEnHitsvsY = new TH1F("sensorenhitsvsy",";y(mm); Energy [MeV]", 100,0, 250);               sensorEnHitsvsY->Sumw2(); sensorEnHitsvsY->SetDirectory(0);

  TFile *outputFile = TFile::Open("SimHistos_"+gun+".root","RECREATE");

  /*
  //prepare output

  TH2F *p_xy[nGenEn][nLayers];
  TH1F *p_Etot[nGenEn][nLayers];
  double Emax[nGenEn];
  */

  for (unsigned iE(0); iE<nGenEn; ++iE){

    std::cout << "[Starting " << genEn[iE] << " GeV] " << std::endl;

    //check inputs
    TString fName(inputDir+version+"/"+gun+"/e_"); fName += genEn[iE]; fName += "/PFcal.root";
    TFile *inputFile = TFile::Open(fName);
    if (!inputFile || inputFile->IsZombie()) {
      std::cout << "[Error] input file " << fName << " cannot be opened or corrupted" << std::endl;
      continue;
    }
    TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
    if (!lTree){
      std::cout << "[Error] ntuple HGCSSTree can't be found" << std::endl;
      inputFile->Close();
      continue;
    }

    //key is (Si pad)
    std::map< unsigned, TH2F *> map_sensorXY;
    std::map< unsigned, TH1F *> map_sensorEnHits, map_sensorTotEn, map_sensorEnHitsvsR,  map_sensorEnHitsvsX,  map_sensorEnHitsvsY;

    //decode tree
    Float_t event, volNb, volX0, volX0trans,  den;
    std::vector<HGCSSSimHit> * hitvec = 0;
    lTree->SetBranchAddress("event",&event);
    lTree->SetBranchAddress("volNb",&volNb);
    lTree->SetBranchAddress("volX0",&volX0);
    lTree->SetBranchAddress("volX0trans",&volX0trans);
    lTree->SetBranchAddress("den",&den);
    lTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);

    std::cout << "   Processing = " << lTree->GetEntriesFast()  << " entries" << std::endl;
    lTree->GetEntry(0);
    Float_t curEvent(event),x0trans(event);
    Int_t evtCtr(1);
    for (Int_t ievt(0); ievt<lTree->GetEntriesFast(); ++ievt){//loop on entries
      
      lTree->GetEntry(ievt);      
      unsigned siLayer = volNb;

      x0trans += volX0trans;
      if(curEvent!=event)
       	{
	  evtCtr++;
	  curEvent=event;
	  x0trans=0;
	}

      //check if histos are there
      if(map_sensorXY.find(siLayer)==map_sensorXY.end())
	{
	  TString pf("_"); pf+=siLayer;
	  map_sensorXY[siLayer]        = (TH2F *)sensorXY->Clone("sensorxy_"+pf);               map_sensorXY[siLayer]->SetDirectory(0);
	  map_sensorEnHits[siLayer]    = (TH1F *)sensorEnHits->Clone("sensorhitsen_"+pf);       map_sensorEnHits[siLayer]->SetDirectory(0);
	  map_sensorTotEn[siLayer]     = (TH1F *)sensorTotEn->Clone("sensoren_"+pf);            map_sensorTotEn[siLayer]->SetDirectory(0);
	  map_sensorEnHitsvsR[siLayer] = (TH1F *)sensorEnHitsvsR->Clone("sensorenhitsvsr_"+pf); map_sensorEnHitsvsR[siLayer]->SetDirectory(0);
	  map_sensorEnHitsvsX[siLayer] = (TH1F *)sensorEnHitsvsX->Clone("sensorenhitsvsx_"+pf); map_sensorEnHitsvsX[siLayer]->SetDirectory(0);
	  map_sensorEnHitsvsY[siLayer] = (TH1F *)sensorEnHitsvsY->Clone("sensorenhitsvsy_"+pf); map_sensorEnHitsvsY[siLayer]->SetDirectory(0);
	}
      map_sensorTotEn[siLayer]->Fill(den);

      Float_t totEnInHits(0);
      for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*hitvec)[iH];
	lHit.layer(siLayer);
	double hitEn = lHit.energy();
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	totEnInHits += hitEn;
	map_sensorXY[siLayer]->Fill(posx,posy,hitEn);
	map_sensorEnHits[siLayer]->Fill(hitEn);
	map_sensorEnHitsvsR[siLayer]->Fill(TMath::Sqrt(TMath::Power(posx,2)+TMath::Power(posy,2)),hitEn);
	map_sensorEnHitsvsX[siLayer]->Fill(fabs(posx),hitEn);
	map_sensorEnHitsvsY[siLayer]->Fill(fabs(posy),hitEn);


	if (debug) {
	  std::cout << " --  Hit " << iH << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}	
      }
      
    }//loop on entries


    outputFile->cd();
    outputFile->mkdir("E"+genEn[iE])->cd();
    for(std::map< unsigned, TH2F *>::iterator it=map_sensorXY.begin(); it!=map_sensorXY.end(); it++)
      {
	unsigned key(it->first);
	map_sensorXY[key]->Scale(1./evtCtr);        map_sensorXY[key]->Write();
	map_sensorEnHits[key]->Scale(1./evtCtr);    map_sensorEnHits[key]->Write();
	map_sensorTotEn[key]->Scale(1./evtCtr);     map_sensorTotEn[key]->Write();
	map_sensorEnHitsvsR[key]->Scale(1./evtCtr); map_sensorEnHitsvsR[key]->Write();
	map_sensorEnHitsvsX[key]->Scale(1./evtCtr); map_sensorEnHitsvsX[key]->Write();
	map_sensorEnHitsvsY[key]->Scale(1./evtCtr); map_sensorEnHitsvsY[key]->Write();
      }
  }


  outputFile->Close();

    /*
      

    for (unsigned iL=0; iL<nLayers; ++iL){
      TString pf("pad_"); pf+= iL;
      



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
    
    
    
    const unsigned nEvts = (pNevts > lTree->GetEntries()/30. || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()/30.) : pNevts;


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
*/
  
}//main
