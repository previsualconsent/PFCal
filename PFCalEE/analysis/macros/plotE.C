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

int plotE(){//main


  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  TString scenario[nS] = {"0"};

  const unsigned nV = 1;
  TString version[nV] = {"18"};

  const double Emip = 0.0559;//in MeV

  const unsigned MAX = 8;
  TString type[MAX];
  Float_t sigmaStoch[nS][MAX];
  Float_t sigmaStochErr[nS][MAX];
  Float_t sigmaConst[nS][MAX];
  Float_t sigmaConstErr[nS][MAX];
  //Float_t sigmaNoise[nS][MAX];
  //Float_t sigmaNoiseErr[nS][MAX];

  std::ostringstream saveName;
 
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/scenario_"+scenario[iS]+"/";
      TString plotDir = "../PLOTS/version_"+version[iV]+"/";

      TFile *inputFile = TFile::Open(plotDir+"CalibHistos.root");
      if (!inputFile) {
	std::cout << " -- Error, input file " << plotDir << "/CalibHistos.root cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      
      const unsigned nLayers = 30;
      
      unsigned genEn[]={5,10,25,50,75,100,150,200,300,500};
      //unsigned genEn[]={10};
      const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
      
      unsigned rebin[nGenEn] = {2,2,4,4,8,8,16,16,20,20};
      
      TH1F *p_Efrac[nGenEn][nLayers];
      TH1F *p_Etotal[nGenEn];
      TH1F *p_Ereco[nGenEn];
      
      TH2F *p_meanFrac = new TH2F("p_meanFrac",";layer;gen E (GeV);<E_{layer}>/E_{tot}",30,0,30,nGenEn,0,nGenEn);
      TH2F *p_rmsFrac = new TH2F("p_rmsFrac",";layer;gen E (GeV);#sigma(E_{layer}/E_{tot})",30,0,30,nGenEn,0,nGenEn);
    
      
      TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
      mycL->Divide(6,5);
      gStyle->SetOptStat(0);

      bool isG4File = false;
      for (unsigned iE(0); iE<nGenEn; ++iE){
	
	std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	
	TString eStr;
	eStr += genEn[iE];
	p_meanFrac->GetYaxis()->SetBinLabel(iE+1,eStr);
	p_rmsFrac->GetYaxis()->SetBinLabel(iE+1,eStr);
	
	std::ostringstream lName;
	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_Efrac_" << genEn[iE] << "_" << iL;
	  p_Efrac[iE][iL] = (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!p_Efrac[iE][iL]) {
	    std::cout << " -- ERROR, pointer for histogram Efrac is null for layer: " << iL << ". Exiting..." << std::endl;
	    return 1;
	  }
	  p_meanFrac->SetBinContent(iL+1,iE+1,p_Efrac[iE][iL]->GetMean());
	  p_rmsFrac->SetBinContent(iL+1,iE+1,p_Efrac[iE][iL]->GetRMS());

	  if (iE==5) {
	    mycL->cd(iL+1);
	    p_Efrac[iE][iL]->Draw();
	  }

	}

	if (iE==5){
	  saveName.str("");
	  saveName << plotDir << "/SimEfraction_100GeV";
	  mycL->Update();
	  mycL->Print((saveName.str()+".png").c_str());
	  mycL->Print((saveName.str()+".pdf").c_str());
	}

	lName.str("");
	lName << "p_Etotal_" << genEn[iE];
	p_Etotal[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!p_Etotal[iE]){
	  std::cout << " -- ERROR, pointer for histogram Etotal is null. Exiting..." << std::endl;
	  return 1;
	}
    
	std::cout << " --- Sim E = entries " << p_Etotal[iE]->GetEntries() 
		  << " mean " << p_Etotal[iE]->GetMean() 
		  << " rms " << p_Etotal[iE]->GetRMS() 
		  << " overflows " << p_Etotal[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
		  << std::endl;
	p_Etotal[iE]->Rebin(rebin[iE]);

	lName.str("");
	lName << "p_Ereco_" << genEn[iE];
	p_Ereco[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!p_Ereco[iE]){
	  std::cout << " -- ERROR, pointer for histogram Ereco is null. Running on G4 file before digitizer." << std::endl;
	  isG4File = true;
	}
	else {
	  std::cout << " --- Reco E = entries " << p_Ereco[iE]->GetEntries() 
		    << " mean " << p_Ereco[iE]->GetMean() 
		    << " rms " << p_Ereco[iE]->GetRMS() 
		    << " overflows " << p_Ereco[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
		    << std::endl;
	  p_Ereco[iE]->Rebin(rebin[iE]);
	}


      }//loop on energies

      //draw energy fractions
      TCanvas *mycF = new TCanvas("mycF","mycF",1500,750);
      gStyle->SetOptStat(0);
      mycF->Divide(2,1);
      mycF->cd(1);
      p_meanFrac->Draw("colz");
      mycF->cd(2);
      p_rmsFrac->Draw("colz");


      saveName.str("");
      saveName << plotDir << "/SimEfractionIntegrated";
      mycF->Update();
      mycF->Print((saveName.str()+".png").c_str());
      mycF->Print((saveName.str()+".pdf").c_str());

      return 1;

      TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
      mycE->Divide(5,2);
  
      gStyle->SetOptStat(0);

      //draw calibration curves
      TGraphErrors *calib = new TGraphErrors();
      calib->SetName("calib");
      calib->SetMarkerStyle(21);
      calib->SetTitle("");
      TGraphErrors *reso = (TGraphErrors *) calib->Clone("reso");
      TGraphErrors *calibFit = (TGraphErrors *) calib->Clone("calibFit");
      TGraphErrors *resoFit = (TGraphErrors *) calib->Clone("resoFit");
      TGraphErrors *calibReco = (TGraphErrors *) calib->Clone("calibReco");
      calibReco->SetMarkerStyle(22);
      calibReco->SetMarkerColor(6);
      calibReco->SetLineColor(6);
      TGraphErrors *resoReco = (TGraphErrors *) calibReco->Clone("resoReco");
      TGraphErrors *calibRecoFit = (TGraphErrors *) calibReco->Clone("calibRecoFit");
      TGraphErrors *resoRecoFit = (TGraphErrors *) calibReco->Clone("resoRecoFit");
  
      //simhits
      for (unsigned iE(0); iE<nGenEn; ++iE){
	std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	//plot total E
	mycE->cd(iE+1);
	gStyle->SetOptFit(0);
	p_Etotal[iE]->GetXaxis()->SetRangeUser(p_Etotal[iE]->GetMean()-5*p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMean()+5*p_Etotal[iE]->GetRMS());
	p_Etotal[iE]->Draw();
	p_Etotal[iE]->Fit("gaus");
	TF1 *fitResult = p_Etotal[iE]->GetFunction("gaus");
	TLatex lat;
	char buf[500];
	sprintf(buf,"<E> = %3.1f MeV",p_Etotal[iE]->GetMean());
	lat.DrawLatex(p_Etotal[iE]->GetMean()-5*p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMaximum()*0.9,buf);
	sprintf(buf,"RMS = %3.1f MeV",p_Etotal[iE]->GetRMS());
	lat.DrawLatex(p_Etotal[iE]->GetMean()-5*p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMaximum()*0.8,buf);
	sprintf(buf,"<Efit> = %3.1f MeV",fitResult->GetParameter(1));
	lat.DrawLatex(p_Etotal[iE]->GetMean()+p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMaximum()*0.9,buf);
	sprintf(buf,"RMSfit = %3.1f MeV",fitResult->GetParameter(2));
	lat.DrawLatex(p_Etotal[iE]->GetMean()+p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMaximum()*0.8,buf);


	Int_t np=calib->GetN();
	calib->SetPoint(np,genEn[iE],p_Etotal[iE]->GetMean());
	calib->SetPointError(np,0.0,p_Etotal[iE]->GetMeanError());
	reso->SetPoint(np,1/sqrt(genEn[iE]),p_Etotal[iE]->GetRMS()/p_Etotal[iE]->GetMean());
	reso->SetPointError(np,0,p_Etotal[iE]->GetRMSError()/p_Etotal[iE]->GetMean());
	calibFit->SetPoint(np,genEn[iE],fitResult->GetParameter(1));
	calibFit->SetPointError(np,0.0,fitResult->GetParError(1));
	resoFit->SetPoint(np,1/sqrt(genEn[iE]),fitResult->GetParameter(2)/fitResult->GetParameter(1));
	resoFit->SetPointError(np,0,fitResult->GetParError(2)/fitResult->GetParameter(1));

      }//loop on energies

      saveName.str("");
      saveName << plotDir << "/SimG4Etotal";
      mycE->Update();
      mycE->Print((saveName.str()+".png").c_str());
      mycE->Print((saveName.str()+".pdf").c_str());

      if (!isG4File){
	//recohits
	for (unsigned iE(0); iE<nGenEn; ++iE){
	  //plot reco E
	  std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	  mycE->cd(iE+1);
	  gStyle->SetOptFit(0);
	  p_Ereco[iE]->GetXaxis()->SetRangeUser(p_Ereco[iE]->GetMean()-5*p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMean()+5*p_Ereco[iE]->GetRMS());
	  p_Ereco[iE]->Draw();
	  p_Ereco[iE]->Fit("gaus");
	  TF1 *fitResult = p_Ereco[iE]->GetFunction("gaus");
	  TLatex lat;
	  char buf[500];
	  sprintf(buf,"<E> = %3.1f MIPs",p_Ereco[iE]->GetMean());
	  lat.DrawLatex(p_Ereco[iE]->GetMean()-5*p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.9,buf);
	  sprintf(buf,"RMS = %3.1f MIPs",p_Ereco[iE]->GetRMS());
	  lat.DrawLatex(p_Ereco[iE]->GetMean()-5*p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.8,buf);
	  sprintf(buf,"<Efit> = %3.1f MIPs",fitResult->GetParameter(1));
	  lat.DrawLatex(p_Ereco[iE]->GetMean()+p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.9,buf);
	  sprintf(buf,"RMSfit = %3.1f MIPs",fitResult->GetParameter(2));
	  lat.DrawLatex(p_Ereco[iE]->GetMean()+p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.8,buf);
      
      
	  Int_t np=calibReco->GetN();
	  calibReco->SetPoint(np,genEn[iE],p_Ereco[iE]->GetMean()*Emip);
	  calibReco->SetPointError(np,0.0,p_Ereco[iE]->GetMeanError()*Emip);
	  resoReco->SetPoint(np,1/sqrt(genEn[iE]),p_Ereco[iE]->GetRMS()/p_Ereco[iE]->GetMean());
	  resoReco->SetPointError(np,0,p_Ereco[iE]->GetRMSError()/p_Ereco[iE]->GetMean());
	  calibRecoFit->SetPoint(np,genEn[iE],fitResult->GetParameter(1)*Emip);
	  calibRecoFit->SetPointError(np,0.0,fitResult->GetParError(1)*Emip);
	  resoRecoFit->SetPoint(np,1/sqrt(genEn[iE]),fitResult->GetParameter(2)/fitResult->GetParameter(1));
	  resoRecoFit->SetPointError(np,0,fitResult->GetParError(2)/fitResult->GetParameter(1));
      
	}//loop on energies

	saveName.str("");
	saveName << plotDir << "/DigiEreco";
	mycE->Update();
	mycE->Print((saveName.str()+".png").c_str());
	mycE->Print((saveName.str()+".pdf").c_str());
      }

      //draw calib
      const unsigned nCanvas = 4;  
      TCanvas *myc[nCanvas];
      for (unsigned iC(0);iC<nCanvas;++iC){
	std::ostringstream lName;
	lName << "myc" << iC;
	myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
      }


      const unsigned imax = isG4File ? 4 : 8;

      for(unsigned i=0; i<imax ; i++)
	{
	  myc[i%4]->cd();
	  type[i] = i==0 ? "calib" : i==1 ? "calibFit" : i==2 ? "reso" : "resoFit";
	  if (!isG4File && i>3) type[i] = i==4 ? "calibReco" : i==5 ? "calibRecoFit" : i==6 ? "resoReco" : "resoRecoFit";

	  std::cout << "- Processing type : " << type[i] << std::endl;

	  TGraphErrors * gr =( i==0 ? calib : i==1 ? calibFit : i==2 ? reso : resoFit);
	  if (!isG4File && i>3) gr = i==4 ? calibReco : i==5 ? calibRecoFit : i==6 ? resoReco : resoRecoFit;

	  gr->Draw(i<4? "ap" : "p");
	  gr->GetYaxis()->SetRangeUser(0,gr->GetYaxis()->GetXmax());
	    
	  if(i<2|| (!isG4File && (i==4 || i==5))) { 
	    gr->GetXaxis()->SetTitle("Beam energy [GeV]");
	    gr->GetYaxis()->SetTitle("Average energy deposited [MeV]"); }
	  else { 
	    gr->GetXaxis()->SetTitle("1/#sqrt{Beam energy} [1/#sqrt{GeV}]"); 
	    gr->GetYaxis()->SetTitle("Relative energy resolution");     
	  }
	  char buf[500];
	  if(i<2 || (!isG4File && (i==4 || i==5))) {
	    TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
	    if (i<4) fitFunc->SetLineColor(1);
	    else fitFunc->SetLineColor(6);
	    gr->Fit(fitFunc,"RME");
	    TLatex lat;
	    if (i>3) lat.SetTextColor(6);
	    else lat.SetTextColor(1);
	    sprintf(buf,"<E> #propto a + b #times E ");
	    if (i<2) lat.DrawLatex(70,gr->GetYaxis()->GetXmax()*0.9,buf);
	    sprintf(buf,"a = %3.3f #pm %3.3f ",fitFunc->GetParameter(0),fitFunc->GetParError(0));
	    lat.DrawLatex(70+i/2*150,gr->GetYaxis()->GetXmax()*(0.8-i/2*0.25),buf);
	    sprintf(buf,"b = %3.3f #pm %3.3f",fitFunc->GetParameter(1),fitFunc->GetParError(1));
	    lat.DrawLatex(70+i/2*150,gr->GetYaxis()->GetXmax()*(0.7-i/2*0.25),buf);
	  }
	  else
	    {
	      //TF1 *fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1]+[2]*x*x*x*x)",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
	      TF1 *fitFunc2;
	      fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
	      //if (i<4) fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
	      //else fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),0.2);//gr->GetXaxis()->GetXmax());
	      fitFunc2->SetParameter(0,0.2);
	      fitFunc2->SetParLimits(0,0,1);
	      fitFunc2->SetParameter(1,0);
	      fitFunc2->SetParLimits(1,0,1);
	      if (i<4) {
		//fitFunc2->SetParameter(2,0.);
		//fitFunc2->SetParLimits(2,0,0);
		fitFunc2->SetLineColor(1);
	      }
	      else {
		//fitFunc2->SetLineColor(6);
		//fitFunc2->SetParameter(2,0.);
		fitFunc2->SetParLimits(2,0,0);
	      }
	      gr->Fit(fitFunc2,"RME");
	      sigmaStoch[iS][i] = sqrt(fitFunc2->GetParameter(0));
	      sigmaStochErr[iS][i] = fitFunc2->GetParError(0)/(2*sigmaStoch[iS][i]);
	      sigmaConst[iS][i] = sqrt(fitFunc2->GetParameter(1));
	      sigmaConstErr[iS][i] = fitFunc2->GetParError(1)/(2*sigmaConst[iS][i]);
	      //sigmaNoise[iS][i] = sqrt(fitFunc2->GetParameter(2));
	      //sigmaNoiseErr[iS][i] = fitFunc2->GetParError(2)/(2*sigmaNoise[iS][i]);
	      TLatex lat;
	      if (i>3) lat.SetTextColor(6);
	      else lat.SetTextColor(1);
	      sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");
	      if (i<4) lat.DrawLatex(0.05,gr->GetYaxis()->GetXmax()*0.9,buf);
	      sprintf(buf,"s=%3.3f #pm %3.3f",sigmaStoch[iS][i],sigmaStochErr[iS][i]);
	      lat.DrawLatex(i/2*0.125-0.075,gr->GetYaxis()->GetXmax()*(0.9-i/2*0.1-i/4*0.3),buf);
	      sprintf(buf,"c=%3.3f #pm %3.3f",sigmaConst[iS][i],sigmaConstErr[iS][i]);
	      lat.DrawLatex(i/2*0.125-0.075,gr->GetYaxis()->GetXmax()*(0.8-i/2*0.1-i/4*0.3),buf);
	      if (i>3){
		//sprintf(buf,"n=%3.3f #pm %3.3f",sigmaNoise[iS][i],sigmaNoiseErr[iS][i]);
		//lat.DrawLatex(i/2*0.125-0.075,gr->GetYaxis()->GetXmax()*(0.7-i/2*0.1-i/4*0.3),buf);
	      }
	    }
	  myc[i%4]->Update();
	  myc[i%4]->Print(plotDir+"/"+type[i]+".pdf");
	  myc[i%4]->Print(plotDir+"/"+type[i]+".png");
	}


    }//loop on scenarios

  }//loop on versions

  for (unsigned iV(0); iV<nV;++iV){
    std::cout << "version " << version[iV] << std::endl;
    std::cout << "scenario & type & sigmaStoch & sigmaConst"
      //<< " & sigmaNoise"
	      << " \\\\ \n" 
	      << "\\hline\n";
    for (unsigned iS(0); iS<nS;++iS){
      for (unsigned i(2); i<MAX;++i){
	if (i==4 || i==5) continue;
	std::cout << scenario[iS] << " & " << type[i] << " & " <<  std::setprecision(3)
		  << "$" << sigmaStoch[iS][i] << "\\pm" << sigmaStochErr[iS][i] << "$ & "
		  << "$" << sigmaConst[iS][i] << "\\pm" << sigmaConstErr[iS][i] << "$"
	  //<< " & $" << sigmaNoise[iS][i] << "\\pm" << sigmaNoiseErr[iS][i] << "$"
		  << "\n";
      }
      std::cout <<"\\hline\n";
    }
  }
  
  return 0;
  
  
}//main
