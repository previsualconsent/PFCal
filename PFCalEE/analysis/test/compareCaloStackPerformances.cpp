#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TMath.h"
#include "TSpline.h"

#include "Math/WrappedTF1.h"
#include "Math/BrentMinimizer1D.h"

#include <map>

#include "HGCSSCaloProperties.hh"

//
int main()
{
  setStyle();

  //size_t versions[]={0,1,2,3,4,5,6,13};
  //size_t versions[]={7,8,9,10,11,12};
  //size_t versions[]={0,20,21};
  //size_t versions[]={20};
  size_t versions[]={100};
  // size_t versions[]={14,15,16,17,18,19};
  size_t nversions( sizeof(versions)/sizeof(size_t) );

  std::vector<TH1F *> stochGr, constGr;
  for(size_t i=0; i<4; i++)
    {
      TString pf(""); pf+=i;
      TH1F *h=new TH1F("stochgr"+pf,";Detector version;(#sigma/E)_{stochastic}",nversions,0,nversions); h->SetMarkerStyle(20+i); h->SetDirectory(0); stochGr.push_back(h);
      h=new TH1F("constgr"+pf,";Detector version;(#sigma/E)_{cte}",nversions,0,nversions);              h->SetMarkerStyle(20+i); h->SetDirectory(0); constGr.push_back(h);
    }  

  for(size_t iv=0; iv<nversions; iv++){
    size_t i=versions[iv];
   
    TString binLabel("CALICE");
    if(i==1) binLabel="Pb CALICE";
    if(i==2) binLabel="Uniform (1.0X_{0})";
    if(i==3) binLabel="Uniform (0.8X_{0})";
    if(i==4) binLabel="Uniform (0.5X_{0})";
    if(i==5) binLabel="Uniform (0.3X_{0})";
    if(i==6) binLabel="JV";
    if(i==7) binLabel ="Si 60#mu";
    if(i==8) binLabel ="Si 80#mu";
    if(i==9) binLabel ="Si 120#mu";
    if(i==10) binLabel="Si 200#mu";
    if(i==11) binLabel="Si 300#mu";
    if(i==12) binLabel="Si 500#mu";
    if(i==13) binLabel="VJ";
    if(i==14) binLabel ="Si 60#mu";
    if(i==15) binLabel ="Si 80#mu";
    if(i==16) binLabel ="Si 120#mu";
    if(i==17) binLabel="Si 200#mu";
    if(i==18) binLabel="Si 300#mu";
    if(i==19) binLabel="Si 500#mu";
    if(i==20) binLabel="HGCal EE";
    if(i==21) binLabel="HGCal EE Si 500#mu";

    TString ver("version_"); ver+=i;

    CaloProperties props(ver);
    props.characterizeCalo();

    for(size_t ialgo=0; ialgo<3; ialgo++)
      {
	stochGr[ialgo]->SetTitle(props.resCurve_[ialgo]->GetTitle());
	stochGr[ialgo]->GetXaxis()->SetBinLabel(iv+1,binLabel);
	stochGr[ialgo]->SetBinContent          (iv+1,props.stochTerms_[ialgo].first);
	stochGr[ialgo]->SetBinError            (iv+1,props.stochTerms_[ialgo].second);
	constGr[ialgo]->SetTitle(props.resCurve_[ialgo]->GetTitle());
	constGr[ialgo]->GetXaxis()->SetBinLabel(iv+1,binLabel);
	constGr[ialgo]->SetBinContent          (iv+1,props.constTerms_[ialgo].first);
	constGr[ialgo]->SetBinError            (iv+1,props.constTerms_[ialgo].second);
      }
  }

  TCanvas *csum=new TCanvas("csum","csum",1200,600);
  csum->Divide(2,1);
  csum->cd(1);
  TLegend *leg=new TLegend(0.15,0.85,0.9,0.95);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetNColumns(3);
  for(Int_t i=0; i<2; i++)
    {
      TPad *p=(TPad *)csum->cd(i+1);
      p->SetTopMargin(0.05);
      p->SetBottomMargin(0.15);

      for(size_t ialgo=0; ialgo<3; ialgo++)
	{	  
	  TH1F *gr=(i==0 ? constGr[ialgo] : stochGr[ialgo]);
	  gr->Draw(ialgo==0 ? "e1" : "e1same");
	  gr->GetYaxis()->SetTitleOffset(1.0);
	  gr->GetYaxis()->SetTitleSize(0.05);
	  gr->GetYaxis()->SetLabelSize(0.04);
	  gr->GetYaxis()->SetRangeUser(0,0.3);
	  gr->GetYaxis()->SetNdivisions(10);
	  gr->GetXaxis()->SetTitleSize(0.05);
	  gr->GetXaxis()->SetLabelSize(0.04);
	  gr->GetXaxis()->SetTitleOffset(1.1);
	  if(i==0) leg->AddEntry(gr,gr->GetTitle(),"p");
	  if(i==0 && ialgo==0) drawHeader();
	}
    }
  csum->cd(1);
  leg->Draw();
  csum->cd();
  csum->Modified();
  csum->Update();
  csum->SaveAs("PLOTS/CaloPerformanceSummary.png");

  return 0;

}//main
