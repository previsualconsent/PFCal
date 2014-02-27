#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"

#include "TMath.h"
#include "Math/Vector3D.h"
#include "TFitResult.h"

#include "HGCSSSimHit.hh"
#include "HGCSSParameters.hh"


typedef ROOT::Math::XYZVectorF Position;

struct hit
{
    Position pos;
    double energy;
    hit(Position new_pos,double new_energy):
        pos(new_pos),
        energy(new_energy)
    {}
};

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

class CaloDimen
{
    public:

        CaloDimen(){}
        CaloDimen(std::string version)
        {

            if (version == "version_20")
            {
                m_layers.push_back(10);
                m_layers.push_back(10);
                m_layers.push_back(10);
                m_layer_thick.push_back(6.8f);
                m_layer_thick.push_back(8.5f);
                m_layer_thick.push_back(10.8f);
            }
        }

        double getZ(unsigned layer)
        {
            float z = 0;

            for(int i=0; i < m_layers.size(); i++)
                z+=m_layers[i]*m_layer_thick[i];
            z= -z/2.0f;
            int sect;
            for(sect = 0; layer > m_layers[sect] && sect < m_layers.size(); sect++)
            {
                z+= m_layers[sect]*m_layer_thick[sect];
                layer-=m_layers[sect];
            }
            z+=m_layer_thick[sect]*layer;
            return z;
        }
        double getThick(unsigned layer)
        {
            int sect;
            for(sect = 0; layer > m_layers[sect] && sect < m_layers.size();sect++)
                layer-=m_layers[sect];
            return  m_layer_thick[sect];
        }
    private:
        std::vector<int> m_layers;
        std::vector<float> m_layer_thick;


};
class EventShowerCone
{
    public:
        EventShowerCone(Position max, Position vertex,std::string version, int event, int genEn, float radius)
        {
            m_version = version;
            m_dim = CaloDimen(m_version);

            m_event = event;
            m_genEn = genEn;
            m_radius = radius;

            m_max = layer2z(max);
            m_vertex = vertex;
            m_slope = m_vertex - m_max;
            m_slope /= m_slope.Z();
            //printPos(max,"max");
            //printPos(m_max,"fixed max");
            //printPos(vertex,"vertex");
            //printPos(m_vertex,"fixed vertex");
            //printPos(m_slope,"slope");


            for(int i=0; i<N_LAYERS; i++)
            {
                std::ostringstream lName;
                lName << "layer";
                if(i<10) lName << "0";
                lName << i;
                m_plots.push_back(new TH2F(lName.str().c_str(),(lName.str()+";x (mm);y (mm)").c_str(),80,-100,100,80,-100,100));
                lName.str("");
            }
        }
        ~EventShowerCone()
        {
            for(std::vector<TH2F*>::iterator it =  m_plots.begin();
                    it != m_plots.end();
                    it++)
                delete (*it);
        }
        void SetHits(std::vector<hit> hits)
        {
            for( std::vector<hit>::iterator it=hits.begin(); it!=hits.end(); it++)
            {
                if(goodhit(*it))
                {
                    m_plots[(int)it->pos.Z()]->Fill(it->pos.X(),it->pos.Y(),it->energy);
                }
            }
        }

        void FitEvent(float angle[2], bool print)
        {
            TGraph * xgraph;
            TGraph * ygraph;

            int nPoints = -1;
            float xavg[N_LAYERS];
            float yavg[N_LAYERS];
            float xerr[N_LAYERS];
            float yerr[N_LAYERS];
            float layers[N_LAYERS];
            for (int iL(0); iL<N_LAYERS; ++iL)
            {
                if(m_plots[iL]->GetEntries() > 0)
                {
                    nPoints++;
                    //std::cout << "supposedly valid: " << layers[iL] << ":" << xavg[iL] << ":" << xerr[iL] << std::endl;;
                }
                else continue;

                layers[nPoints]=m_dim.getZ(iL);
                xavg[nPoints]=m_plots[iL]->GetMean(1);
                yavg[nPoints]=m_plots[iL]->GetMean(2);

                xerr[nPoints]=1.0f/m_plots[iL]->Integral();
                yerr[nPoints]=1.0f/m_plots[iL]->Integral();

                //std::cout << sumw2[i][iL]  << ":"
                //          << p_xy[i][iL]->GetSumOfWeights() << ":"
                //          << p_xy[i][iL]->Integral() << ":"
                //          << p_xy[i][iL]->GetEntries() << ":"
                //          << std::endl;

                //Error for each hit is based solely on cell size
                //xerr[nPoints]=200.0f/(80.0f*TMath::Sqrt(12))*TMath::Sqrt(sumw2[i][iL]/TMath::Power(p_xy[i][iL]->Integral(),2));
                //yerr[nPoints]=xerr[nPoints];

            }
            float * layererr = 0;

            if(nPoints >5)
            {
                TCanvas * c = new TCanvas("asalkfj","asfka");

                xgraph = new TGraphErrors(nPoints,layers,xavg,layererr,xerr);
                xgraph->Draw("AP");
                xgraph->Fit("pol1","Q");

                TF1 * fit = xgraph->GetFunction("pol1");
                float xangle =  fit->Derivative(m_dim.getZ(0));

                std::ostringstream title;
                title << "angle: " << xangle;
                xgraph->SetTitle(title.str().c_str());

                if(print)
                {
                    std::ostringstream saveName;
                    saveName << "PLOTS/Fit_" << m_event << "_" << m_genEn << "GeV" << ".png";
                    c->SaveAs(saveName.str().c_str());
                }
                delete c;



                ygraph = new TGraphErrors(nPoints,layers,yavg,layererr,yerr);
                ygraph->Fit("pol1","Q");
                fit = ygraph->GetFunction("pol1");
                float yangle =  fit->Derivative(m_dim.getZ(0));
                
                angle[0]=xangle;
                angle[1]=yangle;
            }
        }

        void WritePlots()
        {
            TCanvas *myc = new TCanvas("myc","myc",1);
            std::ostringstream saveName;

            for(int i = 0; i<N_LAYERS; i++)
            {
                m_plots[i]->Draw("colz");
                Position line(0,0,i);
                //printPos(line,"ellipse");
                line = layer2z(line);
                //printPos(line,"ellipse");
                line+=m_max;
                line = m_slope*line.Z();

                //printPos(line,"ellipse");
                TEllipse el(line.X(),line.Y(),m_radius,m_radius);
                el.SetFillStyle(0);
                el.Draw();

                myc->Update();
                saveName.str("");
                saveName << "PLOTS/cone_"  << m_event << "_";
                if (i<10) saveName << "0";
                saveName << i << "_" << m_genEn << "GeV.png";
                myc->Print(saveName.str().c_str());
            }
            delete myc;
        }

    private:
        bool goodhit(hit thishit)
        {
            Position pos =layer2z(thishit.pos)-m_max;
            //printPos(pos,"pos");
            //printPos(m_slope*pos.Z(),"line");
            Position dist = pos - m_slope*pos.Z();
            //printPos(dist,"dist");
            return (dist.R() < m_radius);
        }
        void printPos(Position pos,std::string name)
        {
            std::cout << name << ": " << pos.X() << " " << pos.Y() << " " << pos.Z() << std::endl;
        }
        Position layer2z(Position pos)
        {
            return pos.SetZ(m_dim.getZ(pos.Z()));
        }

        Position m_max;
        Position m_vertex;
        Position m_slope;
        CaloDimen m_dim;
        std::string m_version;
        std::vector<TH2F*> m_plots;
        int m_event;
        int m_genEn;
        float m_radius;
}; //class EventShowerCone

int main(int argc, char** argv){//main  

    if (argc < 5) {
        std::cout << " Usage: " 
            << argv[0] << " <nEvts to process (0=all)>"
            << " <radius>" 
            << " <input_dir>" 
            << " <angle>" 
            << " <optional: debug (default=0)>"
            << std::endl;
        return 1;
    }

    const unsigned pNevts = atoi(argv[1]);
    const float radius = atof(argv[2]);
    const std::string inDir = argv[3];
    int start = inDir.find("version");
    const std::string version = inDir.substr(start,inDir.find("/",start)-start);
    bool debug = false;
    std::string angle(argv[4]);
    if (argc >5) debug = atoi(argv[5]);

    TFile *outputFile = TFile::Open("vectorfit.root","RECREATE");

    const unsigned nLayers = N_LAYERS;

    //unsigned genEn[]={5,10,25,50,75,100,150};
    unsigned genEn[]={10,50,100};
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
        input << inDir << "e_" << genEn[iE] << "_" << angle << "/PFcal.root";
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

        TH2F *p_xy[nEvts][nLayers];
        TH3F *p_xyz[nEvts];
        std::vector<hit> hits[nEvts];

        float sumw2[nEvts][nLayers];
        for (unsigned i(0); i<nEvts; ++i){
            for (unsigned iL(0); iL<nLayers; ++iL){
                lName << "p_xy_" << i << "_" << iL;
                p_xy[i][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100);
                lName.str("");

                sumw2[i][iL] = 0;
            }
            lName << "p_xyz_" << i;
            p_xyz[i] = new TH3F(lName.str().c_str(),";x(mm);y(mm)",80,-100,100,80,-100,100,nLayers,0,nLayers);
            lName.str("");

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
                    p_xyz[int(event)]->Fill(posx,posy,layer,weightedE);
                    sumw2[int(event)][layer]+=TMath::Power(weightedE,2);
                    hits[int(event)].push_back(hit(Position(posx,posy,layer),weightedE));

                }
            }//loop on hits

        }//loop on entries

        for (unsigned i(0); i<nEvts; ++i)
        {
            if (i%10 == 0) std::cout << "... Fitting event: " << i << std::endl;

            float showermax_layer = p_xyz[i]->ProjectionZ()->GetMaximumBin();
            Position showermax(p_xy[i][(int)showermax_layer]->GetMean(1),
                    p_xy[i][(int)showermax_layer]->GetMean(2),
                    showermax_layer);
            Position vertex(0.0f,0.0f,-143.55f);
            EventShowerCone cone(showermax,vertex,version,i,genEn[iE], radius);
            cone.SetHits(hits[i]);
            if(i == 5) cone.WritePlots();
            float angle[2];
            cone.FitEvent(angle, (event == 5) );
            p_angle[iE]->Fill(angle[0],angle[1]);
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
