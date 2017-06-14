#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TEfficiency.h"

const int nPtBin = 6;
double ptBins[nPtBin+1]  = {10.,20.,35.,50.,90.,150.,500.};
const int nEtaBin= 10;
double etaBins[nEtaBin+1] = {-2.5,-2.0,-1.566,-1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};

void countTrigEff4WP(const char wp[]="Loose"){
    // Not workging fine, why?
    TCanvas *canvas = new TCanvas();

    TH2D *h_passed = new TH2D("h_passed","", nPtBin, ptBins, nEtaBin, etaBins);
    TH2D *h_total  = new TH2D("h_total" ,"", nPtBin, ptBins, nEtaBin, etaBins);

    TChain *ch = new TChain("ntupler/ntuple");
    ch->Add("./data_Run2016*-03Feb2017*.root");

    TString passCut = TString::Format("tag%s > 0 && probe%s > 0", wp, wp);
    TString totalCut= TString::Format("tag%s > 0", wp);
    std::cout << ch->Draw("probePT:probeEta>>h_passed",passCut .Data(),"goff") << std::endl;
    std::cout << ch->Draw("probePT:probeEta>>h_total" ,totalCut.Data(),"goff") << std::endl;

    h_total->Draw("COLZ");
    canvas->Update();
    canvas->Print(TString::Format("tnp_nTotal_%s.pdf",wp));

    h_passed->Draw("COLZ");
    canvas->Update();
    canvas->Print(TString::Format("tnp_nPassed_%s.pdf",wp));

    TEfficiency *tEf = new TEfficiency(*h_passed,*h_total);

    tEf->Draw("COLZ");
    canvas->Update();
    canvas->Print(TString::Format("tnp_trigEff_%s.pdf",wp));

    // TODO: Print tables

    delete tEf;
    delete ch;
    delete h_total;
    delete h_passed;
    return;
}

void countTrigEff4WPFull(const char wp[]="Loose"){
    TCanvas *canvas = new TCanvas();

    TH2D *h_passed = new TH2D("h_passed","", nPtBin, ptBins, nEtaBin, etaBins);
    TH2D *h_total  = new TH2D("h_total" ,"", nPtBin, ptBins, nEtaBin, etaBins);

    TChain *ch = new TChain("ntupler/ntuple");
    ch->Add("./data_Run2016*-03Feb2017*.root");
    double tagWP, probeWP;
    double probePT, probeEta;

    ch->SetBranchAddress(TString::Format("tag%s"  ,wp).Data(), &tagWP  );
    ch->SetBranchAddress(TString::Format("probe%s",wp).Data(), &probeWP);
    ch->SetBranchAddress("probePT",&probePT);
    ch->SetBranchAddress("probeEta",&probeEta);

    //printf("DEBUT\t: nEntries=%lld\n",ch->GetEntries());
    for(int entry=1; entry<=ch->GetEntries(); entry++){
        ch->GetEntry(entry);
        //printf("DEBUG\t: tagWP=%f\n",tagWP);
        //printf("DEBUG\t: probeWP=%f\n",probeWP);
        //printf("DEBUG\t: probePT=%f\n",probePT);
        //printf("DEBUG\t: probeEta=%f\n",probeEta);
        if (tagWP > 0){
            h_total->Fill(probePT,probeEta);
            if (probeWP > 0){
                h_passed->Fill(probePT,probeEta);
            }
        }
    }

    h_total->Draw("COLZ");
    canvas->Update();
    canvas->Print(TString::Format("tnp_nTotal_%s.pdf",wp));

    h_passed->Draw("COLZ");
    canvas->Update();
    canvas->Print(TString::Format("tnp_nPassed_%s.pdf",wp));

    TEfficiency *tEf = new TEfficiency(*h_passed,*h_total);
    tEf->Draw("COLZ");
    canvas->Update();
    auto h_tEf = tEf->GetPaintedHistogram(); 
    h_tEf->GetZaxis()->SetRangeUser(0.5,1.);

    tEf->Draw("COLZ");
    canvas->Update();
    canvas->Print(TString::Format("tnp_trigEff_%s.pdf",wp));

    // TODO: Print tables

    delete tEf;
    delete ch;
    delete h_total;
    delete h_passed;
    return;
}

void countTrigEff(){
    countTrigEff4WPFull("Loose");   
    countTrigEff4WPFull("CB");   
    countTrigEff4WPFull("CBTriChg");   
    countTrigEff4WPFull("MVA90Iso16");   
    countTrigEff4WPFull("MVA90Iso15");   
    countTrigEff4WPFull("MVA80Iso16");   
    countTrigEff4WPFull("MVA80Iso15");   

    return;
}

int main(){
    countTrigEff();
    return 0;
}
