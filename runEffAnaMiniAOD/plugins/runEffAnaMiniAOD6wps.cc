// -*- C++ -*-
//
// Package:    EleEfficiencyHWW/runEffAnaMiniAOD
// Class:      runEffAnaMiniAOD
// 
/**\class runEffAnaMiniAOD runEffAnaMiniAOD.cc EleEfficiencyHWW/runEffAnaMiniAOD/plugins/runEffAnaMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Po-Hsun Chen
//         Created:  Wed, 15 Feb 2017 08:20:27 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Load here all the dataformat that we will need
//#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Common/interface/ValueMap.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

//#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
//#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/LumiUtils.h"

//L1 EM particles
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>
#include <TMath.h>

#include "MyPatElectronMVAIdSelector.cc"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

const int nLegs=4;
std::string ControlFilterLeg = "hltEle27WPTightGsfTrackIsoFilter";
std::string legsFilter[nLegs] = {
    // DoubleEle Trigger
    "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",  // Leg 1 of HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2
    "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",  // Leg 2 of HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2
    //Singleectr Trigger
    "hltEle27WPTightGsfTrackIsoFilter",     //HLT_Ele27_WPTight_Gsf_v*
    "hltEle25erWPTightGsfTrackIsoFilter"  //HLT_Ele25_eta2p1_WPTight_Gsf_v*
};
class runEffAnaMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>
{//{{{
   public:
      explicit runEffAnaMiniAOD(const edm::ParameterSet&);
      ~runEffAnaMiniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::Service<TFileService> fs;

      // ----------member data ---------------------------
      void      reset();
      bool      passFilter(pat::TriggerObjectStandAlone&, std::string&);
      double    DeltaRtrig(const pat::Electron& , std::vector < pat::TriggerObjectStandAlone > );
      bool      hasZMother(const reco::GenParticle);

      MyPatElectronMVAIdSelector electronMVATightIdExtra;

      // Tokens
          // Data members that are the same for AOD and miniAOD
      edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>triggerPrescale_; 
      //edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<double> rhoToken_;     
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
                                              
      // MiniAOD case data members            
      edm::EDGetToken muonsMiniAODToken_;     
      edm::EDGetToken jetToken_;              
      edm::EDGetToken electronsMiniAODToken_; 
      edm::EDGetToken photonsMiniAODToken_;   
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;
                                              
      // ID decisions objects                 
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > elenontrigMVAlooseMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > elenontrigMVAtightMapToken_;
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genToken_;

      // Tags to be used in data processing
      bool    isMC;
      double  xsec;
      TString dtag;

      // Tree and branches
      TTree *ntuple;
      double Nvtx;
      double ZMass;
      double tagPT;
      double tagEta;
      double tagPhi;
      double tagMVA80Iso15;
      double tagMVA80Iso16;
      double tagMVA90Iso15;
      double tagMVA90Iso16;
      double tagCB;
      double tagCBTriChg;
      double probeMVA80Iso15;
      double probeMVA80Iso16;
      double probeMVA90Iso15;
      double probeMVA90Iso16;
      double probeCB;  
      double probeCBTriChg;  
      double probePT;
      double probeEta;
      double probePhi;
      //double nPU;
      //double nPUTrue;
      double PUWeight;
      double lumi_weight;
      double passProbePTLegs [nLegs];
      double passProbeEtaLegs[nLegs];
      double passProbePhiLegs[nLegs];
      double counting;
      double count_tag;
      double count_probe;
      double count_passprobeleg1;
      double count_passprobeleg2;
      double count_passprobeleg3;
      double passL1PTNorm;
      double passL1PTIso;
      double passL1PTER;
      double passL1PT;
      double passL1EtaNorm;
      double passL1EtaIso;
      double passL1EtaER;
      double passL1Eta;

      TBranch *b_Nvtx                ;
      TBranch *b_ZMass               ;
      TBranch *b_tagMVA80Iso15       ;
      TBranch *b_tagMVA80Iso16       ;
      TBranch *b_tagMVA90Iso15       ;
      TBranch *b_tagMVA90Iso16       ;
      TBranch *b_tagCB               ;
      TBranch *b_tagCBTriChg             ;
      TBranch *b_probeMVA80Iso15       ;
      TBranch *b_probeMVA80Iso16       ;
      TBranch *b_probeMVA90Iso15       ;
      TBranch *b_probeMVA90Iso16       ;
      TBranch *b_probeCB             ;
      TBranch *b_probeCBTriChg             ;
      TBranch *b_tagPT               ;
      TBranch *b_tagEta              ;
      TBranch *b_tagPhi              ;
      TBranch *b_probePT        ;
      TBranch *b_probeEta       ;
      TBranch *b_probePhi       ;
      TBranch *b_passProbePTLegs [nLegs];
      TBranch *b_passProbeEtaLegs[nLegs];
      TBranch *b_passProbePhiLegs[nLegs];
      //TBranch *b_nPU            ;
      //TBranch *b_nPUTrue        ;
      TBranch *b_PUWeight            ;
      TBranch *b_lumi_weight         ;
      TBranch *b_counting            ;
      TBranch *b_count_tag           ;
      TBranch *b_count_probe         ;
      TBranch *b_count_paasprobeleg1 ;
      TBranch *b_count_paasprobeleg2 ;
      TBranch *b_passL1PTNorm        ;
      TBranch *b_passL1PTIso         ;
      TBranch *b_passL1PTER          ;
      TBranch *b_passL1PT            ;
      TBranch *b_passL1EtaNorm       ;
      TBranch *b_passL1EtaIso        ;
      TBranch *b_passL1EtaER         ;
      TBranch *b_passL1Eta           ;
};//}}}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
runEffAnaMiniAOD::runEffAnaMiniAOD(const edm::ParameterSet & iConfig)
{                               //{{{
    eleVetoIdMapToken_           = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("eleVetoIdMap"));
    eleLooseIdMapToken_          = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("eleLooseIdMap"));
    eleMediumIdMapToken_         = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("eleMediumIdMap"));
    eleTightIdMapToken_          = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("eleTightIdMap"));
    eleHEEPIdMapToken_           = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("eleHEEPIdMap"));
    elenontrigMVAlooseMapToken_  = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("elenontrigMVAlooseIdMap"));
    elenontrigMVAtightMapToken_  = consumes < edm::ValueMap < bool > >(iConfig.getParameter < edm::InputTag > ("elenontrigMVAtightIdMap"));
    mvaValuesMapToken_           = consumes < edm::ValueMap < float > >(iConfig.getParameter < edm::InputTag > ("mvaValuesMap"));
    mvaCategoriesMapToken_       = consumes < edm::ValueMap <int> >(iConfig.getParameter < edm::InputTag > ("mvaCategoriesMap"));

    triggerToken_                = mayConsume < edm::TriggerResults > (iConfig.getParameter < edm::InputTag > ("trigger"));
    triggerObjects_              = consumes < pat::TriggerObjectStandAloneCollection >(iConfig.getParameter < edm::InputTag > ("objects"));
    triggerPrescale_             = consumes < pat::PackedTriggerPrescales > (iConfig.getParameter < edm::InputTag > ("prescale"));
    //pileupToken_                 = consumes < edm::View < PileupSummaryInfo > >(iConfig.getParameter < edm::InputTag > ("pileup"));
    rhoToken_                    = consumes < double > (iConfig.getParameter < edm::InputTag > ("rho"));
    beamSpotToken_               = consumes < reco::BeamSpot > (iConfig.getParameter < edm::InputTag > ("beamSpot"));

    muonsMiniAODToken_           = mayConsume < edm::View < pat::Muon > >(iConfig.getParameter < edm::InputTag > ("muonsMiniAOD"));
    electronsMiniAODToken_       = mayConsume < edm::View < pat::Electron > >(iConfig.getParameter < edm::InputTag > ("electronsMiniAOD"));
    conversionsMiniAODToken_     = mayConsume < reco::ConversionCollection >(iConfig.getParameter < edm::InputTag > ("conversionsMiniAOD"));
    vtxMiniAODToken_             = mayConsume < reco::VertexCollection >(iConfig.getParameter < edm::InputTag > ("verticesMiniAOD"));

    //now do what ever initialization is needed
    usesResource("TFileService");
    isMC = iConfig.getUntrackedParameter < bool > ("isMC");
    xsec = iConfig.getUntrackedParameter < double >("xsec");
    dtag = iConfig.getUntrackedParameter < std::string > ("dtag");

}//}}}

runEffAnaMiniAOD::~runEffAnaMiniAOD()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
runEffAnaMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   reset();

   bool passOnlyEle(false);
   double rho;

   //load all the objects we will need to access
   edm::Handle < edm::TriggerResults > triggerBits;
   iEvent.getByToken(triggerToken_,triggerBits);
   const edm::TriggerNames & names = iEvent.triggerNames(*triggerBits);


   edm::Handle < edm::View< pat::Electron > >electronsHandle;
   iEvent.getByToken(electronsMiniAODToken_, electronsHandle);

   reco::VertexCollection vtx;
   edm::Handle < reco::VertexCollection > vtxHandle;
   iEvent.getByToken(vtxMiniAODToken_,vtxHandle);
   if (vtxHandle.isValid()) {
       vtx = *vtxHandle;
   }

   int nPV = vtxHandle->size();
   int firstGoodVertexIdx = 0;

   reco::VertexCollection::const_iterator firstGoodVertex = vtxHandle->end();
   for (reco::VertexCollection::const_iterator vtxIter = vtxHandle->begin(); vtxIter != vtxHandle->end();
        ++vtxIter, ++firstGoodVertexIdx) {
       firstGoodVertex = vtxIter;
       break;
   }
   electronMVATightIdExtra.setFGV(firstGoodVertex);

   edm::Handle < edm::View<pat::Muon> > muonsHandle;
   iEvent.getByToken(muonsMiniAODToken_,muonsHandle);

   edm::Handle < reco::BeamSpot > theBeamSpotH;
   iEvent.getByToken(beamSpotToken_,theBeamSpotH);
   electronMVATightIdExtra.setBS(theBeamSpotH);

   //edm::Handle <edm::View<PileupSummaryInfo> > pileupHandle;
   //iEvent.getByToken(pileupToken_,pileupHandle);
   //for( auto & puInfoElement : *pileupHandle){
   //   if( puInfoElement.getBunchCrossing() == 0 ){
   //      nPU     = puInfoElement.getPU_NumInteractions();
   //      nPUTrue = puInfoElement.getTrueNumInteractions();
   //   }
   //}

   edm::Handle < double >rhoH;
   iEvent.getByToken(rhoToken_,rhoH);
   rho = *rhoH;
   electronMVATightIdExtra.setRho(rho);

   edm::Handle < reco::ConversionCollection > conversionsH;
   iEvent.getByToken(conversionsMiniAODToken_, conversionsH);
   electronMVATightIdExtra.setConversion(conversionsH);

     // Get the electron ID data from the event stream.
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
   edm::Handle<edm::ValueMap<bool> > mva_loose_id_decisions;
   edm::Handle<edm::ValueMap<bool> > mva_tight_id_decisions;
                                          
   iEvent.getByToken(elenontrigMVAlooseMapToken_,mva_loose_id_decisions);
   iEvent.getByToken(elenontrigMVAtightMapToken_,mva_tight_id_decisions);
   iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
   iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
   iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
   iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);
   iEvent.getByToken(eleHEEPIdMapToken_ ,heep_id_decisions);
                                          
   // Get MVA values and categories       
   edm::Handle<edm::ValueMap<float> > mvaValues;
   edm::Handle<edm::ValueMap<int> > mvaCategories;
   iEvent.getByToken(mvaValuesMapToken_,mvaValues);
   iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

   //edm::Handle < l1extra::L1EmParticleCollection > theL1extraParticles;
   //iEvent.getByLabel("l1extraParticles:NonIsolated",theL1extraParticles);

   //edm::Handle < l1extra::L1EmParticleCollection > theL1extraParticlesIsolated;
   //iEvent.getByLabel("l1extraParticles:Isolated",theL1extraParticlesIsolated);

   edm::Handle < pat::TriggerObjectStandAloneCollection > triggerObjects;
   iEvent.getByLabel("selectedPatTrigger",triggerObjects);

   std::vector < pat::TriggerObjectStandAlone > TagtriggerObj;
   std::vector < pat::TriggerObjectStandAlone > ProbetriggerObj;
   std::vector < pat::TriggerObjectStandAlone > passLegsObj[nLegs];

   ///////////////////////
   ///                 ///
   /// LEPTON ANALYSIS ///
   ///                 ///
   ///////////////////////
   std::vector < patUtils::GenericLepton > leptons;
   for (size_t l = 0; l < electronsHandle->size(); l++) {
       leptons.push_back(patUtils::GenericLepton((*(electronsHandle->ptrAt(l)))));
   }
   for (size_t l = 0; l < muonsHandle->size(); l++) {
       leptons.push_back(patUtils::GenericLepton(*(muonsHandle->ptrAt(l))));
   }
   std::sort(leptons.begin(), leptons.end(), utils::sort_CandidatesByPt);

   std::vector < patUtils::GenericLepton > selLeptons;
   //Request exactly two leptons

   //PRE-SELECTION based to pt value
   for (unsigned int j = 0; j < leptons.size(); j++) {
       if (leptons[j].pt() > 8)
           selLeptons.push_back(leptons[j]);
   }

   std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);
   if (selLeptons.size() != 2)
       return;

   // opposite charge condition
   if (!(selLeptons[0].charge() * selLeptons[1].charge() < 0)){
       return;
   }

   if (abs(selLeptons[0].pdgId()) == 11 && abs(selLeptons[1].pdgId()) == 11){
       passOnlyEle = true;
   }

    /// ELE TAG ///
    //Logical tag for Tag and Probe Ele
    bool passTagEleKin(false)    , passProbeEleKin(false);
    bool TagTightIdSTD(false)    , ProbeTightIdSTD(false);
    bool TagTightIdTriChgSTD(false)    , ProbeTightIdTriChgSTD(false);
    bool Tag80XMVAwp80Iso15(false) , Probe80XMVAwp80Iso15(false);
    bool Tag80XMVAwp80Iso16(false) , Probe80XMVAwp80Iso16(false);
    bool Tag80XMVAwp90Iso15(false) , Probe80XMVAwp90Iso15(false);
    bool Tag80XMVAwp90Iso16(false) , Probe80XMVAwp90Iso16(false);
    bool TagEle(false)           , ProbeEle(false);
    bool ZPick(false);

    bool passControlLeg1(false), passControlLeg2(false);
    bool passProbeLegs[nLegs];
    for (int i = 0; i < nLegs; i++) {
        passProbeLegs[i]=false;
    }

    if (passOnlyEle) {//{{{ ELECTRON EFFICENCY
        int first  = rand() % 2;
        int second = (first + 1) % 2;

        Nvtx = vtx.size();
        double ptf        = selLeptons[first ].pt();
        double pts        = selLeptons[second].pt();
        double etaf       = selLeptons[first ].el.superCluster()->eta();
        double etas       = selLeptons[second].el.superCluster()->eta();
        double phif       = selLeptons[first ].phi();
        double phis       = selLeptons[second].phi();
        double etafSC     = selLeptons[first ].el.superCluster()->eta();
        double etasSC     = selLeptons[second].el.superCluster()->eta();

        double ecalPFIsof = selLeptons[first ].el.ecalPFClusterIso();
        double ecalPFIsos = selLeptons[second].el.ecalPFClusterIso();
        double hcalPFIsof = selLeptons[first ].el.hcalPFClusterIso();
        double hcalPFIsos = selLeptons[second].el.hcalPFClusterIso();
        double trackIsof  = selLeptons[first ].el.trackIso();
        double trackIsos  = selLeptons[second].el.trackIso();
        double detaseedf  = selLeptons[first ].el.deltaEtaSeedClusterTrackAtCalo();
        double detaseeds  = selLeptons[second].el.deltaEtaSeedClusterTrackAtCalo();

        double d0f        = selLeptons[first ].el.gsfTrack()->dxy(vtx[firstGoodVertexIdx].position());
        double d0s        = selLeptons[second].el.gsfTrack()->dxy(vtx[firstGoodVertexIdx].position());
        double dzf        = selLeptons[first ].el.gsfTrack()->dz (vtx[firstGoodVertexIdx].position());
        double dzs        = selLeptons[second].el.gsfTrack()->dz (vtx[firstGoodVertexIdx].position());
        
        double misshit    = selLeptons[second].el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        double dphiin     = selLeptons[second].el.deltaPhiSuperClusterTrackAtVtx();
        double chi2       = selLeptons[second].el.gsfTrack()->normalizedChi2();

        //cout << "detain seed = " << selLeptons[second].el.deltaEtaSeedClusterTrackAtCalo() << endl;
        //cout <<" mising hits = " << selLeptons[second].el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) << endl;

        //if (abs(etafSC) <= 1.479) {
        //    passTagISO = abs(d0f) < 0.05 && abs(dzf) < 0.10;
        //} else if (abs(etafSC) > 1.479 && abs(etafSC) < 2.5) {
        //    passTagISO = abs(d0f) < 0.10 && abs(dzf) < 0.20;
        //}

        //if (abs(etasSC) <= 1.479) {
            //passProbeISO = ((ecalPFIsos-rho*0.165)/pts) < 0.160 && ((hcalPFIsos-rho*0.060)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(detaseeds) < 0.004 && misshit < 1 && abs(dphiin) < 0.020;
            //passProbeISO = ((ecalPFIsos-rho*0.165)/pts) < 0.160 && ((hcalPFIsos-rho*0.060)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(detaseeds) < 0.004 && abs(dphiin) < 0.020;        
            //passProbeISO = abs(d0s) < 0.05 && abs(dzs) < 0.10;

        //} else if (abs(etasSC) > 1.479 && abs(etasSC) < 2.5) {
            //passProbeISO = ((ecalPFIsos-rho*0.132)/pts) < 0.120 && ((hcalPFIsos-rho*0.131)/pts) < 0.120 && (trackIsos/pts) < 0.08 && misshit < 1 && abs(chi2) < 3;
            //passProbeISO = ((ecalPFIsos-rho*0.132)/pts) < 0.120 && ((hcalPFIsos-rho*0.131)/pts) < 0.120 && (trackIsos/pts) < 0.08 && abs(chi2) < 3;
            //passProbeISO = abs(d0s) < 0.10 && abs(dzs) < 0.20;
        //}

//       passTagISO = ecalPFIsof < 0.45 && hcalPFIsof < 0.25 && trackIsof < 0.2;
//       passProbeISO = ecalPFIsos < 0.45 && hcalPFIsos < 0.25 && trackIsos < 0.2;

        for (pat::TriggerObjectStandAlone obj:*triggerObjects) {
            obj.unpackPathNames(names);

            // ======================= Triggers for Tag ================================  
            passControlLeg1 = passFilter(obj, ControlFilterLeg);
            if (passControlLeg1) {
                TagtriggerObj.push_back(obj);
            }

            for(int iLeg=0; iLeg<nLegs; iLeg++){
                passProbeLegs[iLeg] = passFilter(obj, legsFilter[iLeg]);
                if (passProbeLegs[iLeg]) {
                    passLegsObj[iLeg].push_back(obj);
                }
            }
        }

        // FIX ME!
        passTagEleKin = (((abs(etaf) >= 0 && abs(etaf) <= 1.4442) || (abs(etaf) >= 1.5660 && abs(etaf) <= 2.5)) && ptf > 35);

        TagTightIdSTD   = (*tight_id_decisions)[electronsHandle->ptrAt(first)];
        ProbeTightIdSTD = (*tight_id_decisions)[electronsHandle->ptrAt(second)];
        TagTightIdTriChgSTD   = (*tight_id_decisions)[electronsHandle->ptrAt(first)] && electronMVATightIdExtra.passTriChg(selLeptons[first].el);
        ProbeTightIdTriChgSTD = (*tight_id_decisions)[electronsHandle->ptrAt(second)] && electronMVATightIdExtra.passTriChg(selLeptons[second].el);
        Tag80XMVAwp80Iso15    = (*mva_tight_id_decisions)[electronsHandle->ptrAt(first)]   && electronMVATightIdExtra.passId(selLeptons[first].el ) && electronMVATightIdExtra.passIso(selLeptons[first].el ,15);
        Probe80XMVAwp80Iso15  = (*mva_tight_id_decisions)[electronsHandle->ptrAt(second)]  && electronMVATightIdExtra.passId(selLeptons[second].el) && electronMVATightIdExtra.passIso(selLeptons[second].el,15);
        Tag80XMVAwp80Iso16    = (*mva_loose_id_decisions)[electronsHandle->ptrAt(first)]   && electronMVATightIdExtra.passId(selLeptons[first].el ) && electronMVATightIdExtra.passIso(selLeptons[first].el ,15);
        Probe80XMVAwp80Iso16  = (*mva_loose_id_decisions)[electronsHandle->ptrAt(second)]  && electronMVATightIdExtra.passId(selLeptons[second].el) && electronMVATightIdExtra.passIso(selLeptons[second].el,15);
        Tag80XMVAwp90Iso15    = (*mva_tight_id_decisions)[electronsHandle->ptrAt(first)]   && electronMVATightIdExtra.passId(selLeptons[first].el ) && electronMVATightIdExtra.passIso(selLeptons[first].el ,16);
        Probe80XMVAwp90Iso15  = (*mva_tight_id_decisions)[electronsHandle->ptrAt(second)]  && electronMVATightIdExtra.passId(selLeptons[second].el) && electronMVATightIdExtra.passIso(selLeptons[second].el,16);
        Tag80XMVAwp90Iso16    = (*mva_loose_id_decisions)[electronsHandle->ptrAt(first)]   && electronMVATightIdExtra.passId(selLeptons[first].el ) && electronMVATightIdExtra.passIso(selLeptons[first].el ,16);
        Probe80XMVAwp90Iso16  = (*mva_loose_id_decisions)[electronsHandle->ptrAt(second)]  && electronMVATightIdExtra.passId(selLeptons[second].el) && electronMVATightIdExtra.passIso(selLeptons[second].el,16);

        bool passTagdRCut(false);

        double dRTag = DeltaRtrig(selLeptons[first].el, TagtriggerObj);

        if (dRTag < 0.1) {
            passTagdRCut = true;
        }

        if (passTagEleKin && TagTightIdSTD && passTagdRCut){
            TagTightIdSTD = true;
        }
        if (passTagEleKin && TagTightIdTriChgSTD && passTagdRCut){
            TagTightIdTriChgSTD = true;
        }
        if (passTagEleKin && Tag80XMVAwp80Iso15 && passTagdRCut){
            Tag80XMVAwp80Iso15 = true;
        }
        if (passTagEleKin && Tag80XMVAwp80Iso16 && passTagdRCut){
            Tag80XMVAwp80Iso16 = true;
        }
        if (passTagEleKin && Tag80XMVAwp90Iso15 && passTagdRCut){
            Tag80XMVAwp90Iso15 = true;
        }
        if (passTagEleKin && Tag80XMVAwp90Iso16 && passTagdRCut){
            Tag80XMVAwp90Iso16 = true;
        }

        TagEle = Tag80XMVAwp80Iso15 || Tag80XMVAwp80Iso16 || Tag80XMVAwp90Iso15 || Tag80XMVAwp90Iso16 || TagTightIdSTD || TagTightIdTriChgSTD;
        if (!TagEle){
            return;
        }
        
        tagPT = ptf;
        tagEta = etaf;
        tagPhi = phif;
        tagMVA80Iso15 = Tag80XMVAwp80Iso15;
        tagMVA80Iso16 = Tag80XMVAwp80Iso16;
        tagMVA90Iso15 = Tag80XMVAwp90Iso15;
        tagMVA90Iso16 = Tag80XMVAwp90Iso16;
        tagCB  = TagTightIdSTD;
        tagCBTriChg = TagTightIdTriChgSTD;

        //Selection of the Probe
        bool ProbeEle = false;
        passProbeEleKin = ((abs(etas) >= 0 && abs(etas) <= 2.5) && pts > 0);    //Fix it       
        if ((ProbeTightIdSTD && passProbeEleKin )){
            ProbeTightIdSTD = true;
        }
        if ((ProbeTightIdTriChgSTD && passProbeEleKin)){
            ProbeTightIdTriChgSTD = true;
        }
        if ((Probe80XMVAwp80Iso15 && passProbeEleKin )){
            Probe80XMVAwp80Iso15 = true;
        }
        if ((Probe80XMVAwp80Iso16 && passProbeEleKin )){
            Probe80XMVAwp80Iso16 = true;
        }
        if ((Probe80XMVAwp90Iso15 && passProbeEleKin )){
            Probe80XMVAwp90Iso15 = true;
        }
        if ((Probe80XMVAwp90Iso16 && passProbeEleKin )){
            Probe80XMVAwp90Iso16 = true;
        }

        //ProbeEle = ProbeTightMVAIdSTD || ProbeTightIdSTD;
        ProbeEle = Probe80XMVAwp80Iso15 || Probe80XMVAwp80Iso16 || Probe80XMVAwp90Iso15 || Probe80XMVAwp90Iso16 || ProbeTightIdSTD || ProbeTightIdTriChgSTD;
        probeMVA80Iso15 = Probe80XMVAwp80Iso15;
        probeMVA80Iso16 = Probe80XMVAwp80Iso16;
        probeMVA90Iso15 = Probe80XMVAwp90Iso15;
        probeMVA90Iso16 = Probe80XMVAwp90Iso16;
        probeCB  = ProbeTightIdSTD;
        probeCBTriChg = ProbeTightIdTriChgSTD;

        TLorentzVector lep1(selLeptons[first].px(),
                            selLeptons[first].py(),
                            selLeptons[first].pz(),
                            selLeptons[first].energy());
        TLorentzVector lep2(selLeptons[second].px(),
                            selLeptons[second].py(),
                            selLeptons[second].pz(),
                            selLeptons[second].energy());
        ZMass = (lep1 + lep2).M();
        ZPick = ZMass > 60 && ZMass < 120;
        if (!ZPick){
            return;
        }

        if (ProbeEle) {//{{{

            probePT  = pts;
            probeEta = etas;
            probePhi = phis;

//now lon the L1 particles

            //l1extra::L1EmParticleCollection::const_iterator itrEm;
            //float maxL1matched = -1;
            //for (itrEm = theL1extraParticles->begin(); itrEm != theL1extraParticles->end(); ++itrEm) {
            //    float deltaR = sqrt(pow(itrEm->eta() - etas, 2) + pow(acos(cos(itrEm->phi() - phis)), 2));
            //    if (deltaR < 0.5) {
            //        if (itrEm->pt() > maxL1matched)
            //            maxL1matched = itrEm->pt();
            //    }
            //}

            //l1extra::L1EmParticleCollection::const_iterator itrEmIso;
            //for (itrEmIso = theL1extraParticlesIsolated->begin(); itrEmIso != theL1extraParticlesIsolated->end();
            //     ++itrEmIso) {
            //    float deltaR = sqrt(pow(itrEmIso->eta() - etas, 2) + pow(acos(cos(itrEmIso->phi() - phis)), 2));
            //    if (deltaR < 0.5) {
            //        if (itrEmIso->pt() > maxL1matched)
            //            maxL1matched = itrEmIso->pt();
            //    }
            //}

            for (int iLeg = 0; iLeg < nLegs; iLeg++) {
                if (DeltaRtrig(selLeptons[second].el, passLegsObj[iLeg]) < 0.1){
                    passProbePTLegs [iLeg] = pts;
                    passProbeEtaLegs[iLeg] = etas;
                    passProbePhiLegs[iLeg] = phis;
                }
            }

            // ======== L1 Trigger Efficiency
            float maxL1MatchedNorm = -1;
            float maxL1MatchedIso  = -1;
            float maxL1MatchedER   = -1;
            bool L1Norm(false), L1Iso(false), L1ER(false);
            //for (int ibx = dmxegs->getFirstBX(); ibx <= dmxegs->getLastBX(); ++ibx) {
            //    for (auto itr = dmxegs->begin(ibx); itr != dmxegs->end(ibx); ++itr) {
            //        float deltaR = sqrt(pow(itr->eta() - etas, 2) + pow(acos(cos(itr->phi() - phis)), 2));
            //        if (deltaR < 0.5) {
            //            //cout << "ET =  : " << itr->et() << " pt=" << itr->pt() << " eta=" << itr->eta() <<  " phi=" << itr->phi() << std::endl;
            //            //cout << "Iso = " << itr->hwIso() << endl; 
            //            if (itr->pt() > maxL1MatchedNorm)
            //                maxL1MatchedNorm = itr->pt();
            //            if (itr->hwIso() == 1 && itr->pt() > maxL1MatchedIso)
            //                maxL1MatchedIso = itr->pt();
            //            if (itr->hwIso() == 1 && fabs(itr->eta()) < 2.1 && itr->pt() > maxL1MatchedER)
            //                maxL1MatchedER = itr->pt();
            //        }
            //    }
            //}

            //if (maxL1MatchedNorm > 40){
            //    L1Norm = true;
            //    passL1PTNorm = pts;
            //    passL1EtaNorm = etas;
            //}
            //if (maxL1MatchedIso > 24){
            //    L1Iso = true;
            //    passL1PTIso = pts;
            //    passL1EtaIso = etas;
            //}
            //if (maxL1MatchedER > 22){
            //    L1ER = true;
            //    passL1PTER = pts;
            //    passL1EtaER = etas;
            //}
            //if (L1Norm == 1 || L1Iso == 1 || L1ER == 1) {
            //    passL1PT = pts;
            //    passL1Eta = etas;
            //}
        } //}}}If you find Probe Electron
        ntuple->Fill();
    }//}}}passOnlyEle
}


// ------------ method called once each job just before starting event loop  ------------
void 
runEffAnaMiniAOD::beginJob()
{//{{{
ntuple = fs->make<TTree>("ntuple", "Efficiency Tree");

b_Nvtx                = ntuple->Branch("Nvtx"              , &Nvtx              , "Nvtx/D");
b_ZMass               = ntuple->Branch("ZMass"             , &ZMass             , "ZMass/D");
b_tagMVA80Iso15              = ntuple->Branch("tagMVA80Iso15"            ,&tagMVA80Iso15             , "tagMVA80Iso15/D");
b_tagMVA80Iso16              = ntuple->Branch("tagMVA80Iso16"            ,&tagMVA80Iso16             , "tagMVA80Iso16/D");
b_tagMVA90Iso15              = ntuple->Branch("tagMVA90Iso15"            ,&tagMVA90Iso15             , "tagMVA90Iso15/D");
b_tagMVA90Iso16              = ntuple->Branch("tagMVA90Iso16"            ,&tagMVA90Iso16             , "tagMVA90Iso16/D");
b_tagCB               = ntuple->Branch("tagCB"             ,&tagCB              , "tagCB/D");
b_tagCBTriChg               = ntuple->Branch("tagCBTriChg"             ,&tagCBTriChg              , "tagCBTriChg/D");
b_probeMVA80Iso15            = ntuple->Branch("probeMVA80Iso15"          ,&probeMVA80Iso15           , "probeMVA80Iso15/D");
b_probeMVA80Iso16            = ntuple->Branch("probeMVA80Iso16"          ,&probeMVA80Iso16           , "probeMVA80Iso16/D");
b_probeMVA90Iso15            = ntuple->Branch("probeMVA90Iso15"          ,&probeMVA90Iso15           , "probeMVA90Iso15/D");
b_probeMVA90Iso16            = ntuple->Branch("probeMVA90Iso16"          ,&probeMVA90Iso16           , "probeMVA90Iso16/D");
b_probeCB             = ntuple->Branch("probeCB"           ,&probeCB            , "probeCB/D");
b_probeCBTriChg             = ntuple->Branch("probeCBTriChg"           ,&probeCBTriChg            , "probeCBTriChg/D");
b_tagPT               = ntuple->Branch("tagPT"             , &tagPT             , "tagPT/D");
b_tagEta              = ntuple->Branch("tagEta"            , &tagEta            , "tagEta/D");
b_tagPhi              = ntuple->Branch("tagPhi"            , &tagPhi            , "tagPhi/D");
b_probePT             = ntuple->Branch("probePT"      , &probePT      , "probePT/D");
b_probeEta            = ntuple->Branch("probeEta"     , &probeEta     , "probeEta/D");
b_probePhi            = ntuple->Branch("probePhi"     , &probePhi     , "probePhi/D");
//b_nPU                 = ntuple->Branch("nPU"            , &nPU            , "nPU/D");
//b_nPUTrue             = ntuple->Branch("nPUTrue"            , &nPUTrue            , "nPUTrue/D");
b_PUWeight            = ntuple->Branch("PUWeight"            , &PUWeight            , "PUWeight/D");
b_lumi_weight         = ntuple->Branch("lumi_weight"         , &lumi_weight         , "lumi_weight/D");

for(int iLeg=0; iLeg<nLegs; iLeg++){
    b_passProbePTLegs [iLeg] = ntuple->Branch(TString::Format("passProbePTLeg%02d" ,iLeg).Data(), &passProbePTLegs [iLeg],   TString::Format("passProbePTLeg%02d/D" ,iLeg).Data());
    b_passProbeEtaLegs[iLeg] = ntuple->Branch(TString::Format("passProbeEtaLeg%02d",iLeg).Data(), &passProbeEtaLegs[iLeg],   TString::Format("passProbeEtaLeg%02d/D",iLeg).Data());
    b_passProbePhiLegs[iLeg] = ntuple->Branch(TString::Format("passProbePhiLeg%02d",iLeg).Data(), &passProbePhiLegs[iLeg],   TString::Format("passProbePhiLeg%02d/D",iLeg).Data());
}

b_counting            = ntuple->Branch("counting"            , &counting);
b_count_tag           = ntuple->Branch("count_tag"           , &count_tag);
b_count_probe         = ntuple->Branch("count_probe"         , &count_probe);
b_count_paasprobeleg1 = ntuple->Branch("count_passprobeleg1" , &count_passprobeleg1 , "count_passprobeleg1/D");
b_count_paasprobeleg2 = ntuple->Branch("count_passprobeleg2" , &count_passprobeleg2 , "count_passprobeleg2/D");
b_passL1PTNorm        = ntuple->Branch("passL1PTNorm"        , &passL1PTNorm        , "passL1PTNorm/D");
b_passL1PTIso         = ntuple->Branch("passL1PTIso"         , &passL1PTIso         , "passL1PTIso/D");
b_passL1PTER          = ntuple->Branch("passL1PTER"          , &passL1PTER          , "passL1PTER/D");
b_passL1PT            = ntuple->Branch("passL1PT"            , &passL1PT            , "passL1PT/D");
b_passL1EtaNorm       = ntuple->Branch("passL1EtaNorm"       , &passL1EtaNorm       , "passL1EtaNorm/D");
b_passL1EtaIso        = ntuple->Branch("passL1EtaIso"        , &passL1EtaIso        , "passL1EtaIso/D");
b_passL1EtaER         = ntuple->Branch("passL1EtaER"         , &passL1EtaER         , "passL1EtaER/D");
b_passL1Eta           = ntuple->Branch("passL1Eta"           , &passL1Eta           , "passL1Eta/D");

counting            = 0;
count_tag           = 0;
count_probe         = 0;
count_passprobeleg1 = 0;
count_passprobeleg2 = 0;

}//}}}

// ------------ method called once each job just after ending the event loop  ------------
void 
runEffAnaMiniAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
runEffAnaMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
edm::ParameterSetDescription desc;
desc.setUnknown();
descriptions.addDefault(desc);
}

void runEffAnaMiniAOD::reset()
{//{{{
int BIGNUMBER = -999;
Nvtx              = BIGNUMBER;
ZMass             = BIGNUMBER;
tagPT             = BIGNUMBER;
tagEta            = BIGNUMBER;
tagPhi            = BIGNUMBER;
tagMVA80Iso15     = BIGNUMBER;
tagMVA80Iso16     = BIGNUMBER;
tagMVA90Iso15     = BIGNUMBER;
tagMVA90Iso16     = BIGNUMBER;
tagCB             = BIGNUMBER;
tagCBTriChg       = BIGNUMBER;
probePT           = BIGNUMBER;
probeEta          = BIGNUMBER;
probePhi          = BIGNUMBER;
probeMVA80Iso15   = BIGNUMBER;
probeMVA80Iso16   = BIGNUMBER;
probeMVA90Iso15   = BIGNUMBER;
probeMVA90Iso16   = BIGNUMBER;
probeCB           = BIGNUMBER;
probeCBTriChg     = BIGNUMBER;
for (int i = 0; i < nLegs; i++) {
    passProbePTLegs [i] = BIGNUMBER;
    passProbeEtaLegs[i] = BIGNUMBER;
    passProbePhiLegs[i] = BIGNUMBER;
}
passL1PTNorm      = BIGNUMBER;
passL1PTIso       = BIGNUMBER;
passL1PTER        = BIGNUMBER;
passL1PT          = BIGNUMBER;
passL1EtaNorm     = BIGNUMBER;
passL1EtaIso      = BIGNUMBER;
passL1EtaER       = BIGNUMBER;
passL1Eta         = BIGNUMBER;
//nPU               = BIGNUMBER;
//nPUTrue           = BIGNUMBER;
PUWeight          = 1.;
lumi_weight       = 1.;
}//}}}

bool runEffAnaMiniAOD::passFilter(pat::TriggerObjectStandAlone & Triggerobj, std::string & FilterName)
{//{{{
for (unsigned h = 0; h < Triggerobj.filterLabels().size(); ++h) {
//cout << Triggerobj.filterLabels()[h] << endl;
    if (FilterName.compare(Triggerobj.filterLabels()[h]) == 0)
        return true;
}
return false;
}//}}}

double runEffAnaMiniAOD::DeltaRtrig(const pat::Electron & el, std::vector < pat::TriggerObjectStandAlone > object)
{//{{{

std::vector < double >coll;
for (unsigned int i = 0; i < object.size(); i++) {
    float dR = deltaR(el.superCluster()->eta(), el.phi(), object[i].eta(), object[i].phi());

    coll.push_back(dR);
}
std::sort(coll.begin(), coll.end());
if (coll.size() <= 0)
    return 999.0;
else
    return coll.at(0);
}//}}}

bool runEffAnaMiniAOD::hasZMother(const reco::GenParticle p)
{//{{{
bool foundZ(false);
const reco::Candidate * part = (p.mother());
// loop on the mother particles to check if is has a W has mother
while ((part->numberOfMothers() > 0)) {
    const reco::Candidate * MomPart = part->mother();
    if ((fabs(MomPart->pdgId()) == 23)) {
        foundZ = true;
        break;
    }
    part = MomPart;
}
return foundZ;
}//}}}

//define this as a plug-in
DEFINE_FWK_MODULE(runEffAnaMiniAOD);
