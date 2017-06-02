#include <DataFormats/Common/interface/Ptr.h>
#include "DataFormats/Common/interface/Handle.h"
//#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" 


class MyPatElectronMVAIdSelector
{
public:
    MyPatElectronMVAIdSelector();
    MyPatElectronMVAIdSelector(const edm::ParameterSet&);

    void setRho(double rho);
    void setBS(edm::Handle<reco::BeamSpot>&);
    void setFGV(reco::VertexCollection::const_iterator&);
    void setConversion(edm::Handle<reco::ConversionCollection>&);
    void setPassMVA80Xwp80(edm::Handle<edm::ValueMap<bool> >&);

    bool passId      (const edm::Ptr<const pat::Electron>);
    bool passId      (const pat::Electron&);
    bool passTriChg      (const edm::Ptr<const pat::Electron>);
    bool passTriChg      (const pat::Electron&);
    bool passIso      (const edm::Ptr<const pat::Electron> , int IsoVer);
    bool passIso      (const pat::Electron& , int IsoVer);
    bool passHLTId (const edm::Ptr<const pat::Electron>);
    bool passHLTId (const pat::Electron&);

private:
    /* data */
    //double isolationUpCut = 0 ;

    bool tripleCharge;
    double trkPFIso ;
    double ptElec ;
    double etaElec ;
    double hOverE ;
    double full5x5sigmaIetaIeta ;
    double dEtaSeed ;
    double dPhiIn ;
    double ecalPFIso ;
    double hcalPFIso ;
    double chisq ;
    double jetRho ;
    double ooEmooP ;

    double effArea ;
    bool passConversionVeto ;
    double expectedMissingInnerHits ;
    edm::Handle<reco::BeamSpot> theBeamSpot; 
    edm::Handle<reco::ConversionCollection> conversions; 
    reco::VertexCollection::const_iterator firstGoodVertex; 
    edm::Handle<edm::ValueMap<bool> > mvaTightIdDecision;

    double isolation ;
    double d0 ;
    double dz ;

    // functions
    void initEle_(const pat::Electron&);
    bool passHLTId_();
    bool passId_();
    bool passIso_(int IsoVer);

};

MyPatElectronMVAIdSelector::MyPatElectronMVAIdSelector(){
}
MyPatElectronMVAIdSelector::MyPatElectronMVAIdSelector(const edm::ParameterSet &ps){
    //isolationUpCut = ps.getParameter<double>("isolationUpCut");
}

void MyPatElectronMVAIdSelector::initEle_(const pat::Electron &el){

    
    tripleCharge         = el.chargeInfo().isGsfCtfScPixConsistent;
    trkPFIso             = el.trackIso();
    ptElec               = el.pt();
    etaElec              = el.eta();
    hOverE               = el.hcalOverEcal();
    full5x5sigmaIetaIeta = el.full5x5_sigmaIetaIeta();
    dEtaSeed             = el.deltaEtaSeedClusterTrackAtCalo();
    dPhiIn               = el.deltaPhiSuperClusterTrackAtVtx();
    ecalPFIso            = el.ecalPFClusterIso();
    hcalPFIso            = el.hcalPFClusterIso();
    chisq                = el.gsfTrack()->normalizedChi2();
    if( el.ecalEnergy() == 0 ){
        ooEmooP = ( 1e30 );
    }else if( !std::isfinite(el.ecalEnergy())){
        ooEmooP = ( 1e30 );
    }else{
        ooEmooP = ( fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() ) );
    }

    edm::FileInPath eaConstantsFile("RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt");
    EffectiveAreas  effectiveAreas(eaConstantsFile.fullPath());
    effArea = effectiveAreas.getEffectiveArea(fabs(el.superCluster()->eta()));

    passConversionVeto       = !ConversionTools::hasMatchedConversion(el, conversions, theBeamSpot->position());
    expectedMissingInnerHits = (el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
    isolation                = (el.pfIsolationVariables().sumChargedHadronPt+std::max(el.pfIsolationVariables().sumNeutralHadronEt + el.pfIsolationVariables().sumPhotonEt-jetRho*effArea,0.))/el.pt();
    d0                       = (-1) * el.gsfTrack()->dxy(firstGoodVertex->position());
    dz                       = el.gsfTrack()->dz(firstGoodVertex->position());
}

void MyPatElectronMVAIdSelector::setFGV(reco::VertexCollection::const_iterator &fGV){
    firstGoodVertex = fGV;
}

void MyPatElectronMVAIdSelector::setBS(edm::Handle<reco::BeamSpot> &bs){
    theBeamSpot = bs; 
}

void MyPatElectronMVAIdSelector::setConversion(edm::Handle<reco::ConversionCollection> &cv){
    conversions = cv;
}

void MyPatElectronMVAIdSelector::setRho(double rho){
    jetRho = rho;
}

bool MyPatElectronMVAIdSelector::passId_(){
    if( /*MVAId cut*/
        abs(etaElec) < 2.5
        &&((abs(etaElec) <= 1.479
        //&& isolation < 0.571 //Iso16
        //&& isolation < 0.0354 //Iso15
        //&& (*mvaTightIdDecision)[elPtr] == 1
        )||(abs(etaElec) > 1.479 && abs(etaElec) < 2.5
        //&& isolation < 0.588 //Iso16
        //&& isolation < 0.0646 //Iso15
        //&& (*mvaTightIdDecision)[elPtr] == 1
        ))){
        return true;
    }else{   
        return false;
    }        
}

bool MyPatElectronMVAIdSelector::passTriChg(const pat::Electron &el){
    initEle_(el);
    return tripleCharge;
}
bool MyPatElectronMVAIdSelector::passTriChg(const edm::Ptr<const pat::Electron> elPtr){
    initEle_(*elPtr);
    return tripleCharge;
}
//PF isolation 2016/2017
bool MyPatElectronMVAIdSelector::passIso_(int IsoVer){
    if( abs(etaElec) < 2.5
        &&((abs(etaElec) <= 1.479
        && (isolation < 0.571 || IsoVer == 16) && (isolation < 0.0354 || IsoVer ==15)
        )||(abs(etaElec) > 1.479 && abs(etaElec) < 2.5
        && (isolation < 0.588 || IsoVer == 16) && (isolation < 0.0646 || IsoVer ==15)
        ))){
        return true;
    }else{   
        return false;
    }        
}

bool MyPatElectronMVAIdSelector::passId(const pat::Electron &el){
    if (!passHLTId(el)) return false;
    return passId_();
}

bool MyPatElectronMVAIdSelector::passId(const edm::Ptr<const pat::Electron> elPtr){
    if (!passHLTId(elPtr)) return false;
    return passId_();
}

bool MyPatElectronMVAIdSelector::passIso(const pat::Electron &el, int IsoVer){
    return passIso_(IsoVer);
}

bool MyPatElectronMVAIdSelector::passIso(const edm::Ptr<const pat::Electron> elPtr, int IsoVer){
    return passIso_(IsoVer);
}
//HLT safe trigger cut
bool MyPatElectronMVAIdSelector::passHLTId_(){
    if(/*loose cut*/
        abs(etaElec) < 2.5
        && trkPFIso/ptElec < 0.08
        &&(((abs(etaElec) <= 1.479)
        &&( hOverE < 0.060        
        && full5x5sigmaIetaIeta < 0.011
        && abs(dEtaSeed)  < 0.004         
        && abs(dPhiIn) < 0.020               
        && abs(ooEmooP) < 0.013                
        && (ecalPFIso - jetRho * 0.165)/ptElec < 0.160
        && (hcalPFIso - jetRho * 0.060)/ptElec < 0.120
        && abs(dz)      < 0.1
        && abs(d0)< 0.05         
        && passConversionVeto   
        && expectedMissingInnerHits< 1  
        ))||((abs(etaElec) > 1.479 && abs(etaElec) < 2.5)
        && (hOverE              < 0.065           
        && full5x5sigmaIetaIeta    < 0.031      
        && abs(ooEmooP) < 0.013                    
        && (ecalPFIso -jetRho * 0.132)/ptElec < 0.120 
        && (hcalPFIso -jetRho * 0.131)/ptElec < 0.120  
        && abs(chisq) < 3      
        && abs(dz)      < 0.2            
        && abs(d0)< 0.1            
        && passConversionVeto     
        && expectedMissingInnerHits<1 
        )))){
        return true;
    }else{   
        return false;
    }        
}

bool MyPatElectronMVAIdSelector::passHLTId(const pat::Electron &el){
    initEle_(el);
    return passHLTId_();
}

bool MyPatElectronMVAIdSelector::passHLTId(const edm::Ptr<const pat::Electron> elPtr){
    initEle_(*elPtr);
    return passHLTId_();
}

