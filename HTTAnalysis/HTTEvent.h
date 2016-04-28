#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include <map>
#include <vector>
#include <bitset> 

//
// class declaration
//
class Wevent{

    int run_ = 0;
    float lumi_ = 0;
    long int event_ = 0;
    int nup_ = -999, npu_ = -999, npv_ = -999; 
    int paircount_ = -1;
    float genevtweight_ = 1;
    int sample_ = -1;
    int bosonId_ = 0;

    ///Boson (H, Z, W) decay mode
    int decayModeBoson_ = -1;
    
    ///Tau decay modes (if available)
    int decModeMinus_, decModePlus_;

    ///PCA analysis data members
    ///Primary Vertices recontructed with different methods

    //Generated PV position 
    TVector3 genPV_;

    //PV stored in miniAOD
    TVector3 thePV_;

    //Reco PV closes to the Gen PV    
    TVector3 bestPV_;

    ///PV recontructed from PF candidates, refitted
    TVector3 refitPfPV_;

    ///PV recontructed from PF candidates, refitted without beamspot (BS) requirement
    TVector3 refitPfPVNoBS_;

    ///Flag marking if refit was successfull
    bool isRefit_;

    ///Number of tracks used in the refit
    int nTracksInRefit_;

  public:

    enum processenum {DATA = 0, DY = 1, WJets=2, TTbar=3, h=3, H=4, A=5, QCD=6};
    
    Wevent();
    ~Wevent();
    void run(int x){run_ = x;}
    void lumi(float x){lumi_ = x;}
    void event(long int x){event_ = x;}
    void nup(int x){nup_ = x;}
    void npu(int x){npu_ = x;}
    void npv(int x){npv_ = x;}
    void paircount(int x){paircount_ = x;}
    void genevtweight(float x){genevtweight_ = x;}
    void sample(int x){sample_ = x;}
    void bosonId(int x){bosonId_ = x;}
    void decModeMinus(int x){decModeMinus_ = x;}
    void decModePlus(int x){decModePlus_ = x;}
    void decayModeBoson(int x){decayModeBoson_ = x;}

    ///Set generated PV
    void genPV(const TVector3 & aPV) {genPV_ = aPV;}

    ///Set PV stored in miniAOD
    void thePV(const TVector3 & aPV) {thePV_ = aPV;}

    //Reco PV closest to the Gen PV
    void bestPV(const TVector3 & aPV) {bestPV_ = aPV;}

    //Set PV refitted using BS
    void refitPfPV(const TVector3 & aPV) {refitPfPV_ = aPV;}

    //Set PV refitted without BS
    void refitPfPVNoBS(const TVector3 & aPV) {refitPfPVNoBS_ = aPV;}

    ///Set isRefit bool.
    void isRefit(bool aBit){isRefit_ = aBit;};

    ///Set number of trascks used in refit.
    void nTracksInRefit(const int & nTracks) {nTracksInRefit_ = nTracks;};

    ///Reset class data members
    void clear();

    int run()const{return run_;}
    float lumi()const{return lumi_;}
    long int event()const{return event_;}
    int nup()const{return nup_;}
    int npu()const{return npu_;}
    int npv()const{return npv_;}
    int paircount()const{return paircount_;}
    float genevtweight()const{return genevtweight_;}
    int sample()const{return sample_;}
    int bosonId()const{return bosonId_;}
    int decModeMinus()const{return decModeMinus_;}
    int decModePlus()const{return decModePlus_;}
    int decayModeBoson()const{return decayModeBoson_;}

    ///Get generated PV 
    const TVector3 & genPV() const {return genPV_;}

    ///Get PV stored in miniAOD
    const TVector3 & thePV() const {return thePV_;}

    //Get reco PV closest to the Gen PV    
    const TVector3 & bestPV() const {return bestPV_;}

    //Get PV refitted using BS
    const TVector3 & refitPfPV() const {return refitPfPV_;}

    //Get PV refitted without BS
    const TVector3 & refitPfPVNoBS() const {return refitPfPVNoBS_;}    

    ///Set isRefit bool.
    const bool isRefit() const {return isRefit_;};

    ///Set number of trascks used in refit.
    const int nTracksInRefit() const {return nTracksInRefit_;}

};

enum tauidenum {
     decayModeFinding, 
     decayModeFindingNewDMs, 
     byCombinedIsolationDeltaBetaCorrRaw3Hits,
     byLooseCombinedIsolationDeltaBetaCorr3Hits,
     byMediumCombinedIsolationDeltaBetaCorr3Hits,
     byTightCombinedIsolationDeltaBetaCorr3Hits,
     chargedIsoPtSum,
     neutralIsoPtSum,
     puCorrPtSum,
     againstMuonLoose3,
     againstMuonTight3,
     againstElectronVLooseMVA5,
     againstElectronLooseMVA5,
     againstElectronMediumMVA5,
     againstElectronTightMVA5,
     againstElectronVTightMVA5,
     byIsolationMVA3newDMwoLTraw,
     byIsolationMVA3oldDMwoLTraw,
     byIsolationMVA3newDMwLTraw,
     byIsolationMVA3oldDMwLTraw,

};
class Wtau{
    
    float pt_ = -999;
    float phi_ = -999;
    float eta_ = -999;
    float mass_ = -999;
    int charge_ = -999;
    int decayMode_ = -999;
    float mt_ = -999;
    float d0_ = -999;
    float dz_ = -999;
    float idSF_ = 1.0;
    int   trackPattern_ = 0;
    float triggerSF_ = 1.0;

    ///Leading tau track four momemntum.
    TLorentzVector leadingTk_;

    ///Secondary vertex position (from GEN)
    TVector3 sv_;

    ///PCA vector (vector from PV to PCA)
    TVector3 nPCA_;

    ///PCA vectors calculated using different PV estimates
    TVector3 nPCAAODvx_, nPCAGenvx_, nPCARefitvx_;
    
    std::vector<float> tauID_ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 

  public:
    Wtau();
    ~Wtau();
    void clear();

    void pt(float x){pt_ = x;}
    void phi(float x){phi_ = x;}
    void eta(float x){eta_ = x;}
    void mass(float x){mass_ = x;}
    void charge(float x){charge_ = x;}
    void decayMode(int x){decayMode_ = x;}
    void mt(float x){mt_ = x;}
    void tauID(tauidenum y, float x){tauID_[y] = x;}
    void d0(float x){d0_ = x;}
    void dz(float x){dz_ = x;}
    void trackPattern(int x){trackPattern_ = x;}
    void sv(const TVector3 & x){sv_ = x;}
    void idSF(float x){idSF_ = x;}
    void triggerSF(float x){triggerSF_ = x;}

    ///Set tau leading charged track.
    void leadingTk(const TLorentzVector & a4v) {leadingTk_ = a4v;};

    ///Set PCA vector calculated using PV stored in AOD
    void nPCA(const TVector3 & a3v) {nPCA_ = a3v;};

    ///Set PCA vector calculated using PV selected using PF weights
    ///relaculated from miniAOD
    void nPCAAODvx(const TVector3 & a3v) {nPCAAODvx_ = a3v;};

    ///Set PCA vector calculated using generated PV 
    void nPCAGenvx(const TVector3 & a3v) {nPCAGenvx_ = a3v;};

    ///Set PCA vector calculated using refitted PV 
    void nPCARefitvx(const TVector3 & a3v) {nPCARefitvx_ = a3v;};

    float pt()const{return pt_;}
    float phi()const{return phi_;}
    float eta()const{return eta_;}
    float mass()const{return mass_;}
    float charge()const{return charge_;}
    int decayMode()const{return decayMode_;}
    float mt()const{return mt_;}
    float d0()const{return d0_;}
    float dz()const{return dz_;}
    int trackPattern()const{return trackPattern_;}
    TVector3 sv()const{return sv_;}
    float idSF()const{return idSF_;}
    float triggerSF()const{return triggerSF_;}

    float tauID(tauidenum y){return tauID_[y];}

    ///Get leading charged track
    const TLorentzVector & leadingTk() const {return leadingTk_;};

    ///Get 4-momentum
    TLorentzVector p4() const {TLorentzVector a4v; a4v.SetPtEtaPhiM(pt_, eta_, phi_, mass_); return a4v;};

    ///Get PCA vector calculated using PV stored in AOD
    const TVector3 & nPCA() {return nPCA_;};

    ///Get PCA vector calculated using PV selected using PF weights
    ///relaculated from miniAOD
    const TVector3 & nPCAAODvx() {return nPCAAODvx_;};

    ///Get PCA vector calculated using generated PV 
    const TVector3 & nPCAGenvx() {return nPCAGenvx_;};

    ///Get PCA vector calculated using refitted PV 
    const TVector3 & nPCARefitvx() {return nPCARefitvx_;};

};
typedef std::vector<Wtau> WtauCollection;


class Wmu{
    
    float pt_ = -999;
    float phi_ = -999;
    float eta_ = -999;
    float mass_ = -999;
    int charge_ = -999;
    float mt_ = -999;
    float d0_ = -999;
    float dz_ = -999;
    int trackPattern_=0;
    float isLooseMuon_ = -999;
    float isTightMuon_ = -999;
    float isHighPtMuon_ = -999;
    float isMediumMuon_ = -999;
    float isTightnovtxMuon_ = -999;
    float iso_ = -999;
    float idSF_ = 1.0;
    float triggerSF_ = 1.0;
     

    ///Secondary vertex position (from GEN)
    TVector3 sv_;

    ///PCA vector (vector from PV to PCA)
    TVector3 nPCA_;

    ///PCA vectors calculated using different PV estimates
    TVector3 nPCAAODvx_, nPCAGenvx_, nPCARefitvx_;

  public:
    Wmu();
    ~Wmu();
    void clear();

    void pt(float x){pt_ = x;}
    void phi(float x){phi_ = x;}
    void eta(float x){eta_ = x;}
    void mass(float x){mass_ = x;}
    void charge(float x){charge_ = x;}
    void mt(float x){mt_ = x;}
    void d0(float x){d0_ = x;}
    void dz(float x){dz_ = x;}
    void trackPattern(int x){trackPattern_ = x;}
    void isLooseMuon(float x){isLooseMuon_ = x;}
    void isTightMuon(float x){isTightMuon_ = x;}
    void isHighPtMuon(float x){isHighPtMuon_ = x;}
    void isMediumMuon(float x){isMediumMuon_ = x;}
    void isTightnovtxMuon(float x){isTightnovtxMuon_ = x;}
    void iso(float x){iso_ = x;}

    void idSF(float x){idSF_ = x;}
    void triggerSF(float x){triggerSF_ = x;}

    ///Set PCA vector calculated using PV stored in AOD
    void nPCA(const TVector3 & a3v) {nPCA_ = a3v;};

    ///Set PCA vector calculated using PV selected using PF weights
    ///relaculated from miniAOD
    void nPCAAODvx(const TVector3 & a3v) {nPCAAODvx_ = a3v;};

    ///Set PCA vector calculated using generated PV 
    void nPCAGenvx(const TVector3 & a3v) {nPCAGenvx_ = a3v;};

    ///Set PCA vector calculated using refitted PV 
    void nPCARefitvx(const TVector3 & a3v) {nPCARefitvx_ = a3v;};

    float pt()const{return pt_;}
    float phi()const{return phi_;}
    float eta()const{return eta_;}
    float mass()const{return mass_;}
    float charge()const{return charge_;}
    float mt()const{return mt_;}
    float d0()const{return d0_;}
    float dz()const{return dz_;}
    int trackPattern()const {return trackPattern_;}
    float isLooseMuon()const{return isLooseMuon_;}
    float isTightMuon()const{return isTightMuon_;}
    float isHighPtMuon()const{return isHighPtMuon_;}
    float isMediumMuon()const{return isMediumMuon_;}
    float isTightnovtxMuon()const{return isTightnovtxMuon_;}
    float iso()const{return iso_;}
    float idSF()const{return idSF_;}
    float triggerSF()const{return triggerSF_;}

    ///
    ///Get leading charged track
    TLorentzVector leadingTk() const {TLorentzVector a4v; a4v.SetPtEtaPhiM(pt_, eta_, phi_, 0.105658); return a4v;};
    
    ///Get 4-momentum
    TLorentzVector p4() const {TLorentzVector a4v; a4v.SetPtEtaPhiM(pt_, eta_, phi_, 0.105658); return a4v;};

    ///Get PCA vector calculated using PV stored in AOD
    const TVector3 & nPCA() {return nPCA_;};

    ///Get PCA vector calculated using PV selected using PF weights
    ///relaculated from miniAOD
    const TVector3 & nPCAAODvx() {return nPCAAODvx_;};

    ///Get PCA vector calculated using generated PV 
    const TVector3 & nPCAGenvx() {return nPCAGenvx_;};

    ///Get PCA vector calculated using refitted PV 
    const TVector3 & nPCARefitvx() {return nPCARefitvx_;};


};
typedef std::vector<Wmu> WmuCollection;

class Welectron{

    float pt_ = -999;
    float phi_ = -999;
    float eta_ = -999;
    float mass_ = -999;
    int charge_ = -999;

    ///Secondary vertex position (from GEN)
    TVector3 sv_;
    
    ///PCA vector (vector from PV to PCA)
    TVector3 nPCA_;
    
    ///PCA vectors calculated using different PV estimates
    TVector3 nPCAAODvx_, nPCAGenvx_, nPCARefitvx_;

  public:
    Welectron();
    ~Welectron();
    void clear();
    
    void pt(float x){pt_ = x;}
    void phi(float x){phi_ = x;}
    void eta(float x){eta_ = x;}
    void mass(float x){mass_ = x;}
    void charge(float x){charge_ = x;}

    ///Set PCA vector calculated using PV stored in AOD
    void nPCA(const TVector3 & a3v) {nPCA_ = a3v;};

    ///Set PCA vector calculated using PV selected using PF weights
    ///relaculated from miniAOD
    void nPCAAODvx(const TVector3 & a3v) {nPCAAODvx_ = a3v;};

    ///Set PCA vector calculated using generated PV 
    void nPCAGenvx(const TVector3 & a3v) {nPCAGenvx_ = a3v;};

    ///Set PCA vector calculated using refitted PV 
    void nPCARefitvx(const TVector3 & a3v) {nPCARefitvx_ = a3v;};
 
    float pt()const{return pt_;}
    float phi()const{return phi_;}
    float eta()const{return eta_;}
    float mass()const{return mass_;}
    float charge()const{return charge_;}

    ///Get leading charged track
    TLorentzVector leadingTk() const {TLorentzVector a4v; a4v.SetPtEtaPhiM(pt_, eta_, phi_, 0.); return a4v;};
    
    ///Get 4-momentum
    TLorentzVector p4() const {TLorentzVector a4v; a4v.SetPtEtaPhiM(pt_, eta_, phi_, 0.); return a4v;};

    ///Get PCA vector calculated using PV stored in AOD
    const TVector3 & nPCA() {return nPCA_;};

    ///Get PCA vector calculated using PV selected using PF weights
    ///relaculated from miniAOD
    const TVector3 & nPCAAODvx() {return nPCAAODvx_;};

    ///Get PCA vector calculated using generated PV 
    const TVector3 & nPCAGenvx() {return nPCAGenvx_;};

    ///Get PCA vector calculated using refitted PV 
    const TVector3 & nPCARefitvx() {return nPCARefitvx_;};

};
typedef std::vector<Welectron> WelectronCollection;

class Wmet{


    float metpx_ = -999;
    float metpt_ = -999;
    float metphi_ = -999;
    float mvacov00_ = -999;
    float mvacov01_ = -999;
    float mvacov10_ = -999;
    float mvacov11_ = -999;
 
  public:

    void metpx(float x){metpx_ = x;}
    void metpt(float x){metpt_ = x;}
    void metphi(float x){metphi_ = x;}
    void mvacov00(float x){mvacov00_ = x;}
    void mvacov01(float x){mvacov01_ = x;}
    void mvacov10(float x){mvacov10_ = x;}
    void mvacov11(float x){mvacov11_ = x;}

    float metpx()const{return metpx_;}
    float metpt()const{return metpt_;}
    float metphi()const{return metphi_;}
    float mvacov00()const{return mvacov00_;}
    float mvacov01()const{return mvacov01_;}
    float mvacov10()const{return mvacov10_;}
    float mvacov11()const{return mvacov11_;}
};
typedef std::vector<Wmet> WmetCollection;


enum triggersenum {
    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
    HLT_IsoMu17_eta2p1_LooseIsoPFTau20,
    HLT_IsoMu24_eta2p1,
    HLT_IsoMu27,
    HLT_IsoMu18,
    HLT_IsoMu22,
    HLT_IsoMu17_eta2p1,
    HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1

};

class Wpair{
    
     float svfit_ = -999;
     float diq_ = -999;
     float pth_ = -999;
     float ptvis_ = -999;
     float m_vis_ = -999;
    
    bool ChannelSelector_ = 0;
    bool PATPairSelector_ = 0;
    bool PairBaselineSelection_ = 0;
    bool PostSynchSelection_ = 0;

    std::vector<bool> triggers_ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
   
  public:
     Wpair();
     ~Wpair();

     ///Reset class data members
     void clear();
     
     void svfit(float x){svfit_ = x;}
     void diq(float x){diq_ = x;}
     void pth(float x){pth_ = x;}
     void ptvis(float x){ptvis_ = x;}
     void m_vis(float x){m_vis_ = x;}
    void ChannelSelector(float x){ChannelSelector_ = x;}
    void PATPairSelector(float x){PATPairSelector_ = x;}
    void PairBaselineSelection(float x){PairBaselineSelection_ = x;}
    void PostSynchSelection(float x){PostSynchSelection_ = x;}
    void trigger(triggersenum y, bool x){triggers_[y] = x;}


     float svfit()const{return svfit_;}
     float diq()const{return diq_;}
     float pth()const{return pth_;}
     float ptvis()const{return ptvis_;}
     float m_vis()const{return m_vis_;}
    bool ChannelSelector()const{return ChannelSelector_;}
    bool PATPairSelector()const{return PATPairSelector_;}
    bool PairBaselineSelection()const{return PairBaselineSelection_;}
    bool PostSynchSelection()const{return PostSynchSelection_;}
    bool trigger(triggersenum y){return triggers_[y];}
};
typedef std::vector<Wpair> WpairCollection;

class Wjet{

    float pt_ = -999;
    float eta_ = -999;
    float phi_ = -999;
    float id_ = -999;
    float bptag_ = -999;
    float csvtag_ = -999;
    float bjet_ = -999;
    float jecfactor_ = -999;
    float jetlooseID_ = -999;
    float pujetetaid_ = -999;
  public:

    void pt(float x){pt_ = x;}
    void eta(float x){eta_ = x;}
    void phi(float x){phi_ = x;}
    void id(float x){id_ = x;}
    void bptag(float x){bptag_ = x;}
    void csvtag(float x){csvtag_ = x;}
    void bjet(float x){bjet_ = x;}
    void jecfactor(float x){jecfactor_ = x;}
    void jetlooseID(float x){jetlooseID_ = x;}
    void pujetetaid(float x){pujetetaid_ = x;}

    float pt()const{return pt_;}
    float eta()const{return eta_;}
    float phi()const{return phi_;}
    float id()const{return id_;}
    float bptag()const{return bptag_;}
    float csvtag()const{return csvtag_;}
    float bjet()const{return bjet_;}
    float jecfactor()const{return jecfactor_;}
    float jetlooseID()const{return jetlooseID_;}
    float pujetetaid()const{return pujetetaid_;}
};
typedef std::vector<Wjet> WjetCollection;









#endif

