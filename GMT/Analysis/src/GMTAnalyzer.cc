#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>
#include <cmath> 
#include <iostream>
#include "TCanvas.h"
#include "TreeAnalyzer.h"
#include "GMTAnalyzer.h"
#include "EventProxyOMTF.h"
#include "EventProxyOMTFNANOAOD.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TMath.h"
#include "TROOT.h"
using namespace TMath;
using namespace std;
std::vector<double> ptRanges = {0.,   1.,   2.,   3.,   4.,   
                                  5.,   6.,   7.,   8.,   9., 
                                  10.,  11.,  12.,  13.,  14.,
                                  15.,  16.,  17.,  18.,  19.,
                                  20.,  21.,  22.,  23.,  24.,
                                  25.,  26.,  28.,  30.,  32.,  34.,
                                  36.,  38.,  40.,  50.,  60.,
                                  70.,  80.,  90.,  100., 200., 99999};
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
GMTAnalyzer::GMTAnalyzer(const std::string & aName):Analyzer(aName){

  selectionFlavours_.push_back(aName);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
GMTAnalyzer::~GMTAnalyzer(){
  
  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::initialize(TDirectory* aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;
  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new GMTHistograms(aDir, selectionFlavours_);

}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
Analyzer* GMTAnalyzer::clone() const{

  GMTAnalyzer* clone = new GMTAnalyzer(name());
  clone->setHistos(myHistos_);
  return clone;

};
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::finalize(){

  myHistos_->finalizeHistograms();

}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
bool GMTAnalyzer::passQuality(const L1Obj & aL1Cand,
			       const std::string & sysType,
			       const std::string & selType){
  
  bool lowPtVeto = false;

   if(sysType.find("uGMT")!=std::string::npos){
     
     return aL1Cand.type==L1Obj::uGMT && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;
   }
   else if(sysType.find("OMTF")!=std::string::npos){
     return true;
   }   
   return false;
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillTurnOnCurve(const MuonObj & aMuonCand,
                                  const int & iPtCut,
				                          const std::string & sysType,
				                          const std::string & selType){
  
 //int is important for histo name construction
  int ptCut = GMTHistograms::ptBins[iPtCut];
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2DGmt"+selType;
  if(sysType=="OMTF") {   
    hName = "h2DOMTF"+selType;
  }
  if(sysType=="uGMT") {   
    hName = "h2DuGMT"+selType;
  }
  if(sysType=="EMTF") {   
    hName = "h2DEMTF"+selType;
  }

  ///Find the best matching L1 candidate
  double deltaR = 0.15;
  L1Obj selectedCand;

   
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);
    if(!pass) continue;    
    
    double phiValue = aCand.phiValue();
    if(phiValue>M_PI) phiValue-=2*M_PI;
    
    double dEta = std::abs(aMuonCand.l1eta()-aCand.etaValue());
    double dPhi = std::abs(aMuonCand.l1phi()-phiValue);
    if(dPhi>2*M_PI) dPhi=-2*M_PI;
    double delta = sqrt(dEta*dEta + dPhi*dPhi);
    if(delta<deltaR && selectedCand.ptValue()<aCand.ptValue()){
      deltaR = delta;
      selectedCand = aCand;      
    }    
  }

  bool passPtCut = selectedCand.ptValue()>=ptCut && selectedCand.ptValue()>0;
   
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtOMTF";
  myHistos_->fill2DHistogram(tmpName, aMuonCand.pt(), selectedCand.ptValue());
  
  //Generic eff vs selected variable calculated for muons on plateau
  if(!selType.size() && aMuonCand.pt()<ptCut+20) return;
  tmpName = hName+"EtauGMT"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.eta(), passPtCut);

  tmpName = hName+"PhiuGMT"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.phi(), passPtCut);
  


}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
 
void GMTAnalyzer::fillRateHisto(const MuonObj & aRecoMuon,
                                const std::string & sysType,
				                        const std::string & selType){

  //Generator level information is not available for the neutrino sample
  if(name()=="NU_RATEAnalyzer" && aRecoMuon.pt()>0.0) return;

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2D"+sysType+"Rate"+selType;

  L1Obj selectedCand;
  for(auto aCand: myL1Coll){

    bool pass = passQuality(aCand ,sysType, selType);    
    if(pass && selectedCand.ptValue()<aCand.ptValue()) selectedCand = aCand;

  }

  bool pass = selectedCand.ptValue()>=20;
  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon.pt(),selectedCand.ptValue());
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon.pt(),pass*aRecoMuon.eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon.pt(),pass*aRecoMuon.pt()+(!pass)*(-100));
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillHistosForRecoMuon(const MuonObj & aRecoMuon){   

  bool isGMTAcceptance = fabs(aRecoMuon.eta())<2.4;
  if(!isGMTAcceptance) return;

  bool isOMTFAcceptance = fabs(aRecoMuon.eta())>0.83 && fabs(aRecoMuon.eta())<1.24;
  if(!isOMTFAcceptance) return;
  
  myHistos_->fill1DHistogram("h1DPtProbe", aRecoMuon.pt());
  myHistos_->fill1DHistogram("h1DAbsEtaProbe", std::abs(aRecoMuon.eta()));
  
  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurve(aRecoMuon, iCut, "OMTF", selType);
      fillTurnOnCurve(aRecoMuon, iCut, "uGMT", selType);
  }

  int iCut = 18;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
    
    if(iType==0) pass = aRecoMuon.pt()>ptCut + 20;
    else if(iType==1) pass = aRecoMuon.pt()>ptCut && aRecoMuon.pt()<(ptCut+5);
    else if(iType==2) pass = aRecoMuon.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(aRecoMuon, iCut, "OMTF", selType);
    fillTurnOnCurve(aRecoMuon, iCut, "uGMT", selType);
  }
}
//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
double GMTAnalyzer::customDeltaR(TLorentzVector T1, TLorentzVector T2){
		double dphi = T1.Phi() - T2.Phi();
                double deta = T1.Eta() - T2.Eta();
                if(dphi > TMath::Pi()) dphi = 2.0*TMath::Pi() - dphi;
                if(dphi < 0) dphi = -1.0*dphi;
                double delR =  TMath::Sqrt( TMath::Power(deta,2) + TMath::Power(dphi,2));
                return delR ;

}
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool GMTAnalyzer::analyze(const EventProxyBase& iEvent){
   
   clear();
   
  bool useNanoAOD = false; 

if (useNanoAOD) {
    const EventProxyOMTFNANOAOD& myProxy = dynamic_cast<const EventProxyOMTFNANOAOD&>(iEvent);
    const_cast<EventProxyOMTFNANOAOD&>(myProxy).fillL1ObjColl();
    const_cast<EventProxyOMTFNANOAOD&>(myProxy).fillMuonObjColl();
    myMuonObjColl = const_cast<EventProxyOMTFNANOAOD&>(myProxy).getRecoMuonObjColl();
    myL1ObjColl = const_cast<EventProxyOMTFNANOAOD&>(myProxy).getL1ObjColl();
/*
 //cannot perform a static_cast from a const EventProxyBase& to EventProxyOMTFNANOAOD& because the source object (iEvent) is const, while the destination type (EventProxyOMTFNANOAOD&) is non-const.
 //To resolve this issue, you can use const_cast to remove the const qualifier temporarily. 

    EventProxyOMTFNANOAOD& myProxy = static_cast< EventProxyOMTFNANOAOD&>(iEvent);
    myProxy.fillL1ObjColl();  
    myProxy.fillMuonObjColl();
    myMuonObjColl = myProxy.getRecoMuonObjColl();
    myL1ObjColl = myProxy.getL1ObjColl();*/
} else {
    const EventProxyOMTF& myProxy = static_cast<const EventProxyOMTF&>(iEvent);
    myEventId = myProxy.getEventId();
    myMuonObjColl = myProxy.getRecoMuonObjColl();
    myL1ObjColl = myProxy.getL1ObjColl();
}

  const std::vector<MuonObj> & myMuonColl = myMuonObjColl->getMuonObjs();
  

  if(myMuonColl.size() < 2 )return false;

  MuonObj aTagCand =  myMuonColl.at(0);
  bool tagPass = aTagCand.pt()>20 && aTagCand.matchedisohlt();
  if(!tagPass) return true;
  
  
  tagFourVector.SetPtEtaPhiM(aTagCand.pt(), aTagCand.eta(), aTagCand.phi(), nominalMuonMass);
  myHistos_->fill1DHistogram("h1DPtTag", aTagCand.pt());
  myHistos_->fill1DHistogram("h1DAbsEtaTag", std::abs(aTagCand.eta()));

  MuonObj aProbeCand;
  double m_Z = 91.1876;
  double deltaM_Z = 20;
  double tmpDelta = 2*deltaM_Z;
  for (auto aMuonCand: myMuonColl){   
      randomMuonLeg.SetPtEtaPhiM(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), nominalMuonMass);
      tmpDelta = std::abs((tagFourVector+randomMuonLeg).M()-m_Z);
      if(aMuonCand.tightID() && tmpDelta<deltaM_Z){
      deltaM_Z = tmpDelta;
      aProbeCand = aMuonCand;   
    }
  }
  if(aProbeCand.pt()<1) return true;

  probeFourVector.SetPtEtaPhiM(aProbeCand.pt(), aProbeCand.eta(), aProbeCand.phi(), nominalMuonMass);
  myHistos_->fill1DHistogram("h1DDiMuonMassTagProbe",(tagFourVector+probeFourVector).M());   
  fillHistosForRecoMuon(aProbeCand);

  double deltaRCut = 0.6;
  for (auto aMuonCand: myMuonColl){   
      TLorentzVector puMuonLeg;
      puMuonLeg.SetPtEtaPhiM(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), nominalMuonMass);
      if(aMuonCand.tightID() 
               && customDeltaR(puMuonLeg,tagFourVector) > deltaRCut 
               && customDeltaR(puMuonLeg,probeFourVector) > deltaRCut){
        fillHistosForRecoMuon(aMuonCand);
      }
    }

    
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
