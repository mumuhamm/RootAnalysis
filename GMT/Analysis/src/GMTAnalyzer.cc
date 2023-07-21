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

   
   if (aName.find("_nanoAOD") != std::string::npos) {
        inputType = "nanoAOD";
    } else if (aName.find("_omtfTree") != std::string::npos) {
        inputType = "omtfTree";
    } else {
        std::cout << "Invalid process name!" << std::endl;
        return;
    }
    selectionFlavours_.push_back(aName);
    std::cout << "Input type: " << inputType << std::endl;
 
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
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool GMTAnalyzer::passQuality(const L1Obj & aL1Cand,
			       const std::string & sysType,
			       const std::string & selType){
  
  bool lowPtVeto = false;

   if(sysType.find("uGMT")!=std::string::npos){
     
     return aL1Cand.type==L1Obj::uGMT && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;
   }
   else if(sysType.find("OMTF")!=std::string::npos){
    
     return aL1Cand.type==L1Obj::OMTF && aL1Cand.q>=12 && aL1Cand.bx==0;
   }   
   return false;

}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillTurnOnCurve( const TLorentzVector & aMuonCand4Vector,
                                  const int & iPtCut,
				                          const std::string & sysType,
				                          const std::string & selType){
  
 //int is important for histo name construction

  int ptCut = GMTHistograms::ptBins[iPtCut];
 /* bool useNanoAOD = (inputType == "nanoAOD");
  if (useNanoAOD) {
     const EventProxyOMTFNANOAOD& myProxy = dynamic_cast<const EventProxyOMTFNANOAOD&>(myEvent);
     const_cast<EventProxyOMTFNANOAOD&>(myProxy).fillnanoL1ObjColl();
     myL1ObjColl = const_cast<EventProxyOMTFNANOAOD&>(myProxy).getL1ObjColl();
 }*/
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
  double deltaR = 0.4;
  L1Obj selectedCand;

   
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);    
   // std::cout<<"the quality and the bunch crossing:  "<< aCand.q << "\t " << aCand.bx << "\t"<< "\n";
    if(!pass) continue;
    double phiValue = aCand.phiValue();
    if(phiValue>M_PI) phiValue-=2*M_PI;
    
    double dEta = std::abs(aMuonCand4Vector.Eta()-aCand.etaValue());
    double dPhi = std::abs(aMuonCand4Vector.Phi()-phiValue);
    if(dPhi>2*M_PI) dPhi=-2*M_PI;
    double delta = sqrt(dEta*dEta + dPhi*dPhi);
    if(delta<deltaR && selectedCand.ptValue()<aCand.ptValue()){
      deltaR = delta;
      selectedCand = aCand;      
    }    
  }

  bool passPtCut = selectedCand.ptValue()>=ptCut && selectedCand.ptValue()>0;
   
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand4Vector.Pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand4Vector.Pt(), passPtCut);

  tmpName = hName+"RecoMuonPtVsL1Pt";
  myHistos_->fill2DHistogram(tmpName, aMuonCand4Vector.Pt(), selectedCand.ptValue());
  
  //Generic eff vs selected variable calculated for muons on plateau
  if(!selType.size() && aMuonCand4Vector.Pt()<ptCut+20) return;
  tmpName = hName+"Eta"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand4Vector.Eta(), passPtCut);

  tmpName = hName+"Phi"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand4Vector.Phi(), passPtCut);
  


}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
 
void GMTAnalyzer::fillRateHisto(const TLorentzVector & aRecoMuon4Vector,
                                const std::string & sysType,
				                        const std::string & selType){

  //Generator level information is not available for the neutrino sample
 

  if(name()=="NU_RATEAnalyzer" && aRecoMuon4Vector.Pt()>0.0) return;

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2D"+sysType+"Rate"+selType;

  L1Obj selectedCand;
  for(auto aCand: myL1Coll){

    bool pass = passQuality(aCand ,sysType, selType);    
    if(pass && selectedCand.ptValue()<aCand.ptValue()) selectedCand = aCand;

  }
  bool pass = selectedCand.ptValue()>=20;
  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon4Vector.Pt(),selectedCand.ptValue());
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon4Vector.Pt(),pass*aRecoMuon4Vector.Eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon4Vector.Pt(),pass*aRecoMuon4Vector.Pt()+(!pass)*(-100));


}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillHistosForRecoMuon( const TLorentzVector & aRecoMuon4Vector){   

 

  bool isGMTAcceptance = fabs(aRecoMuon4Vector.Eta())<2.4;
  if(!isGMTAcceptance) return;

  bool isOMTFAcceptance = fabs(aRecoMuon4Vector.Eta())>0.83 && fabs(aRecoMuon4Vector.Eta())<1.24;
  if(!isOMTFAcceptance) return;
  myHistos_->fill1DHistogram("h1DPtProbe", aRecoMuon4Vector.Pt());
  myHistos_->fill1DHistogram("h1DAbsEtaProbe", std::abs(aRecoMuon4Vector.Eta()));
  
  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurve(aRecoMuon4Vector, iCut, "OMTF", selType);
      fillTurnOnCurve(aRecoMuon4Vector, iCut, "uGMT", selType);
  }

  int iCut = 18;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
    
    if(iType==0) pass = aRecoMuon4Vector.Pt()>ptCut + 20;
    else if(iType==1) pass = aRecoMuon4Vector.Pt()>ptCut && aRecoMuon4Vector.Pt()<(ptCut+5);
    else if(iType==2) pass = aRecoMuon4Vector.Pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(aRecoMuon4Vector, iCut, "OMTF", selType);
    fillTurnOnCurve(aRecoMuon4Vector, iCut, "uGMT", selType);
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
   
  
  bool useNanoAOD = (inputType == "nanoAOD");

if (useNanoAOD) {
    const EventProxyOMTFNANOAOD& myProxy = dynamic_cast<const EventProxyOMTFNANOAOD&>(iEvent);
    const_cast<EventProxyOMTFNANOAOD&>(myProxy).fillnanoL1ObjColl();
    const_cast<EventProxyOMTFNANOAOD&>(myProxy).fillnanoMuonObjColl();
    myMuonObjColl = const_cast<EventProxyOMTFNANOAOD&>(myProxy).getRecoMuonObjColl();
    myL1ObjColl = const_cast<EventProxyOMTFNANOAOD&>(myProxy).getL1ObjColl();
   //cannot perform a static_cast from a const EventProxyBase& to EventProxyOMTFNANOAOD& because the source object (iEvent) is const, while the destination type (EventProxyOMTFNANOAOD&) is non-const.
   //To resolve this issue, you can use const_cast to remove the const qualifier temporarily. 
} 
else {
    const EventProxyOMTF& myProxy = static_cast<const EventProxyOMTF&>(iEvent);
    myEventId = myProxy.getEventId();
    myMuonObjColl = myProxy.getRecoMuonObjColl();
    myL1ObjColl = myProxy.getL1ObjColl();
    myL1PhaseIIObjColl = myProxy.getL1PhaseIIObjColl();
}

/*const std::vector<L1PhaseIIObj> & myL1PhaseIIColl  = myL1PhaseIIObjColl->getL1PhaseIIObjs();

	for(auto aPhaseIICand : myL1PhaseIIColl){
		std::cout<< " phase-II pt, eta and phi & charge : "<< aPhaseIICand.ptValue() <<"\t"<<  aPhaseIICand.etaValue() <<"\t"<< aPhaseIICand.phiValue() <<"\t"<<  aPhaseIICand.chargeValue()<<"\n";
}*/

const std::vector<MuonObj> & myMuonColl = myMuonObjColl->getMuonObjs();
  if(myMuonColl.empty()) return false;
  if(myMuonColl.size() < 2 )return false;

   MuonObj aTagCand =  myMuonColl.at(0);
  bool tagPass = aTagCand.pt()>10 && aTagCand.matchedisohlt();
  if(!tagPass) return true;
 
  

  tagFourVector.SetPtEtaPhiM(aTagCand.pt(), aTagCand.eta(), aTagCand.phi(), nominalMuonMass);
  myHistos_->fill1DHistogram("h1DPtTag", tagFourVector.Pt());
  myHistos_->fill1DHistogram("h1DAbsEtaTag", std::abs(tagFourVector.Eta()));

  MuonObj aProbeCand;
  double m_Z = 91.1876;
  double deltaM_Z = 20;
  double tmpDelta = 2*deltaM_Z;
  for (auto aMuonCand: myMuonColl){   
      randomMuonLeg.SetPtEtaPhiM(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), nominalMuonMass);
      fillRateHisto(randomMuonLeg, "OMTF","Tot");
      fillRateHisto(randomMuonLeg, "OMTF","Tot");
      tmpDelta = std::abs((tagFourVector+randomMuonLeg).M()-m_Z);
      if(aMuonCand.tightID() && tmpDelta<deltaM_Z){
      deltaM_Z = tmpDelta;
      aProbeCand = aMuonCand;   
    }
  }
  bool isOMTFAcceptanceP = fabs(aProbeCand.eta())>0.83 && fabs(aProbeCand.eta())<1.24;
  if(!isOMTFAcceptanceP) return false;
  if(aProbeCand.pt()<1) return false;// Select only hard scattering proccesses 

  probeFourVector.SetPtEtaPhiM(aProbeCand.pt(), aProbeCand.eta(), aProbeCand.phi(), nominalMuonMass);
  myHistos_->fill1DHistogram("h1DDiMuonMassTagProbe",(tagFourVector+probeFourVector).M());   
  fillHistosForRecoMuon(probeFourVector);

  double deltaRCut = 0.6;
  for (auto aMuonCand: myMuonColl){   
      TLorentzVector puMuonLeg;
      puMuonLeg.SetPtEtaPhiM(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), nominalMuonMass);
      if(aMuonCand.tightID() 
               && customDeltaR(puMuonLeg,tagFourVector) > deltaRCut 
               && customDeltaR(puMuonLeg,probeFourVector) > deltaRCut){
        fillHistosForRecoMuon(puMuonLeg);
      }
    }
    
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
