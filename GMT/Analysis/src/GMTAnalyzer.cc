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
void GMTAnalyzer::fillTurnOnCurve( const TVector3 & instantiatedVector,
                                  const int &  iPtCut, const std::string & sysType,
				                       const std::string & selType){
  
 //int is important for histo name construction

  int ptCut = GMTHistograms::ptBins[iPtCut];
  aRecoMuon3Vector = instantiatedVector;
  
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
  TVector3 selectedL13Vector;
  double deltaR = 0.2;
  bool useNanoAOD = (inputType == "nanoAOD");
  if (useNanoAOD) { 
 for (const auto& obj : aL1Object3VectorCollection) {
    //std::cout<< " pt, eta , phi :"<< obj.Pt() <<"\t"<< obj.Eta() << "\t"<< obj.Phi()<< "\n";
    double phiValue = obj.Phi();
    if (phiValue > M_PI) phiValue -= 2 * M_PI;

    double dEta = std::abs(aRecoMuon3Vector.Eta() - obj.Eta());
    double dPhi = std::abs(aRecoMuon3Vector.Phi() - phiValue);
    if (dPhi > 2 * M_PI) dPhi -= 2 * M_PI;
    
    double delta = sqrt(dEta * dEta + dPhi * dPhi);
    if (delta < 3.0 && selectedL13Vector.Pt() < obj.Pt() ) {
        deltaR = delta;
        selectedL13Vector.SetPtEtaPhi(obj.Pt(), obj.Eta(), obj.Phi());
    }
  }
}
 
else{
  L1Obj selectedCand;
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);    
    if(!pass) continue;
    double phiValue = aCand.phiValue();
    if(phiValue>M_PI) phiValue-=2*M_PI;
    
    double dEta = std::abs(aRecoMuon3Vector.Eta()-aCand.etaValue());
    double dPhi = std::abs(aRecoMuon3Vector.Phi()-phiValue);
    if(dPhi>2*M_PI) dPhi=-2*M_PI;
    double delta = sqrt(dEta*dEta + dPhi*dPhi);
    if(delta<deltaR && selectedCand.ptValue()<aCand.ptValue()){
      deltaR = delta;
      selectedCand = aCand;
      selectedL13Vector.SetPtEtaPhi(aCand.ptValue(), aCand.etaValue(), aCand.phiValue());      
    }    
  }
}
  
  bool passPtCut = selectedL13Vector.Pt()>=ptCut && selectedL13Vector.Pt()>0;
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon3Vector.Pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon3Vector.Pt(), passPtCut);

  tmpName = hName+"RecoMuonPtVsL1Pt";
  myHistos_->fill2DHistogram(tmpName, aRecoMuon3Vector.Pt(), selectedL13Vector.Pt());
  
  //Generic eff vs selected variable calculated for muons on plateau
  if(!selType.size() && aRecoMuon3Vector.Pt()<ptCut+20) return;
  tmpName = hName+"Eta"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon3Vector.Eta(), passPtCut);

  tmpName = hName+"Phi"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon3Vector.Phi(), passPtCut);
  


}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
 
void GMTAnalyzer::fillRateHisto(const TVector3 & instantiatedVector,
                                const std::string & sysType,
				                        const std::string & selType){
  
  bool useNanoAOD = (inputType == "nanoAOD"); 
  std::string hName = "h2D"+sysType+"Rate"+selType; 
  aRecoMuon3Vector = instantiatedVector;
  TVector3 selectedL13Vector;
  //Generator level information is not available for the neutrino sample
  if(name()=="NU_RATEAnalyzer" && aRecoMuon3Vector.Pt()>0.0) return;

 if (useNanoAOD) {
  for (const auto& obj : aL1Object3VectorCollection) {
    if (selectedL13Vector.Pt() < obj.Pt())selectedL13Vector = obj;
  }
 }
 else{
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  L1Obj selectedCand;
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);    
    if(pass && selectedCand.ptValue()<aCand.ptValue()) selectedCand = aCand;
    selectedL13Vector.SetPtEtaPhi(selectedCand.ptValue(), selectedCand.etaValue(), selectedCand.phiValue());   
  }
}
  //double lhcFrequency = 11245.6;  //Hz
  //double collidingBunches = 2544; 
  double aFactor = 1.0;// lhcFrequency * collidingBunches; 
  bool pass = selectedL13Vector.Pt() >= 20;
  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon3Vector.Pt(),selectedL13Vector.Pt()*aFactor);
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon3Vector.Pt(),pass*aRecoMuon3Vector.Eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon3Vector.Pt(),pass*aRecoMuon3Vector.Pt()+(!pass)*(-100));


}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillHistosForObjectVectors( const TVector3 & instantiatedVector){   
          
          

  aRecoMuon3Vector = instantiatedVector; 
  bool isGMTAcceptance = fabs(aRecoMuon3Vector.Eta())<2.4;
  if(!isGMTAcceptance) return;

  bool isOMTFAcceptance = fabs(aRecoMuon3Vector.Eta())>0.83 && fabs(aRecoMuon3Vector.Eta())<1.24;
  if(!isOMTFAcceptance) return;
  myHistos_->fill1DHistogram("h1DPtProbe", aRecoMuon3Vector.Pt());
  myHistos_->fill1DHistogram("h1DAbsEtaProbe", std::abs(aRecoMuon3Vector.Eta()));
  
  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurve(aRecoMuon3Vector, iCut, "OMTF", selType);
      fillTurnOnCurve(aRecoMuon3Vector, iCut, "uGMT", selType);
  }

  int iCut = 18;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
    
    if(iType==0) pass = aRecoMuon3Vector.Pt()>ptCut + 20;
    else if(iType==1) pass = aRecoMuon3Vector.Pt()>ptCut && aRecoMuon3Vector.Pt()<(ptCut+3);
    else if(iType==2) pass = aRecoMuon3Vector.Pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(aRecoMuon3Vector, iCut, "OMTF", selType);
    fillTurnOnCurve(aRecoMuon3Vector, iCut, "uGMT", selType);
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
   //std::cout << " the particle flow and the tracker isolation value " << aTagCand.isTrackIsolated()<< "\t"<< aTagCand.isParticleFlowIsolated()<< "\n";
   bool tagPass = aTagCand.pt()>27 && aTagCand.matchedisohlt() && aTagCand.tightID() && aTagCand.isParticleFlowIsolated();
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
      randomthreeVector.SetPtEtaPhi(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi());
      fillRateHisto(randomthreeVector, "OMTF","Tot");
      fillRateHisto(randomthreeVector, "uGMT","Tot");
      tmpDelta = std::abs((tagFourVector+randomMuonLeg).M()-m_Z);
      
      if(aMuonCand.tightID() && aMuonCand.isParticleFlowIsolated()&& tmpDelta<deltaM_Z  && customDeltaR(randomMuonLeg, tagFourVector) > 0.4){
      deltaM_Z = tmpDelta;
      aProbeCand = aMuonCand;   
    }
  }
  if(aProbeCand.pt()<1) return false;// Select only hard scattering proccesses 

  probeFourVector.SetPtEtaPhiM(aProbeCand.pt(), aProbeCand.eta(), aProbeCand.phi(), nominalMuonMass);
  probethreeVector.SetPtEtaPhi(aProbeCand.pt(), aProbeCand.eta(), aProbeCand.phi());
  myHistos_->fill1DHistogram("h1DDiMuonMassTagProbe",(tagFourVector+probeFourVector).M());   
  fillHistosForObjectVectors(probethreeVector);


  double deltaRCut = 0.6;
  for (auto aMuonCand: myMuonColl){   
      puMuonLeg.SetPtEtaPhiM(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), nominalMuonMass);
      puthreeVector.SetPtEtaPhi(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi());
      if(aMuonCand.tightID() 
               && customDeltaR(puMuonLeg,tagFourVector) > deltaRCut 
               && customDeltaR(puMuonLeg,probeFourVector) > deltaRCut){
        fillHistosForObjectVectors(puthreeVector);
      }
    }
  

  ////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// LEVELONE-OBJECT //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////


if (useNanoAOD) {
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  aL1Object3VectorCollection.clear();
  for(auto levelOneCand: myL1Coll){
               bool localnanoPass = false;   
               localnanoPass = levelOneCand.q == 12 && levelOneCand.bx == 0;
               if(!localnanoPass) continue; 
               double epsilon = 1e-6;
               double eta = levelOneCand.etaValue();
               if (std::isnan(eta) || (fabs(eta) < epsilon)) {continue;}
               aL1Object3Vector.SetPtEtaPhi(levelOneCand.ptValue(), levelOneCand.etaValue(), levelOneCand.nanoPhiValue());
               aL1Object3VectorCollection.push_back(aL1Object3Vector);

/* 
     std::cout << "type of the detector     : " << levelOneCand.type           << "\n"
               << "Quality of the processor : " << levelOneCand.q              << "\n";
	       << "Tranverse momentum       : " << levelOneCand.ptValue()      << "\n"
               << "Pseudorapidity           : " << levelOneCand.etaValue()     << "\n"
               << "Azimuthal angle          : " << levelOneCand.nanoPhiValue() << "\n"
               << "Bunch crossing           : " << levelOneCand.bx             << "\n"
               << "-----------------------------------------------------------"<< "\n";
*/ 
               }
 }


  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
