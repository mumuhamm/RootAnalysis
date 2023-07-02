nclude "EventProxyOMTFNANOAOD.h"
EventProxyOMTFNANOAOD::EventProxyOMTFNANOAOD(){
     myL1ObjColl = new L1ObjColl();
     myMuonObjColl = new MuonObjColl();
}
EventProxyOMTFNANOAOD::~EventProxyOMTFNANOAOD(){
    delete myL1ObjColl;
    delete myMuonObjColl;
}

EventProxyBase* EventProxyOMTFNANOAOD::clone() const{

  return new EventProxyOMTFNANOAOD();
  
}

void EventProxyOMTFNANOAOD::init(std::vector<std::string> const& iFileNames){

  treeName_ = "Events";   
  EventProxyBase::init(iFileNames);
  fChain->SetBranchStatus("*",0);
  fChain->SetMakeClass(0);
}
 void  EventProxyOMTFNANOAOD::fillnanoL1ObjColl() const{
    fChain->SetBranchStatus("nL1Mu", 1);
    fChain->SetBranchStatus("L1Mu_hwCharge", 1);
    fChain->SetBranchStatus("L1Mu_hwDXY", 1);
    fChain->SetBranchStatus("L1Mu_bx", 1);
    fChain->SetBranchStatus("L1Mu_hwQual", 1);
    fChain->SetBranchStatus("L1Mu_eta", 1);
    fChain->SetBranchStatus("L1Mu_etaAtVtx", 1);
    fChain->SetBranchStatus("L1Mu_phi", 1);
    fChain->SetBranchStatus("L1Mu_phiAtVtx", 1);
    fChain->SetBranchStatus("L1Mu_pt", 1);
          for( Int_t i =0; i< nL1Mu; ++i){
              aL1Obj.eta  = L1Mu_eta[i];
	      aL1Obj.phi  = L1Mu_phi[i];
              aL1Obj.pt   = L1Mu_pt[i];
              aL1Obj.charge  = L1Mu_hwCharge[i];
              aL1Obj.bx        = L1Mu_bx[i];
              aL1Obj.q    = L1Mu_hwQual[i];
              myL1ObjColl->push_back(aL1Obj, false, 0.0);
             
            }
}
void EventProxyOMTFNANOAOD::fillnanoMuonObjColl() const{
     fChain->SetBranchStatus("nMuon", 1);
    fChain->SetBranchStatus("Muon_mediumId", 1);
    fChain->SetBranchStatus("Muon_pfIsoId", 1);
    fChain->SetBranchStatus("Muon_softId", 1);
    fChain->SetBranchStatus("Muon_tightId", 1);
    fChain->SetBranchStatus("Muon_charge", 1);
    fChain->SetBranchStatus("Muon_pdgId", 1);
    fChain->SetBranchStatus("Muon_eta", 1);
    fChain->SetBranchStatus("Muon_phi", 1);
    fChain->SetBranchStatus("Muon_pt", 1);
    fChain->SetBranchStatus("HLT_IsoMu20", 1);
    fChain->SetBranchStatus("HLT_IsoMu24", 1);
    fChain->SetBranchStatus("HLT_IsoMu27", 1);
         for (Int_t i = 0; i < nMuon; ++i) {
            aMuonObj.theCharge = Muon_charge[i];
            aMuonObj.thePt = Muon_pt[i];
            aMuonObj.theEta = Muon_eta[i];
            aMuonObj.thePhi = Muon_phi[i];
            aMuonObj.isTight = Muon_tightId[i];
            aMuonObj.isMedium = Muon_mediumId[i];
            aMuonObj.isMatchedIsoHlt = (HLT_IsoMu20 || HLT_IsoMu24 || HLT_IsoMu27);
            myMuonObjColl->push_back(aMuonObj, false, 0.0);
            
        }
}
std::vector<OMTFHit> EventProxyOMTFNANOAOD::getHits() const {
  std::vector<OMTFHit> hits;
  // TODO: Implement the logic to populate the `hits` vector
  return hits;
}
