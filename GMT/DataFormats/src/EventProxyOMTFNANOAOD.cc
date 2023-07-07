#include "EventProxyOMTFNANOAOD.h"



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
  fChain->SetBranchStatus("*",1);
  fChain->SetMakeClass(0);

   TBranch        *b_nL1Mu;   	
   TBranch        *b_L1Mu_eta;
   TBranch        *b_L1Mu_phi; 
   TBranch        *b_L1Mu_pt;   
   TBranch        *b_L1Mu_hwCharge;   
   TBranch        *b_L1Mu_hwDXY;   
   TBranch        *b_L1Mu_bx;   
   TBranch        *b_L1Mu_hwQual;   
   TBranch        *b_L1Mu_etaAtVtx;   
   TBranch        *b_L1Mu_phiAtVtx; 
  

   TBranch        *b_nMuon;   
   TBranch        *b_Muon_mediumId;   
   TBranch        *b_Muon_pfIsoId;   
   TBranch        *b_Muon_softId;   
   TBranch        *b_Muon_tightId;   
   TBranch        *b_Muon_charge;   
   TBranch        *b_Muon_pdgId;   
   TBranch        *b_Muon_eta;   
   TBranch        *b_Muon_phi;   
   TBranch        *b_Muon_pt;   
   TBranch        *b_HLT_IsoMu20;   
   TBranch        *b_HLT_IsoMu24;   
   TBranch        *b_HLT_IsoMu27;   
   

    fChain->SetBranchAddress("nL1Mu", &nL1Mu, &b_nL1Mu);
   fChain->SetBranchAddress("L1Mu_hwCharge", L1Mu_hwCharge, &b_L1Mu_hwCharge);
   fChain->SetBranchAddress("L1Mu_hwDXY", L1Mu_hwDXY, &b_L1Mu_hwDXY);
   fChain->SetBranchAddress("L1Mu_bx", L1Mu_bx, &b_L1Mu_bx);
   fChain->SetBranchAddress("L1Mu_hwQual", L1Mu_hwQual, &b_L1Mu_hwQual);
   fChain->SetBranchAddress("L1Mu_eta", L1Mu_eta, &b_L1Mu_eta);
   fChain->SetBranchAddress("L1Mu_etaAtVtx", L1Mu_etaAtVtx, &b_L1Mu_etaAtVtx);
   fChain->SetBranchAddress("L1Mu_phi", L1Mu_phi, &b_L1Mu_phi);
   fChain->SetBranchAddress("L1Mu_phiAtVtx", L1Mu_phiAtVtx, &b_L1Mu_phiAtVtx);
   fChain->SetBranchAddress("L1Mu_pt", L1Mu_pt, &b_L1Mu_pt);    

  fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
  fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
  fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
  fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
  fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
  fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
  fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
  fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
  fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
  fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt); 
  fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
  fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
  fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);



}
void  EventProxyOMTFNANOAOD::fillnanoL1ObjColl() {
          for( Int_t i =0; i< nL1Mu; ++i){
              //std::cout<< " the iteration: "<< i << "\t and the eta: "<< L1Mu_eta[i]<< "\n";
              aL1Obj.eta  = L1Mu_eta[i];
              //std::cout<< " now the eta from the asigned object : "<< aL1Obj.eta<< "\n"; 
	      aL1Obj.phi  = L1Mu_phi[i];
              aL1Obj.pt   = L1Mu_pt[i];
              aL1Obj.charge  = L1Mu_hwCharge[i];
              aL1Obj.bx        = L1Mu_bx[i];
              aL1Obj.q    = L1Mu_hwQual[i];
              myL1ObjColl->push_back(aL1Obj, false, 0.0);
             
            }
}
void EventProxyOMTFNANOAOD::fillnanoMuonObjColl()  {
         for (Int_t i = 0; i < nMuon; ++i) {
            aMuonObj.setCharge(Muon_charge[i]);
            aMuonObj.setPt(Muon_pt[i]);
            aMuonObj.setEta(Muon_eta[i]);
            aMuonObj.setPhi(Muon_phi[i]);
            aMuonObj.setTight(Muon_tightId[i]);
            aMuonObj.setMedium(Muon_mediumId[i]);
            aMuonObj.setMatchedIsoHlt(HLT_IsoMu20 || HLT_IsoMu24 || HLT_IsoMu27);
            myMuonObjColl->addMuonObj(aMuonObj);
            
        }

}
std::vector<OMTFHitNano> EventProxyOMTFNANOAOD::getHits() const {
  std::vector<OMTFHitNano> hits;
  // TODO: Implement the logic to populate the `hits` vector
  return hits;
}
