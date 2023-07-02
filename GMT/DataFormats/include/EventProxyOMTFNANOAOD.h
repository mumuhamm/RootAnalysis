#ifndef EVENTPROXYOMTFNANOAOD_h
#define EVENTPROXYOMTFNANOAOD_h

#include <string>
#include <typeinfo>
#include <vector>
#include "boost/shared_ptr.hpp"

#include "EventProxyBase.h"
#include "EventObj.h"
#include "GenObjColl.h"
#include "L1ObjColl.h"
#include "L1PhaseIIObjColl.h"
#include "MuonObjColl.h"
#include "TBranch.h"
#include "L1Obj.h"

struct OMTFHit{
  
public:
  
  int iPhi, iLayer, iHit, iQuality;
  
};

std::ostream& operator<< (std::ostream& stream, const OMTFHit& aHit); 


   class EventProxyOMTFNANOAOD: public EventProxyBase{

   public:

      EventProxyOMTFNANOAOD();
      virtual ~EventProxyOMTFNANOAOD();

      void init(std::vector<std::string> const& iFileNames);

      virtual EventProxyBase* clone() const;

     const EventObj  *getEventId() const {return myEvent;};
     const GenObjColl *getGenObjColl() const {return myGenObjColl;};
     const L1ObjColl  *getL1ObjColl() const { return myL1ObjColl;};
     const L1PhaseIIObjColl  *getL1PhaseIIObjColl() const { return myL1PhaseIIObjColl;};
     const MuonObjColl *getRecoMuonObjColl() const { return myMuonObjColl;}; 
     std::vector<OMTFHit> getHits() const;


   
   L1Obj aL1Obj;
   MuonObj aMuonObj;
   private:
   //Declaration of fill function : 
   void fillnanoL1ObjColl();
   void fillnanoMuonObjColl();
   //Declaration of initiations of branches: 

     Int_t           nL1Mu;
     Float_t         L1Mu_eta[128];
     Float_t         L1Mu_phi[128];
     Float_t         L1Mu_pt[128];
     Short_t         L1Mu_hwCharge[128];
     Short_t         L1Mu_hwDXY[128];
     Short_t         L1Mu_bx[128];
     Int_t           L1Mu_hwQual[128];
     Float_t         L1Mu_etaAtVtx[128]; 
     Float_t         L1Mu_phiAtVtx[128];

     Int_t           nMuon;
     Int_t           Muon_charge[128];
     Int_t           Muon_pdgId[128];
     Float_t         Muon_eta[128];
     Float_t         Muon_phi[128];
     Float_t         Muon_pt[128];
     Bool_t          Muon_tightId[128];
     UChar_t         Muon_pfIsoId[128];
     Bool_t          Muon_mediumId[128];      
     Bool_t          Muon_softId[128];
     Bool_t          HLT_IsoMu20;
     Bool_t          HLT_IsoMu24;
     Bool_t          HLT_IsoMu27;


      //Declaration of leaf types
     EventObj       *myEvent;
     GenObjColl     *myGenObjColl;
     L1ObjColl      *myL1ObjColl;
     L1PhaseIIObjColl *myL1PhaseIIObjColl;
     MuonObjColl *myMuonObjColl;
     std::vector<int> *hits;
     std::vector<int> *hitsQuality;
     

};
#endif
