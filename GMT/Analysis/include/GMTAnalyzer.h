#ifndef GMTAnalyzer_H
#define GMTAnalyzer_H

#include <string>
#include <map>

#include "GMTHistograms.h"
#include "Analyzer.h"

#include "EventObj.h"
#include "GenObjColl.h"
#include "L1ObjColl.h"
#include "L1PhaseIIObjColl.h"
#include "L1PhaseIIObj.h"
#include "MuonObjColl.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
using namespace TMath;
using namespace std;
class TriggerHistograms;
class OMTFHit;

class GMTAnalyzer:public Analyzer{

 public:

  GMTAnalyzer(const std::string & aName);

  ~GMTAnalyzer();

  void initialize(TDirectory* aDir, pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);
  void finalize();

  Analyzer* clone() const;
  double zmass ;
  double nominalMuonMass = 0.1056583;
  double tpdeltaR = 0.4;
  TLorentzVector theMuonLegPositive;
  TLorentzVector theMuonLegNegative;
  TLorentzVector theZResonance;
  TLorentzVector tagFourVector;
  TLorentzVector probeFourVector;
  TLorentzVector randomMuonLeg;
  TLorentzVector puMuonLeg;
  TVector3 tagthreeVector;
  TVector3 probethreeVector;
  TVector3 puthreeVector;
  TVector3 randomthreeVector;
  TVector3 level1threeVector;
  void setHistos(GMTHistograms *histos) { myHistos_ = histos;};
  void parseProcessName(); 
private:
  
  void fillHistosForRecoMuon( const TVector3 & aRecoMuon3Vector, const std::string & ObjectTag);  
  

  void fillTurnOnCurve(  const TVector3 & aMuonCand3Vector,
                      const int & ptCut, const std::string & sysType,
		                  const std::string & selType);

  void fillRateHisto(const TVector3 & aRecoMuon3Vector,
                    const std::string & sysType,
		                const std::string & selType);
  double customDeltaR(TLorentzVector T1, TLorentzVector T2);
  bool passQuality(const L1Obj & aL1Cand,
		              const std::string & sysType,
		              const std::string & selType = "");
  double zResonance(const MuonObj  aRecoMuon);
  double detaTagAndProbe(const MuonObj  aRecoMuon); 
  std::string inputType;
  ///Histograms for this analysis
  GMTHistograms *myHistos_;
  
  
  const EventObj  *myEventId;
  const MuonObjColl  *myMuonObjColl;
  const L1ObjColl  *myL1ObjColl;
  const L1PhaseIIObjColl  *myL1PhaseIIObjColl; 
};

#endif
