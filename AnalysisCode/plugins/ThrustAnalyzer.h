#ifndef ThrustAnalyzer_H
#define ThrustAnalyzer_H


// Header file for Thrust Analyzer. Macro necessary to calculate the Jet Thrust, thrust major, thrust minor
// Built for Fast analysis in Run2 data, CMSSW_7_5_3_patch1, 
// Auther: Raghav Kunnawalkam Elayavalli
//         Oct 13th 2015
//         Rutgers University. email: raghav.k.e at CERN dot CH
//
// The main logic for the code is written by Jennifer Coulter, Rutgers Undergraduate student during the summer of 2015
// https://github.com/jcoulter120/SummerRutgers15/blob/master/thrust_HiForest.C with some additions and changes by Raghav. 
//


#include <memory>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandidateWithRef.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// include the basic jet for the PuPF jets. 
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
// include the pf candidates 
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
// include the voronoi subtraction
#include "DataFormats/HeavyIonEvent/interface/VoronoiBackground.h"
#include "RecoHI/HiJetAlgos/interface/UEParameters.h"
// include the centrality variables
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

// root headers

#include <TH1.h>
#include <TH2.h>
#include "TROOT.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace edm;
using namespace std;

class ThrustAnalyzer : public edm::EDAnalyzer{
 public:
  explicit ThrustAnalyzer(const edm::ParameterSet&);
  ~ThrustAnalyzer() {}

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  // necessary jet collections
  edm::InputTag mInputCollection;
  std::string JetType;
  std::string UEAlgo;
  int radius;
  edm::InputTag tagVtx;
  // edm::InputTag mInputPFCandCollection;
  double mRThreshold;
  double mRecoJetPtThreshold;
  double mReco_SubJetPtThreshold;
  //std::string JetCorrectionService;
  /* std::string mhltPath; */
  /* edm::EDGetTokenT<edm::TriggerResults> triggerBits_; */
  /* edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_; */
  /* edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_; */
  /* edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent_; */
  
  //edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;
  //int triggerBit;    
  
  bool isCaloJet;
  bool isJPTJet;
  bool isPFJet;
  edm::EDGetTokenT<reco::CaloJetCollection> caloJetsToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetsToken_;
  edm::EDGetTokenT<reco::BasicJetCollection> basicJetsToken_;
  edm::EDGetTokenT<reco::JPTJetCollection> jptJetsToken_;
  /* edm::EDGetTokenT<reco::PFCandidateCollection> pfCandToken_;  */

  /* edm::InputTag centralityTag_; */
  /* edm::EDGetTokenT<reco::Centrality> centralityToken; */
  /* edm::Handle<reco::Centrality> centrality_; */
  /* edm::InputTag centralityBinTag_; */
  /* edm::EDGetTokenT<int> centralityBinToken; */
  /* edm::Handle<int>centralityBin_; */

  int event, run, lumi, hiBin;
  double vz, vx, vy;
  
  // variables for the output root files
  edm::Service<TFileService> fout;
  map< string, TH1D* > histos1D;  
  map< string, TH2D* > histos2D;

  //Function to normalize a vector
  static TVector3 Norm(TVector3 v){
    if ( (v(0) == 0) && (v(1) == 0) && (v(2) == 0)) return v; 
    Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
    v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
    return v; 
  }//end normalize

  static float deltaPhi(float phi1, float phi2){
    float dphi = fabs(phi1 - phi2);
    if(dphi > M_PI )dphi -= 2*M_PI;
    return dphi;
  }
  
};

//plane class
class Dplane{
 public:
  TVector3 v1, v2, proj, u1, u2;
  Double_t scalar1, scalar2, mag1, mag2; 
  Dplane(TVector3);

  //returns a projection onto the 2D plane 
  TVector3 Projection(TVector3 jaxis){
    //Find the projection of a jet onto this subspace

    if(v1.Mag() == 0) { scalar1 = 0; }   else { scalar1 = jaxis.Dot(v1)/(v1.Dot(v1)); } 
    if(v2.Mag() == 0) { scalar2 = 0; }   else { scalar2 = jaxis.Dot(v2)/(v2.Dot(v2)); } 
    v1 = scalar1*v1;
    v2 = scalar2*v2;
    proj(0) = v1(0) + v2(0);
    proj(1) = v1(1) + v2(1);
    proj(2) = v1(2) + v2(2); 
    
    return proj;
  }//end of projection
};
  
//plane class constructor
Dplane::Dplane(TVector3 nT){
  //Use TVector3 to find an orthogonal vector and a second vector orthogonal to the first and nT
  v1 = nT.Orthogonal();  v2 = nT.Cross(v1);

  //Normalize, checking for 0 length axes
  if ((v1(0) == 0) && (v1(1) == 0) && (v1(2) == 0)){  v1(0) = 0;    v1(1) = 0;    v1(2) = 0; }
  else { mag1 = v1.Mag();   v1(0) = v1(0)/mag1;    v1(1) = v1(1)/mag1;    v1(2) = v1(2)/mag1; } 
  if ((v2(0) == 0) && (v2(1) == 0) && (v2(2) == 0)){  v2(0) = 0;    v2(1) = 0;    v2(2) = 0; } 
  else { mag2 = v2.Mag();   v2(0) = v2(0)/mag2;    v2(1) = v2(1)/mag2;    v2(2) = v2(2)/mag2; }
}//end plane constructor


#endif
