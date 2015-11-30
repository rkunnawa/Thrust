// -*- C++ -*-
//
// Package:    ThrustAnalyzer
// Class:      ThrustAnalyzer
// 
/**\class ThrustAnalyzer ThrustAnalyzer.cc Thrust/AnalysisCode/src/ThrustAnalyzer.cc

 Description: 

 Implementation:
     
*/
//
// Original Author:  Raghav Kunnawalkam Elayavalli
//                   Rutgers University, NJ 
//         Created:  Wednesday Oct 7 19:23:23 EST 2015
// $Id$
//
//

#include "Thrust/AnalysisCode/plugins/ThrustAnalyzer.h"

using namespace std;
using namespace reco;
using namespace edm;

// declare the constructors 

ThrustAnalyzer::ThrustAnalyzer(const edm::ParameterSet& iConfig) :
  mInputCollection               (iConfig.getParameter<edm::InputTag>       ("jet")),
  JetType                        (iConfig.getUntrackedParameter<std::string>("JetType")),
  UEAlgo                         (iConfig.getUntrackedParameter<std::string>("UEAlgo")),
  radius                         (iConfig.getUntrackedParameter<int>        ("radius")),
  tagVtx                         (iConfig.getParameter<edm::InputTag>       ("vertices")),
  // mInputPFCandCollection         (iConfig.getParameter<edm::InputTag>       ("PFcands")),
  mRThreshold                    (iConfig.getParameter<double>              ("RThreshold")),
  mRecoJetPtThreshold            (iConfig.getParameter<double>              ("mRecoJetPtThreshold")),
  mReco_SubJetPtThreshold        (iConfig.getParameter<double>              ("mReco_SubJetPtThreshold"))
  //JetCorrectionService           (iConfig.getParameter<std::string>         ("JetCorrections")),
  // mhltPath                       (iConfig.getParameter<std::string>         ("hltpath")),
  // triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  // triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  // triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  // triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("hltTrigger"))),
{
  
  std::string inputCollectionLabel(mInputCollection.label());
  
  isCaloJet = (std::string("calo")==JetType);
  isJPTJet  = (std::string("jpt") ==JetType);
  isPFJet   = (std::string("pf")  ==JetType);
  
  if (isCaloJet) caloJetsToken_  = consumes<reco::CaloJetCollection>(mInputCollection);
  if (isJPTJet)  jptJetsToken_   = consumes<reco::JPTJetCollection>(mInputCollection);
  if (isPFJet)   {
    if(std::string("Pu")==UEAlgo) basicJetsToken_    = consumes<reco::BasicJetCollection>(mInputCollection);
    if(std::string("Vs")==UEAlgo) pfJetsToken_    = consumes<reco::PFJetCollection>(mInputCollection);
  }

  // pfCandToken_ = consumes<reco::PFCandidateCollection>(mInputPFCandCollection);
  // centralityTag_ = iConfig.getParameter<InputTag>("centralitycollection");
  // centralityToken = consumes<reco::Centrality>(centralityTag_);

  // centralityBinTag_ = (iConfig.getParameter<edm::InputTag> ("centralitybincollection"));
  // centralityBinToken = consumes<int>(centralityBinTag_);

}


void ThrustAnalyzer::beginJob() {

  TFileDirectory jet = fout->mkdir("ak" + UEAlgo + Form("%d",radius) + JetType + "_Thrust");

  // setup the histograms or tree for the next step.
  histos1D["hThrust"] = jet.make< TH1D >("hThrust", ";#tau;counts", 50, 0, 1);
  histos1D["hThrust_min"] = jet.make< TH1D >("hThrust_min", ";#tau_{minor};counts", 50, 0, 1);
  histos1D["hThrust_maj"] = jet.make< TH1D >("hThrust_maj", ";#tau_{major};counts", 50, 0, 1);
  histos1D["hpT"] = jet.make< TH1D >("hpT", ";Jet p_{T} (GeV/c);counts", 500, 0, 500);
  histos1D["heta"] = jet.make< TH1D >("heta", ";#eta;counts", 60, -3, 3);
  histos1D["hphi"] = jet.make< TH1D >("hphi", ";#phi;counts", 60, -3.2, 3.2);
  histos2D["hThrust_vs_Aj"] = jet.make< TH2D >("hThrust_vs_Aj", ";A_{j};#tau", 50, 0, 1, 50, 0, 1);
  // need to add thrust vs Aj once centrality bin variable goes into the code  
}

void ThrustAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //Get run, event, centrality
  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();

  bool debug = false;
  
  edm::Handle< vector<reco::Vertex> > vertex;
  iEvent.getByLabel(tagVtx,vertex);
  if(vertex->size() > 0){
    vx = vertex->begin()->x();
    vy = vertex->begin()->y();
    vz = vertex->begin()->z();
  } else {
    vx = -1;
    vy = -1;
    vz = -1;
  }

  if(fabs(vz) > 15) return;

  // get the centrality 
  // edm::Handle<reco::Centrality> cent;
  // iEvent.getByToken(centralityToken, cent);  //_centralitytag comes from the cfg
  // mHF->Fill(cent->EtHFtowerSum());
  // Float_t HF_energy = cent->EtHFtowerSum();  
  // edm::Handle<int> cbin;
  // iEvent.getByToken(centralityBinToken, cbin);
  
  // if(!cent.isValid()) return;
  // int hibin = -999;
  // if(cent.isValid())
  //   hibin = *cbin;

  // pf candidates
  // iEvent.getByToken(pfCandToken_, pfCandidates);
  // const reco::PFCandidateCollection *pfCandidateColl = pfCandidates.product();

  // for(unsigned icand=0;icand<pfCandidateColl->size(); icand++){
  //   const reco::PFCandidate pfCandidate = pfCandidateColl->at(icand);
  //   reco::CandidateViewRef ref(pfcandidates_,icand);
  //   if(pfCandidate.pt() < 5) continue;
  // }// pf candidiate loop


  std::vector<Jet> recoJets;
  recoJets.clear();
  
  edm::Handle<CaloJetCollection>  caloJets;
  edm::Handle<JPTJetCollection>   jptJets;
  edm::Handle<PFJetCollection>    pfJets;
  edm::Handle<BasicJetCollection> basicJets;

  if (isCaloJet) iEvent.getByToken(caloJetsToken_, caloJets);
  if (isJPTJet)  iEvent.getByToken(jptJetsToken_, jptJets);
  if (isPFJet) {  
    if(std::string("Pu")==UEAlgo) iEvent.getByToken(basicJetsToken_, basicJets);
    if(std::string("Vs")==UEAlgo) iEvent.getByToken(pfJetsToken_, pfJets);
  }


  if (isCaloJet && !caloJets.isValid()) {
    return;
  }
  if (isJPTJet  && !jptJets.isValid()) {
    return;
  }
  if (isPFJet){
    if(std::string("Pu")==UEAlgo){if(!basicJets.isValid())   return;}
    if(std::string("Vs")==UEAlgo){if(!pfJets.isValid())   return;}
  }

  if (isCaloJet){
    for (unsigned ijet=0; ijet<caloJets->size(); ijet++) {
      recoJets.push_back((*caloJets)[ijet]);
    } 
  }

  if (isJPTJet){
    for (unsigned ijet=0; ijet<jptJets->size(); ijet++) 
      recoJets.push_back((*jptJets)[ijet]);
  }

  if (isPFJet) {
    if(std::string("Pu")==UEAlgo){
      for (unsigned ijet=0; ijet<basicJets->size();ijet++) {
	recoJets.push_back((*basicJets)[ijet]);
      }
    }
    if(std::string("Vs")==UEAlgo){
      for (unsigned ijet=0; ijet<pfJets->size(); ijet++){
	recoJets.push_back((*pfJets)[ijet]);
      }
    }
  }

  // bool isDijetEvent = false;
  // //apply dijet selection
  // if(recoJets[0].pt() > mRecoJetPtThreshold && recoJets[1].pt() > mReco_SubJetPtThreshold) 
  //   isDijetEvent = true;

  // if(!isDijetEvent) return;

  Float_t dot = 0;
  Double_t mag = 0;
  Double_t thrust_temp = 0;
  Double_t thrust_max = 0;
  Double_t dot_maj = 0;
  Double_t dot_min = 0;
  Double_t min_temp = 0;
  Double_t maj_temp = 0;
  Double_t thrust_maj_max =0;
  Double_t thrust_min_max = 0;
  TVector3 max_thrust_axis;
  TVector3 p3Norm;
  // Int_t max_nref;

  
  // get the pt vector for each event which passes you jet selections based on pT and eta
  vector <float> pt_v;
  vector <float> eta_v;
  vector <float> phi_v;

  int nref = recoJets.size();
  
  for(int ij = 0; ij<nref; ++ij){
    
    if(fabs(recoJets[ij].eta()) < 2.0){
      pt_v.push_back(recoJets[ij].pt());	
      eta_v.push_back(recoJets[ij].eta());	
      phi_v.push_back(recoJets[ij].phi());
    }
  }
  
  if (debug) cout<<"Total Number of Jets    : "<<nref<<endl;
  if (debug) cout<<"Number of Selected Jets : "<<pt_v.size()<<endl;

  int NJets_Sel = pt_v.size();
  
  if(NJets_Sel < 2) {
    if (debug) cout<<"This event had only 1 Jet"<<endl;
    return;	
  }
  
  float Aj = (double)(pt_v[0] - pt_v[1])/(pt_v[0]+pt_v[1]);
  
  thrust_max = 0;

  vector <double> px;
  vector <double> py;
  vector <double> pz;
      
  for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){

    histos1D["hpT"]->Fill(pt_v[naxis]);
    histos1D["heta"]->Fill(eta_v[naxis]);
    histos1D["hphi"]->Fill(phi_v[naxis]);
    
    float axis_jet_pt = pt_v[naxis];
    float axis_jet_eta = eta_v[naxis];
    float axis_jet_phi = phi_v[naxis];
    px.push_back((double)axis_jet_pt * TMath::Cos(axis_jet_phi));
    py.push_back((double)axis_jet_pt * TMath::Sin(axis_jet_phi));
    pz.push_back((double)axis_jet_pt * TMath::SinH(axis_jet_eta));

    if(py[naxis] == 0) {
      if(debug) cout << "PYZERO" << endl; 
      if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<axis_jet_pt<<"\n \t eta = "<<axis_jet_eta<<"\n \t phi = "<<axis_jet_phi<<endl;

    }
  }

  //PART 1 ====================================================
  //Runs through all the jets in an event, checking them to see if they are the axis that maximizes thrust
  //max axis finding loop
  for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){

    // float axis_jet_pt = pt_v[naxis];
    // float axis_jet_eta = eta_v[naxis];
    // float axis_jet_phi = phi_v[naxis];

    //reset values for this particular event
    thrust_temp = 0;  // maj_temp = 0;   min_temp = 0;
	
    // px[naxis] = pt[naxis]*TMath::Cos(phi[naxis]);
    // py[naxis] = pt[naxis]*TMath::Sin(phi[naxis]);
    // pz[naxis] = pt[naxis]*TMath::SinH(eta[naxis]);
	
    //calculates axis for this particular jet
    TVector3 nT (px[naxis], py[naxis], pz[naxis]);
    if(debug) cout<<"Test Axis Unnormed = {" << nT(0) << ", " << nT(1) << ", " << nT(2)<< "}" << endl;
    nT = Norm(nT);
	
    if(debug) cout<<"Test Axis = {" << nT(0) << ", " << nT(1) << ", " << nT(2)<< "}" << endl;
    //temp_phi = axis_jet_phi;   temp_eta = axis_jet_eta;
	
    //resets for next jet loop
    dot = 0;
    mag = 0;

    //PART 2 ====================================================
    //Loops through all the jets to find the thrust value for the chosen axis 
    //jet loop
    for(Long64_t njet = 0; njet < NJets_Sel; ++njet){
	  
      // if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)){ continue;}
      float jet_pt = pt_v[njet];
      float jet_eta = eta_v[njet];
      float jet_phi = phi_v[njet];
	
      if(debug) cout<< " \n --------- New Jet (Thrust)--------- " << endl; 
      if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<jet_pt<<"\n \t eta = "<<jet_eta<<"\n \t phi = "<<jet_phi<<endl;
	  
      //calculate px, py, pz
      // px[njet] = pt[njet]*TMath::Cos(phi[njet]);
      // py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
      // pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
      
      //define momentum three vector
      TVector3 p3 (px[njet], py[njet], pz[njet]);
      //TVector3 p3Norm = Norm(p3);
	  
      if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	  
      //dots the two vectors for Thrust, Tmin and Tmaj
      dot += TMath::Abs(p3.Dot(nT)); 
      if(debug) cout<<"dot sum = " << dot << endl;
	  
      //sum the total p from the individual p magnitudes
      mag += TMath::Abs(p3.Mag());
      if(debug) cout<<"mag sum = " << mag << endl;
	  
    }//end jet loop
	
    //calculate the thrust
    thrust_temp = ((dot)/mag);
	
    //Compare to see if this axis is a new maximum 
    if(debug) cout<< "\ntemp thrust = " << thrust_temp << endl; 
	
    if(thrust_temp>thrust_max){
      thrust_max = thrust_temp;
      max_thrust_axis = nT;
      //max_eta = temp_eta;
      //max_phi = temp_phi;
      //max_nref = naxis;
	  
    }
    if(debug) cout<< "max thrust = " << thrust_max << endl;
	
  }//end axis loop
      
  if (debug) cout << "FINAL THRUST VALUE: " << thrust_max << endl; 
      
  // FILL BAD THRUST VALUES TO DEBUG =============================
  // if(thrust_max < 0.47) {
	
  //   //h_TBadpT->Fill(pt[naxis], pThat_weight);
  //   h_TBad[cBin]->Fill(thrust_max);
  //   h_etaBad[cBin]->Fill(eta_v[max_nref]);
  //   h_phiBad[cBin]->Fill(phi_v[max_nref]);
  //   h_nrefBad[cBin]->Fill(max_nref);
  //   h_jetCountBad[cBin]->Fill(NJets_Sel);
  //   //h_weightBad->Fill(pThat_weight);
  //   //h_pthatBad->Fill(pThat);
  //   h_fileNum->Fill(ifile);

  //   if (debug) cout << "______________________________" << endl; 
  //   if (debug) cout << "| X0X : THRUST LESS THAN 0.5 |" << endl;
  //   if (debug) cout << "|  Max Thrust: " << thrust_max << endl;
  //   //cout << "|  Max Thrust: " << thrust_max << endl;
  //   if (debug) cout << "______________________________" << endl; 
  // }

  // fill thrust vs Aj plot
  histos2D["hThrust_vs_Aj"]->Fill(Aj, thrust_max);
    
  //PART 3 ====================================================
  //Begin code to select the Thrust Major and Minor axes
      
  //define the plane perpendicular to this axis in order to calculate Tmaj and Tmin
  Dplane* perp = new Dplane(max_thrust_axis);
      
  //reset maximum values for new axis test
  thrust_maj_max = 0;   thrust_min_max = 0; 
      
  //Thrust maj axis loop
  for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){

    // if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)){ continue;}
    if(debug) cout<< " \n --------- New Test Axis (Min/Maj)--------- " << endl; 
	
    //define the jet axis for this iteration
    //calculate px, py, pz
    // px[naxis] = pt[naxis]*TMath::Cos(phi[naxis]);
    // py[naxis] = pt[naxis]*TMath::Sin(phi[naxis]); 
    // pz[naxis] = pt[naxis]*TMath::SinH(eta[naxis]);
	
    //define momentum three vector
    TVector3 p3 (px[naxis], py[naxis], pz[naxis]);
    if(debug) cout<<"Jet Axis UnNormed = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
    p3 = Norm(p3);
    if(debug) cout<<"Jet Axis Normed = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	
    //define maj_axis and min_axis 
    TVector3 maj_axis = perp->Projection((p3));
    if(debug) cout<<"Maj Axis UnNormed = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
    maj_axis = Norm(maj_axis);
    TVector3 min_axis = max_thrust_axis.Cross(maj_axis);
    min_axis = Norm(min_axis); 
	
    if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
    if(debug) cout<<"Maj Axis = {" << maj_axis(0) << ", " << maj_axis(1) << ", " << maj_axis(2)<< "}" << endl;
    if(debug) cout<<"Min Axis = {" << min_axis(0) << ", " << min_axis(1) << ", " << min_axis(2)<< "}\n" << endl;
	
    //reset for new axis test
    dot_maj = 0;   dot_min = 0;   mag = 0;

    //PART 4 ====================================================
    //Test the axis defined by the above loop to determine if this axis is the maximum
    //jet loop
    for(Long64_t njet = 0; njet < NJets_Sel; ++njet){
	  
      //make a ptcut
      // if((pt[njet] < pT_cut)||(TMath::Abs(eta[njet]) > 2)){ continue;}
      float jet_pt = pt_v[njet];
      float jet_eta = eta_v[njet];
      float jet_phi = phi_v[njet];
	  
      if(debug) cout<< " \n --------- New Jet (Maj/Min)--------- " << endl; 
      //if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<pt[njet]<<"\n \t eta = "<<eta[njet]<<"\n \t phi = "<<phi[njet]<<endl;
      if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<jet_pt<<"\n \t eta = "<<jet_eta<<"\n \t phi = "<<jet_phi<<endl;

      //calculate px, py, pz
      // px[njet] = pt[njet]*TMath::Cos(phi[njet]);
      // py[njet] = pt[njet]*TMath::Sin(phi[njet]); 
      // pz[njet] = pt[njet]*TMath::SinH(eta[njet]);
	  
      //define momentum three vector
      TVector3 p3 (px[njet], py[njet], pz[njet]);
      TVector3 p3Norm = Norm(p3);
	  
      if(debug) cout<<"Jet Axis = {" << p3(0) << ", " << p3(1) << ", " << p3(2)<< "}" << endl;
	  
      //dots the two vectors for Tmin and Tmaj
      dot_maj += TMath::Abs(p3.Dot(maj_axis)); 
      dot_min += TMath::Abs(p3.Dot(min_axis));
      if(debug) cout<<"dot maj sum = " << dot_maj << endl;
      if(debug) cout<<"dot min sum = " << dot_min << endl;
	  
      //sum the total p from the individual p magnitudes
      mag += TMath::Abs(p3.Mag());
      if(debug) cout<<"mag sum = " << mag << endl;
	  
    }//end jet loop

    //calculate the thrust major and minor for this axis
    maj_temp = dot_maj/mag;
    min_temp = dot_min/mag;
	
    //test to to see if this particular Tmaj and Tmin are the new maxima
    if(maj_temp>thrust_maj_max){
      thrust_maj_max = maj_temp;  
      thrust_min_max = min_temp; 
      if(debug) cout << "thrust major max = "<< thrust_maj_max<< endl;
      /*
	if(maj_temp > 0.5) {
	h_TmajBadpT->Fill(pt[naxis]);
	h_TmajBad->Fill(maj_temp);
	} 
	if(min_temp > 0.5) {
	h_TminBadpT->Fill(pt[naxis]);
	h_TminBad->Fill(maj_temp);
	} 
      */
    }   
  }//end of major/minor axis loop
      
  //fill all the maximum values before finishing
  // if(jetCount > 1){

  

  histos1D["hThrust"]->Fill(thrust_max);
  histos1D["hThrust_min"]->Fill(thrust_min_max);
  histos1D["hThrust_maj"]->Fill(thrust_maj_max);
  	
  // if(debug) {
  //   if (thrust_max < 0.5)     {  cout << "FLAG_thrust1: " << thrust_max <<  " , " << jetCount << endl; }
  //   if (thrust_maj_max > 0.5) {  cout << "FLAG_maj: " << thrust_maj_max <<  " , " << jetCount << endl; }
  //   if (thrust_min_max > 0.5) {  cout << "FLAG_min: " << thrust_min_max <<  " , " << jetCount << endl; }
  // }	
  //if(jt80)    h_80->Fill(thrust_max,jt80_pre * pThat_weight);
  //if(jt60)    h_60->Fill(thrust_max,jt60_pre * pThat_weight);
  //if(jt40)    h_40->Fill(thrust_max,jt40_pre * pThat_weight);
  //}

  pt_v.clear();
  eta_v.clear();
  phi_v.clear();

  px.clear();
  py.clear();
  pz.clear();      
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThrustAnalyzer);
