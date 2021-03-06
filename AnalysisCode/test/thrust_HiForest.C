// Jennifer Coulter
// July 22nd 2015
// Rutgers University, jennifer.coulter@cern.ch
//
// Test macro for plotting thrust, an event shape variable.
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "TVector3.h"
#include <TROOT.h>

int findBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  if(bin<10)ibin=0; //! 0-5%
  else if(bin>=10  && bin<20 )ibin=1; //! 5-10%
  else if(bin>=20  && bin<60 )ibin=2;  //! 10-30%
  else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
  else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
  else if(bin>=140 && bin<180)ibin=5;  //! 70-90%
  else if(bin>=180 && bin<200)ibin=6;  //! 90-100%
  return ibin;
}

const int nbins_cent = 6;
double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};

using namespace std;

//plane class
class Plane{
public:
  TVector3 v1, v2, proj, u1, u2;
  Double_t scalar1, scalar2, mag1, mag2; 
  Plane(TVector3);

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
Plane::Plane(TVector3 nT){
  
  //Use TVector3 to find an orthogonal vector and a second vector orthogonal to the first and nT
  v1 = nT.Orthogonal();  v2 = nT.Cross(v1);

  //Normalize, checking for 0 length axes
  if ((v1(0) == 0) && (v1(1) == 0) && (v1(2) == 0)){  v1(0) = 0;    v1(1) = 0;    v1(2) = 0; }
  else { mag1 = v1.Mag();   v1(0) = v1(0)/mag1;    v1(1) = v1(1)/mag1;    v1(2) = v1(2)/mag1; } 
  if ((v2(0) == 0) && (v2(1) == 0) && (v2(2) == 0)){  v2(0) = 0;    v2(1) = 0;    v2(2) = 0; } 
  else { mag2 = v2.Mag();   v2(0) = v2(0)/mag2;    v2(1) = v2(1)/mag2;    v2(2) = v2(2)/mag2; }
}//end plane constructor

//creates histograms in terms of Thrust vs. dN/dT
TH1F* DivideByBinWidth(TH1F * hist, const char * name){

  TH1F* h_return = new TH1F(name, "", hist->GetNbinsX(), 0,1);
  hist->Sumw2(); 
  //loops through all the bins
  for (int i=1;i<=hist->GetNbinsX();++i){
    Float_t bin = hist->GetBinWidth(i);
    Float_t val = hist->GetBinContent(i);
    Float_t valErr = h_return->GetBinError(i);
    val = val/bin;
    valErr= valErr/bin;
    h_return->SetBinError(i,valErr);
    h_return->SetBinContent(i, val); 
  }//end bin loop
  return h_return;
}//end rebin function

//Function to normalize a vector
TVector3 Norm(TVector3 v){
  if ( (v(0) == 0) && (v(1) == 0) && (v(2) == 0)) return v; 
  Double_t mag = TMath::Sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); 
  v(0) = v(0)/mag;    v(1) = v(1)/mag;   v(2) = v(2)/mag;
  return v; 
}//end normalize

//plot thrust
//void thrust_HiForest(Int_t startfile, Int_t endfile, Int_t jobNumber){

void thrust_HiForest(Int_t startfile = 0,
		     Int_t endfile = 1,
		     Int_t radius = 3,
		     std::string kFoname="test_output.root",
		     float ptCut = 30.0,
		     float etaCut = 2.0){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  bool debug = false;
  //Float_t pT_cut = 30;
  //Int_t radius = 3;
  
  //define trees and file
  TFile * file; 
  //TFile * weight_file;
  TTree * t;
  TTree * hiEvt;
  //TTree * hlt;
  TTree * skim;
  //TTree * weight;
  TTree * thrust_tree;

  /*
  //set up the tree to be filled
  thrust_tree->Branch("pthatweight",&pthatweight,"pthatweight/D");
  thrust_tree->Branch("hiBin",&hiBin,"hiBin/I");
  thrust_tree->Branch("evt",&evnt,"evt/I");
  thrust_tree->Branch("lumi",&lumi,"lumi/I");
  thrust_tree->Branch("vz",&vz,"vz/F");
  */


  //Tree variables
  Float_t pt[1000];    Int_t jt80;   Int_t jt80_pre;
  Float_t eta[1000];   Int_t jt40;   Int_t jt60_pre;
  Float_t phi[1000];   Int_t jt60;   Int_t jt40_pre;


  // Double_t px[1000];
  // Double_t py[1000];
  // Double_t pz[1000];
  
  Int_t nref;
  Float_t vz;
  Int_t hiBin;
  Int_t halo;
  Int_t noise;
  //Double_t pThat_weight;
  //Float_t pThat; 

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
  //Float_t max_eta = 0;   Float_t temp_eta = 0;  
  //Float_t max_phi = 0;   Float_t temp_phi = 0;
  Int_t max_nref;
  Int_t jetCount = 0; //in order to sum up the number of jets per event
  Int_t eventCount = 0;//check to see how many of the events in each file are actually being used

  //define instream
  string input_file = "jetRAA_PbPb_data_forest.txt"; 
  ifstream count(input_file.c_str(), ifstream::in);
  Int_t fileCount = 0;
  string * filename = new string[10000];

  //count up the total number of files and save their names to an array of filenames for initialization purposes
  
  string line;
  while(getline(count, line)){
    filename[fileCount] = line;
    if (debug) cout << filename[fileCount] << endl; 
    fileCount++;
  }
  count.close();

  
  TFile * save_File = new TFile(kFoname.c_str(),"RECREATE");

  TH1F * h_thrust[nbins_cent], * h_min[nbins_cent], * h_maj[nbins_cent], * h_pT[nbins_cent], * h_pTcut[nbins_cent], * h_nref[nbins_cent],
    * h_jetCount[nbins_cent], * h_eta[nbins_cent], * h_phi[nbins_cent], * h_TBad[nbins_cent],
    * h_TminBad[nbins_cent], * h_TmajBad[nbins_cent], * h_TBadpT[nbins_cent],
    * h_TmajBadpT[nbins_cent],* h_TminBadpT[nbins_cent], * h_etaBad[nbins_cent], * h_phiBad[nbins_cent], * h_nrefBad[nbins_cent], * h_jetCountBad[nbins_cent];

  //TH1F * h_weight = new TH1F("weighting", "", 500, 0, 500);
  //TH1F * h_weightBad = new TH1F("weightingBad", "", 500, 0, 500);
  //TH1F * h_pthatBad = new TH1F("pthatBad", "", 500, 0, 500);
  TH1F * h_fileNum = new TH1F("fileNum", "", 45, 0, 45);

  TH2F * hThrust_vs_Aj[nbins_cent];


  for(int i = 0; i<nbins_cent; ++i){

    hThrust_vs_Aj[i] = new TH2F(Form("hThrust_vs_Aj_cent%d",i),"",50, 0, 1, 50, 0, 1);
    h_thrust[i] = new TH1F(Form("thrust_unscaled_cent%d",i), "", 50,0,1);
    h_min[i] = new TH1F(Form("thrust_min_cent%d",i), "", 50,0,1);
    h_maj[i] = new TH1F(Form("thrust_maj_cent%d",i), "", 50,0,1);
    h_pT[i] = new TH1F(Form("pT_cent%d",i), "", 100, 0, 120);
    h_pTcut[i] = new TH1F(Form("pTcut_cent%d",i), "", 100, 0, 120);
    //h_40[i] = new TH1F(Form("thrust_40_cent%d",i), "", 50,0,1);
    //h_60[i] = new TH1F(Form("thrust_60_cent%d",i), "", 50,0,1);
    //h_80[i] = new TH1F(Form("thrust_80_cent%d",i), "", 50,0,1);
    h_nref[i] = new TH1F(Form("nref_cent%d",i), "", 12, 0, 12);
    h_jetCount[i] = new TH1F(Form("jetCount_cent%d",i), "", 12, 0, 12);
    h_eta[i] = new TH1F(Form("eta_cent%d",i), "", 60, -2, 2);
    h_phi[i] = new TH1F(Form("phi_cent%d",i), "", 60, -3.15, 3.15);
  
    h_TBad[i] = new TH1F(Form("thrust_bad_cent%d",i), "", 50,0,1);
    h_TminBad[i] = new TH1F(Form("thrust_min_bad_cent%d",i), "", 50,0,1);
    h_TmajBad[i] = new TH1F(Form("thrust_maj_bad_cent%d",i), "", 50,0,1);
    h_TBadpT[i] = new TH1F(Form("TBadpT_cent%d",i), "", 100, 0, 120);
    h_TmajBadpT[i] = new TH1F(Form("TmajBadpT_cent%d",i), "", 100, 0, 120);
    h_TminBadpT[i] = new TH1F(Form("TminBadpT_cent%d",i), "", 100, 0, 120);
    h_etaBad[i] = new TH1F(Form("etaBad_cent%d",i), "", 60, -2, 2);
    h_phiBad[i] = new TH1F(Form("phiBad_cent%d",i), "", 60, -3.15, 3.15);
    h_nrefBad[i] = new TH1F(Form("nrefBad_cent%d",i), "", 12, 0, 12);
    h_jetCountBad[i] = new TH1F(Form("jetCountBad_cent%d",i), "", 12, 0, 12);
  }
  
  // For every file in file list, process trees
  for(int ifile = startfile; ifile < endfile; ifile++){

    string s = filename[ifile];
    //string w = Form("weights_pbpb_%d.root", ifile+1); 
    file = TFile::Open(s.c_str());
    //weight_file = TFile::Open(w.c_str());

    if (debug) cout << "\n **** =========================== New File ================================= **** \n ";
    cout << "File Name: " << filename[ifile] << endl; 
    cout << "File Number: " << endfile-ifile << "/" << endfile-startfile << endl;
    //cout << "Weight File: " << w << endl;

    //define trees and file
    t = (TTree*)file->Get(Form("akPu%dPFJetAnalyzer/t", radius));
    hiEvt = (TTree*)file->Get("hiEvtAnalyzer/HiTree");
    //hlt = (TTree*)file->Get("hltanalysis/HltTree");
    skim = (TTree*)file->Get("skimanalysis/HltTree");
    //weight = (TTree*)weight_file->Get("weights");
    
    //Set branches of the tree 
    t->SetBranchAddress("jtpt", &pt);
    t->SetBranchAddress("jteta", &eta);
    t->SetBranchAddress("jtphi", &phi);
    t->SetBranchAddress("nref", &nref);
    //t->SetBranchAddress("pthat", &pThat);

    hiEvt->SetBranchAddress("vz", &vz);
    hiEvt->SetBranchAddress("hiBin", &hiBin);

    skim->SetBranchAddress("pHBHENoiseFilter", &noise);
    skim->SetBranchAddress("pcollisionEventSelection",&halo);
  
    //hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jt80);
    //hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jt60);
    //hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jt40);
    //hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jt80_pre);
    //hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jt60_pre);
    //hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jt40_pre);

    //weight->SetBranchAddress("pthatweight", &pThat_weight);

    t->AddFriend(hiEvt);
    t->AddFriend(skim);
    //t->AddFriend(hlt);
    //t->AddFriend(weight);
    
    Long64_t nentries = t->GetEntries();
    //nentries = 10000;
    
    cout << "Events in File: " << nentries << endl;
    eventCount = 0;

    save_File->cd();
    
    //event loop
    for(Long64_t nentry = 0; nentry<nentries; ++nentry){
      //for(Long64_t nentry = 6662; nentry<6663; ++nentry){
      
      if(nentry%10000 == 0) cout << nentry << endl;
      
      t->GetEntry(nentry);
      skim->GetEntry(nentry);
      hiEvt->GetEntry(nentry);
      //weight->GetEntry(nentry);

      int cBin = findBin(hiBin);//tells us the centrality of the event. 
      if(cBin==-1 || cBin==nbins_cent) continue;
      
      jetCount = 0;
      bool select = false;

      //make selection cuts
      if(TMath::Abs(vz) > 15 || halo == 0 || noise == 0) continue;

      // apply the pT dijet clean sample selection
      if(pt[0]< 120 || pt[1] <50) continue; 
      
      //fill pThat spectra plot
      //h_weight->Fill(pThat, pThat_weight); 
      
      //if((TMath::Abs(vz) > 15)) {continue;}

      // get the pt vector for each event which passes you jet selections based on pT and eta
      vector <float> pt_v;
      vector <float> eta_v;
      vector <float> phi_v;

      for(int ij = 0; ij<nref; ++ij){

	if(pt[ij] > ptCut && fabs(eta[ij]) < etaCut){
	  pt_v.push_back(pt[ij]);	
	  eta_v.push_back(eta[ij]);	
	  phi_v.push_back(phi[ij]);
	}
      }

      if (debug) cout<<"Total Number of Jets    : "<<nref<<endl;
      if (debug) cout<<"Number of Selected Jets : "<<pt_v.size()<<endl;

      int NJets_Sel = pt_v.size();

      if(NJets_Sel < 2) {
	if (debug) cout<<"This event had only 1 Jet"<<endl;
	continue;	
      }

      double Aj = (double)(pt_v[0] - pt_v[1])/(pt_v[0]+pt_v[1]);
      
      // //intial pT cut
      // for(int k = 0; k < nref; k++){
      // 	if(pt[k] >  30){
      // 	  if((TMath::Abs(eta[k])<2)) select = true; 
      // 	}
      // }

      // if(!select) continue;
      if(debug) cout<< " \n ******* New Event ******** " << endl;
      if(debug) cout<< " ******* " << nentry << " ******** " << endl;
      
      //reset maximum values
      eventCount++;
      thrust_max = 0;

      vector <double> px;
      vector <double> py;
      vector <double> pz;
      
      for(Long64_t naxis = 0; naxis < NJets_Sel; ++naxis){
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

	float axis_jet_pt = pt_v[naxis];
	float axis_jet_eta = eta_v[naxis];
	float axis_jet_phi = phi_v[naxis];
	
	//Cut checks
	//if(debug) cout<< "Jet Variables: "<< "\n \t pT = "<<axis_jet_pt<<"\n \t eta = "<<axis_jet_eta<<"\n \t phi = "<<axis_jet_phi<<endl;
	// h_pT->Fill(pt[naxis]);

	// if((pt[naxis] < pT_cut)||(TMath::Abs(eta[naxis]) > 2)) {
	//   continue;
	// }
	h_pTcut[cBin]->Fill(axis_jet_pt);
	
	// jetCount++; 
	
	if(debug) cout<< " \n --------- New Test Axis (Thrust)--------- " << endl; 
	
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
	dot = 0;   mag = 0;

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
	  max_nref = naxis;
	  
	}
	if(debug) cout<< "max thrust = " << thrust_max << endl;
	
      }//end axis loop
      
      if (debug) cout << "FINAL THRUST VALUE: " << thrust_max << endl; 
      
      // FILL BAD THRUST VALUES TO DEBUG =============================
      if(thrust_max < 0.47) {
	
	//h_TBadpT->Fill(pt[naxis], pThat_weight);
	h_TBad[cBin]->Fill(thrust_max);
	h_etaBad[cBin]->Fill(eta_v[max_nref]);
	h_phiBad[cBin]->Fill(phi_v[max_nref]);
	h_nrefBad[cBin]->Fill(max_nref);
	h_jetCountBad[cBin]->Fill(NJets_Sel);
	//h_weightBad->Fill(pThat_weight);
	//h_pthatBad->Fill(pThat);
	h_fileNum->Fill(ifile);

	if (debug) cout << "______________________________" << endl; 
	if (debug) cout << "| X0X : THRUST LESS THAN 0.5 |" << endl;
	if (debug) cout << "|  Max Thrust: " << thrust_max << endl;
	//cout << "|  Max Thrust: " << thrust_max << endl;
	if (debug) cout << "______________________________" << endl; 
      }

      // fill thrust vs Aj plot
      hThrust_vs_Aj[cBin]->Fill(Aj, thrust_max);

    
      //PART 3 ====================================================
      //Begin code to select the Thrust Major and Minor axes
      
      //define the plane perpendicular to this axis in order to calculate Tmaj and Tmin
      Plane* perp = new Plane(max_thrust_axis);
      
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
      timer.Stop();
      
      //fill all the maximum values before finishing
      // if(jetCount > 1){
	
      h_thrust[cBin]->Fill(thrust_max);
      h_eta[cBin]->Fill(eta_v[max_nref]);
      h_phi[cBin]->Fill(phi_v[max_nref]);
      h_min[cBin]->Fill(thrust_min_max);
      h_maj[cBin]->Fill(thrust_maj_max);
      h_nref[cBin]->Fill(nref);
      h_jetCount[cBin]->Fill(NJets_Sel);
	
      if(debug) {
	if (thrust_max < 0.5)     {  cout << "FLAG_thrust1: " << thrust_max <<  " , " << jetCount << endl; }
	if (thrust_maj_max > 0.5) {  cout << "FLAG_maj: " << thrust_maj_max <<  " , " << jetCount << endl; }
	if (thrust_min_max > 0.5) {  cout << "FLAG_min: " << thrust_min_max <<  " , " << jetCount << endl; }
      }	
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
      
    }//end of event loop


    
    gROOT->GetListOfFiles()->Remove(file);
    //gROOT->GetListOfFiles()->Remove(weight);
    
    cout << "Events Selected: " << eventCount << endl;
    cout << "File Finished" << endl; 
    
  }//end file loop

  //scale the histograms
  Double_t integral;
  
  // integral = h_thrust->Integral();
  // h_thrust->Scale(1/integral);
  //h_thrust->Scale(1./h_thrust->Integral());
  //h_maj->Scale(1./h_maj->Integral());
  //h_min->Scale(1./h_min->Integral());

  // integral = h_maj->Integral(); 
  // h_maj->Scale(1/integral);

  // integral = h_min->Integral();
  // h_min->Scale(1/integral);

  //integral = h_80->Integral();
  //h_80->Scale(integral);

  //integral = h_60->Integral();
  //h_60->Scale(integral);

  // integral = h_40->Integral();
  //h_40->Scale(integral); 
 
  //Create the plot for Thrust vs. dN/dT
  //define histograms
  // TH1F * h_T = DivideByBinWidth(h_thrust[], "thrust_scaled");
  // TH1F * h_Tmaj = DivideByBinWidth(h_maj, "thrust_maj_scaled");
  // TH1F * h_Tmin = DivideByBinWidth(h_min, "thrust_min_scaled");
  
  //Float_t divEntries = 1./(h_thrust->GetBinWidth(1));
  
  //h_T->Scale(divEntries);
  //h_Tmaj->Scale(divEntries);
  //h_Tmin->Scale(divEntries);
  
  //h_40 = DivideByBinWidth(h_40, "thrust_40_new");   h_40->Scale(divEntries);
  //h_60 = DivideByBinWidth(h_60, "thrust_60_new");   h_60->Scale(divEntries);
  //h_80 = DivideByBinWidth(h_80, "thrust_80_new");   h_80->Scale(divEntries);

  // h_T->Print("base");
  // h_Tmaj->Print("base");
  // h_Tmin->Print("base");
  // h_pT->Print("base");
  // h_pTcut->Print("base");
  // //h_40->Print("base");
  // //h_60->Print("base");
  // //h_80->Print("base");
  // h_nref->Print("base");
  // h_jetCount->Print("base");
  // h_eta->Print("base");
  // h_phi->Print("base");
  // h_weight->Print("base"); 
  
  // h_T->Write();
  // h_Tmaj->Write();
  // h_Tmin->Write();
  // h_pT->Write();
  // h_pTcut->Write();
  // //h_40->Write();
  // //h_60->Write();
  // //h_80->Write();
  // h_nref->Write();
  // h_jetCount->Write();
  // h_eta->Write();
  // h_phi->Write();
  // h_weight->Write(); 
  
  save_File->Write();
  save_File->Close();

  cout<<endl<<endl<<endl<<endl<<"GOING TO WRITE BAD OUTPUT FILE"<<endl<<endl<<endl<<endl;
  
  //if(debug) {
    
  // TFile * bad_File = new TFile(kFoname2.c_str(),"RECREATE");  
  // bad_File->cd();
  
  // h_TBad->Print("base");
  // h_TmajBad->Print("base");
  // h_TminBad->Print("base");
  // h_TBadpT->Print("base");
  // h_TmajBadpT->Print("base");
  // h_TminBadpT->Print("base");
  // h_etaBad->Print("base");
  // h_phiBad->Print("base");
  // h_nrefBad->Print("base");
  // h_jetCountBad->Print("base");
  // h_pthatBad->Print("base");
  // h_weightBad->Print("base");
  // h_fileNum->Print("base"); 
    
  // h_TBad->Write();
  // h_TmajBad->Write();
  // h_TminBad->Write();
  // h_TBadpT->Write();
  // h_TmajBadpT->Write();
  // h_TminBadpT->Write();
  // h_etaBad->Write();
  // h_phiBad->Write();
  // h_nrefBad->Write();
  // h_jetCountBad->Write();
  // h_pthatBad->Write();
  // h_weightBad->Write();
  // h_fileNum->Write();
    
  // bad_File->Write();
  // bad_File->Close();
  // // }
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//end of plot thrust
