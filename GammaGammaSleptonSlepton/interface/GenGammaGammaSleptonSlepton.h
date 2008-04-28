#ifndef UserCode_GenGammaGammaSleptonSlepton
#define UserCode_GenGammaGammaSleptonSlepton

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h" 
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h" 
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "UserCode/GammaGammaSleptonSlepton/interface/AcceptanceTableHelper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>


using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;

class GenGammaGammaSleptonSlepton : public edm::EDAnalyzer {
public:
  explicit GenGammaGammaSleptonSlepton(const edm::ParameterSet&);
  ~GenGammaGammaSleptonSlepton();
  
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
private:
  int nEvt;// used to count the number of events

  
  // to be used for root output tree
  std::string rootfilename;
  TFile *thefile;
  TTree *smalltree;
  int nMCPar;
  int MCPARMAX;// used to set maximum of arrays
  int PROTMAX;
  int MUONMAX;
  int ELEMAX;

  int MCPar_status[1000];
  double MCPar_px[1000];
  double MCPar_py[1000];
  double MCPar_pz[1000];
  double MCPar_e[1000];
  double MCPar_phi[1000];
  double MCPar_eta[1000];
  double MCPar_mass[1000];
  int MCPar_pdgid[1000];
  int MCPar_motherid[1000];

  // Generator level
  int nGenMuonCand;
  double GenMuonCand_px[10];
  double GenMuonCand_py[10];
  double GenMuonCand_pz[10];
  double GenMuonCand_p[10];
  double GenMuonCand_eta[10];
  double GenMuonCand_pt[10];
  double GenMuonCand_phi[10];
  double GenMuonCand_e[10];
  int GenMuonCand_charge[10];
  int GenMuonCand_moth[10];

  int nGenEleCand; 
  double GenEleCand_px[10]; 
  double GenEleCand_py[10]; 
  double GenEleCand_pz[10]; 
  double GenEleCand_p[10]; 
  double GenEleCand_eta[10]; 
  double GenEleCand_pt[10]; 
  double GenEleCand_phi[10]; 
  double GenEleCand_e[10]; 
  int GenEleCand_charge[10]; 
  int GenEleCand_moth[10];

  int nGenProtCand;
  double GenProtCand_px[100];
  double GenProtCand_py[100];
  double GenProtCand_pz[100];
  double GenProtCand_p[100];
  double GenProtCand_eta[100];
  double GenProtCand_pt[100];
  double GenProtCand_phi[100];
  double GenProtCand_e[100];
  double GenProtCand_t[100];
  double GenProtCand_xi[100];
  int GenProtCand_charge[100];
  int GenProtCand_tag[100];
  int nGenSmuonCand;
  double GenSmuonCand_px[10];
  double GenSmuonCand_py[10];
  double GenSmuonCand_pz[10];
  double GenSmuonCand_p[10];
  double GenSmuonCand_eta[10];
  double GenSmuonCand_pt[10];
  double GenSmuonCand_phi[10];
  double GenSmuonCand_mass[10];
  int GenSmuonCand_charge[10];
  int nGenNeutCand;
  double GenNeutCand_px[10];
  double GenNeutCand_py[10];
  double GenNeutCand_pz[10];
  double GenNeutCand_p[10];
  double GenNeutCand_eta[10];
  double GenNeutCand_pt[10];
  double GenNeutCand_phi[10];
  double GenNeutCand_mass[10];

  int nGenPhotCand; 
  double GenPhotCand_px[10]; 
  double GenPhotCand_py[10]; 
  double GenPhotCand_pz[10]; 
  double GenPhotCand_p[10]; 
  double GenPhotCand_eta[10]; 
  double GenPhotCand_pt[10]; 
  double GenPhotCand_phi[10]; 
  double GenPhotCand_e[10]; 

  double GenElEl_dpt;
  double GenElEl_dphi;
  double GenElEl_mass;
  double GenMuMu_mass;
  double GenMuMu_dpt;
  double GenMuMu_dphi;
  double MuMu_mass;
  double GenProPro_mass;

  AcceptanceTableHelper helper420beam1;   
  AcceptanceTableHelper helper420beam2;   
  AcceptanceTableHelper helper220beam1;   
  AcceptanceTableHelper helper220beam2;   
  AcceptanceTableHelper helper420a220beam1;   
  AcceptanceTableHelper helper420a220beam2;   
 
  float acc420b1, acc220b1, acc420and220b1, acc420or220b1; // beam 1 (clockwise)   
  float acc420b2, acc220b2, acc420and220b2, acc420or220b2; // beam 2 (anti-clockwise)   

};
#endif
