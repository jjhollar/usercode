
// -*- C++ -*-
//
// Package:    GammaGammaSleptonSlepton
// Class:      GammaGammaSleptonSlepton
// 

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <vector> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/MuonReco/interface/Muon.h" 
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h" 
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "UserCode/GammaGammaSleptonSlepton/interface/GammaGammaSleptonSlepton.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

// Vertexing
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;


GammaGammaSleptonSlepton::GammaGammaSleptonSlepton(const edm::ParameterSet& iConfig)
{
  nEvt=0;
  MCPARMAX=1000;
  ELEMAX=100;
  MUONMAX=100;
  JETMAX=300;
  TRACKMAX=300;
  PHOTMAX=100;
  PROTMAX=100;

  rootfilename       = iConfig.getUntrackedParameter<std::string>("outfilename","gammagamma.anal.root");

  thefile = new TFile(rootfilename.c_str(),"recreate"); 
  thefile->cd(); 
  smalltree= new TTree("ntp1","ntp1");
  smalltree->Branch("nMCPar",&nMCPar,"nMCPar/I");
  //  smalltree->Branch("MCPar_status",&MCPar_status,"MCPar_status[nMCPar]/I");
  //  smalltree->Branch("MCPar_px",MCPar_px,"MCPar_px[nMCPar]/D");
  //  smalltree->Branch("MCPar_py",MCPar_py,"MCPar_py[nMCPar]/D");
  //  smalltree->Branch("MCPar_pz",MCPar_pz,"MCPar_pz[nMCPar]/D");
  //  smalltree->Branch("MCPar_e",MCPar_e,"MCPar_e[nMCPar]/D");
  //  smalltree->Branch("MCPar_phi",MCPar_phi,"MCPar_phi[nMCPar]/D");
  //  smalltree->Branch("MCPar_eta",MCPar_eta,"MCPar_eta[nMCPar]/D");
  //  smalltree->Branch("MCPar_pdgid",MCPar_pdgid,"MCPar_pdgid[nMCPar]/I");
  //  smalltree->Branch("MCPar_mass",MCPar_mass,"MCPar_mass[nMCPar]/D");

//   smalltree->Branch("nEleCand",&nEleCand,"nEleCand/I");
//   smalltree->Branch("EleCand_px",EleCand_px,"EleCand_px[nEleCand]/D");
//   smalltree->Branch("EleCand_py",EleCand_py,"EleCand_py[nEleCand]/D");
//   smalltree->Branch("EleCand_pz",EleCand_pz,"EleCand_pz[nEleCand]/D");
//   smalltree->Branch("EleCand_e",EleCand_e,"EleCand_e[nEleCand]/D");

//   smalltree->Branch("nJetCand",&nJetCand,"nJetCand/I");
//   smalltree->Branch("JetCand_px",JetCand_px,"JetCand_px[nJetCand]/D");
//   smalltree->Branch("JetCand_py",JetCand_py,"JetCand_py[nJetCand]/D");
//   smalltree->Branch("JetCand_pz",JetCand_pz,"JetCand_pz[nJetCand]/D");
//   smalltree->Branch("JetCand_e",JetCand_e,"JetCand_e[nJetCand]/D");
//   smalltree->Branch("JetCand_eta",JetCand_eta,"JetCand_eta[nJetCand]/D");
//   smalltree->Branch("JetCand_phi",JetCand_phi,"JetCand_phi[nJetCand]/D");

  smalltree->Branch("nJetCand",&nJetCand,"nJetCand/I");
  smalltree->Branch("nJetCandE20",&nJetCandE20,"nJetCandE20/I");
  smalltree->Branch("JetCand_px",JetCand_px,"JetCand_px[nJetCand]/D");
  smalltree->Branch("JetCand_py",JetCand_py,"JetCand_py[nJetCand]/D");
  smalltree->Branch("JetCand_pz",JetCand_pz,"JetCand_pz[nJetCand]/D");
  smalltree->Branch("JetCand_e",JetCand_e,"JetCand_e[nJetCand]/D");
  smalltree->Branch("JetCand_eta",JetCand_eta,"JetCand_eta[nJetCand]/D");
  smalltree->Branch("JetCand_phi",JetCand_phi,"JetCand_phi[nJetCand]/D");

  smalltree->Branch("nMuonCand",&nMuonCand,"nMuonCand/I");
  smalltree->Branch("MuonCand_px",MuonCand_px,"MuonCand_px[nMuonCand]/D");
  smalltree->Branch("MuonCand_py",MuonCand_py,"MuonCand_py[nMuonCand]/D");
  smalltree->Branch("MuonCand_pz",MuonCand_pz,"MuonCand_pz[nMuonCand]/D");
  smalltree->Branch("MuonCand_p",MuonCand_p,"MuonCand_p[nMuonCand]/D");
  smalltree->Branch("MuonCand_pt",MuonCand_pt,"MuonCand_pt[nMuonCand]/D");
  smalltree->Branch("MuonCand_eta",MuonCand_eta,"MuonCand_eta[nMuonCand]/D");
  smalltree->Branch("MuonCand_phi",MuonCand_phi,"MuonCand_phi[nMuonCand]/D");
  smalltree->Branch("MuonCand_hcalisor3",MuonCand_hcalisor3,"MuonCand_hcalisor3[nMuonCand]/D");
  smalltree->Branch("MuonCand_ecalisor3",MuonCand_ecalisor3,"MuonCand_ecalisor3[nMuonCand]/D"); 
  smalltree->Branch("MuonCand_trkisor3",MuonCand_trkisor3,"MuonCand_trkisor3[nMuonCand]/D"); 
  smalltree->Branch("MuonCand_hcalisor5",MuonCand_hcalisor5,"MuonCand_hcalisor5[nMuonCand]/D"); 
  smalltree->Branch("MuonCand_ecalisor5",MuonCand_ecalisor5,"MuonCand_ecalisor5[nMuonCand]/D"); 
  smalltree->Branch("MuonCand_trkisor5",MuonCand_trkisor5,"MuonCand_trkisor5[nMuonCand]/D");
  smalltree->Branch("MuonCand_vtxx",MuonCand_vtxx,"MuonCand_vtxx[nMuonCand]/D");
  smalltree->Branch("MuonCand_vtxy",MuonCand_vtxy,"MuonCand_vtxy[nMuonCand]/D"); 
  smalltree->Branch("MuonCand_vtxz",MuonCand_vtxz,"MuonCand_vtxz[nMuonCand]/D"); 

  smalltree->Branch("nTrackCand",&nTrackCand,"nTrackCand/I"); 
  smalltree->Branch("TrackCand_px",TrackCand_px,"TrackCand_px[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_py",TrackCand_py,"TrackCand_py[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_pz",TrackCand_pz,"TrackCand_pz[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_p",TrackCand_p,"TrackCand_p[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_pt",TrackCand_pt,"TrackCand_pt[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_eta",TrackCand_eta,"TrackCand_eta[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_phi",TrackCand_phi,"TrackCand_phi[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_vtxx",TrackCand_vtxx,"TrackCand_vtxx[nTrackCand]/D"); 
  smalltree->Branch("TrackCand_vtxy",TrackCand_vtxy,"TrackCand_vtxy[nTrackCand]/D");  
  smalltree->Branch("TrackCand_vtxz",TrackCand_vtxz,"TrackCand_vtxz[nTrackCand]/D");  

  smalltree->Branch("nGenMuonCand",&nGenMuonCand,"nGenMuonCand/I");
  smalltree->Branch("GenMuonCand_px",GenMuonCand_px,"GenMuonCand_px[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_py",GenMuonCand_py,"GenMuonCand_py[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_pz",GenMuonCand_pz,"GenMuonCand_pz[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_p",GenMuonCand_p,"GenMuonCand_p[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_pt",GenMuonCand_pt,"GenMuonCand_pt[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_eta",GenMuonCand_eta,"GenMuonCand_eta[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_phi",GenMuonCand_phi,"GenMuonCand_phi[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_e",GenMuonCand_e,"GenMuonCand_e[nGenMuonCand]/D");

  smalltree->Branch("nGenProtCand",&nGenProtCand,"nGenProtCand/I");
  smalltree->Branch("GenProtCand_px",GenProtCand_px,"GenProtCand_px[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_py",GenProtCand_py,"GenProtCand_py[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_pz",GenProtCand_pz,"GenProtCand_pz[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_p",GenProtCand_p,"GenProtCand_p[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_pt",GenProtCand_pt,"GenProtCand_pt[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_eta",GenProtCand_eta,"GenProtCand_eta[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_phi",GenProtCand_phi,"GenProtCand_phi[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_e",GenProtCand_e,"GenProtCand_e[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_t",GenProtCand_t,"GenProtCand_t[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_xi",GenProtCand_xi,"GenProtCand_xi[nGenProtCand]/D");
  smalltree->Branch("GenProtCand_tag",GenProtCand_tag,"GenProtCand_tag[nGenProtCand]/I");

  smalltree->Branch("nGenSmuonCand",&nGenSmuonCand,"nGenSmuonCand/I");
  smalltree->Branch("GenSmuonCand_px",GenSmuonCand_px,"GenSmuonCand_px[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_py",GenSmuonCand_py,"GenSmuonCand_py[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_pz",GenSmuonCand_pz,"GenSmuonCand_pz[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_p",GenSmuonCand_p,"GenSmuonCand_p[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_pt",GenSmuonCand_pt,"GenSmuonCand_pt[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_eta",GenSmuonCand_eta,"GenSmuonCand_eta[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_phi",GenSmuonCand_phi,"GenSmuonCand_phi[nGenSmuonCand]/D");
  smalltree->Branch("GenSmuonCand_mass",GenSmuonCand_mass,"GenSmuonCand_mass[nGenSmuonCand]/D");

  smalltree->Branch("nGenNeutCand",&nGenNeutCand,"nGenNeutCand/I");
  smalltree->Branch("GenNeutCand_px",GenNeutCand_px,"GenNeutCand_px[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_py",GenNeutCand_py,"GenNeutCand_py[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_pz",GenNeutCand_pz,"GenNeutCand_pz[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_p",GenNeutCand_p,"GenNeutCand_p[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_pt",GenNeutCand_pt,"GenNeutCand_pt[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_eta",GenNeutCand_eta,"GenNeutCand_eta[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_phi",GenNeutCand_phi,"GenNeutCand_phi[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_phi",GenNeutCand_phi,"GenNeutCand_phi[nGenNeutCand]/D");
  smalltree->Branch("GenNeutCand_mass",GenNeutCand_mass,"GenNeutCand_mass[nGenNeutCand]/D");
  smalltree->Branch("eventWeight",&eventWeight,"eventWeight/D");

  //  smalltree->Branch("nGenPhotCand",&nGenPhotCand,"nGenPhotCand/I"); 
  //  smalltree->Branch("GenPhotCand_px",GenPhotCand_px,"GenPhotCand_px[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_py",GenPhotCand_py,"GenPhotCand_py[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_pz",GenPhotCand_pz,"GenPhotCand_pz[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_p",GenPhotCand_p,"GenPhotCand_p[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_pt",GenPhotCand_pt,"GenPhotCand_pt[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_eta",GenPhotCand_eta,"GenPhotCand_eta[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_phi",GenPhotCand_phi,"GenPhotCand_phi[nGenPhotCand]/D"); 
  //  smalltree->Branch("GenPhotCand_e",GenPhotCand_e,"GenPhotCand_e[nGenPhotCand]/D"); 

  smalltree->Branch("GenMuMu_mass",&GenMuMu_mass,"GenMuMu_mass/D");
  smalltree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  smalltree->Branch("MuMu_vtxx",&MuMu_vtxx,"MuMu_vtxx/D");
  smalltree->Branch("MuMu_vtxy",&MuMu_vtxy,"MuMu_vtxy/D"); 
  smalltree->Branch("MuMu_vtxz",&MuMu_vtxz,"MuMu_vtxz/D"); 
  smalltree->Branch("MuMu_extratracksz10mm",&MuMu_extratracksz10mm,"MuMu_extratracksz10mm/I");
  smalltree->Branch("MuMu_extratracksz5mm",&MuMu_extratracksz5mm,"MuMu_extratracksz5mm/I"); 
  smalltree->Branch("MuMu_extratracksz2mm",&MuMu_extratracksz2mm,"MuMu_extratracksz2mm/I"); 
  smalltree->Branch("MuMu_extratracksz1mm",&MuMu_extratracksz1mm,"MuMu_extratracksz1mm/I"); 

  //  smalltree->Branch("GenProPro_mass",&GenProPro_mass,"GenProPro_mass/D");

}


GammaGammaSleptonSlepton::~GammaGammaSleptonSlepton() 
{  
  thefile->Write();
  thefile->Close();
}

void GammaGammaSleptonSlepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // this analyzer produces a small root file with basic candidates and some MC information
  // some additional print statements
  nEvt++;
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;

  // Trigger results
  //  std::vector<std::string>  hlNames_;
  //  edm::TriggerNames triggerNames_;

  //  edm::Handle<edm::TriggerResults> HLTR ;
  //  iEvent.getByLabel("TriggerResults::HLT",HLTR) ;

  //  triggerNames_.init(*HLTR);
  //  hlNames_=triggerNames_.triggerNames();
  //  const unsigned int n(hlNames_.size());

  //  for (unsigned int i=0; i!=n; ++i) {
  //    cout << hlNames_[i] << endl;
  //    if (HLTR->accept(i))
  //      cout << "\tAccepted" << endl;
  //  }

  
  
  // step 1: fill some basic MC information into the root tree
  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel("source", evt);
 
  nMCPar=0;
  nMuonCand = 0;
  nGenMuonCand = 0;
  nGenProtCand = 0;
  nGenSmuonCand = 0;
  nGenNeutCand = 0;
  nGenPhotCand = 0;
  MuMu_mass = 0.0;
  MuMu_extratracksz1mm = 0; MuMu_extratracksz2mm = 0; MuMu_extratracksz5mm = 0; MuMu_extratracksz10mm = 0; 
  double highesteproton = 0.0;
  double highestptmuon = 0.0;

  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(evt->GetEvent()));

  //  std::cout << "Got GenEvent" << std::endl;
  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end() && nMCPar<MCPARMAX ; ++p ) {
    if ( abs((*p)->pdg_id()) !=0 && (abs((*p)->pdg_id())<30 || abs((*p)->pdg_id() == 2212) || abs((*p)->pdg_id()) > 1000000) && nMCPar<MCPARMAX){
       MCPar_status[nMCPar]=(*p)->status();
       MCPar_px[nMCPar]=(*p)->momentum().x();
       MCPar_py[nMCPar]=(*p)->momentum().y();
       MCPar_pz[nMCPar]=(*p)->momentum().z();
       MCPar_phi[nMCPar]=atan2(MCPar_py[nMCPar],MCPar_px[nMCPar]);
       MCPar_eta[nMCPar] = (*p)->momentum().eta();
       MCPar_pdgid[nMCPar]=(*p)->pdg_id();
       MCPar_mass[nMCPar]=(*p)->momentum().m();
       MCPar_e[nMCPar] = sqrt(MCPar_mass[nMCPar]*MCPar_mass[nMCPar] + (MCPar_px[nMCPar]*MCPar_px[nMCPar] + MCPar_py[nMCPar]*MCPar_py[nMCPar] + MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 

      if((*p)->status() < 3)
	{
	  if((((*p)->pdg_id() == 13) || ((*p)->pdg_id() == -13)) && (nGenMuonCand < MUONMAX))
	    {
	      GenMuonCand_px[nGenMuonCand] = MCPar_px[nMCPar];
	      GenMuonCand_py[nGenMuonCand] = MCPar_py[nMCPar];
	      GenMuonCand_pz[nGenMuonCand] = MCPar_pz[nMCPar];
	      GenMuonCand_phi[nGenMuonCand] = MCPar_phi[nMCPar];
	      GenMuonCand_charge[nGenMuonCand] = 
		(int)(MCPar_pdgid[nMCPar]/13);
	      GenMuonCand_eta[nGenMuonCand] = MCPar_eta[nMCPar];
	      GenMuonCand_e[nGenMuonCand] = MCPar_e[nMCPar];
	      GenMuonCand_p[nGenMuonCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
						 (MCPar_py[nMCPar]*MCPar_py[nMCPar])+
						 (MCPar_pz[nMCPar]*MCPar_pz[nMCPar]));

	      nGenMuonCand++;
	    }
	  //	  cout << MCPar_pz[nMCPar] << ", " << MCPar_pdgid[nMCPar] << endl;
	  if(((*p)->pdg_id() == 2212) && (nGenProtCand < PROTMAX))
	    {
	      GenProtCand_px[nGenProtCand] = MCPar_px[nMCPar];
	      GenProtCand_py[nGenProtCand] = MCPar_py[nMCPar];
	      GenProtCand_pz[nGenProtCand] = MCPar_pz[nMCPar];
	      //	       GenProtCand_eta[nGenProtCand] = MCPar_eta[nMCPar];
	      GenProtCand_phi[nGenProtCand] = MCPar_phi[nMCPar];
	      GenProtCand_charge[nGenProtCand] = 
		(int)(MCPar_pdgid[nMCPar]/2212);
	      GenProtCand_eta[nGenProtCand] = MCPar_eta[nMCPar];
	      GenProtCand_e[nGenProtCand] = MCPar_e[nMCPar];
	      GenProtCand_p[nGenProtCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
						 (MCPar_py[nMCPar]*MCPar_py[nMCPar])+
						 (MCPar_pz[nMCPar]*MCPar_pz[nMCPar]));
	      if(GenProtCand_e[nGenProtCand] > highesteproton)
		highesteproton = GenProtCand_e[nGenProtCand];
	      GenProtCand_tag[nGenProtCand] = 0;
	      
	      // FP420-Totem stuff
	      double pt = sqrt(GenProtCand_px[nGenProtCand]*GenProtCand_px[nGenProtCand] +  
			       GenProtCand_py[nGenProtCand]*GenProtCand_py[nGenProtCand]); 
	      double phi = GenProtCand_phi[nGenProtCand];  
	      double mp = 0.938272029;  
	      // ... compute kinimatical variable   
   
	      float xi  = 1.0;    // fractional momentum loss   
	      if (GenProtCand_pz[nGenProtCand]>0)   
		xi -= GenProtCand_pz[nGenProtCand]/7000.0;   
	      else   
		xi += GenProtCand_pz[nGenProtCand]/7000.0;   
   
	      double t   = (-pt*pt - mp*mp*xi*xi) / (1-xi); // "t"   
 
	      if (xi<0.0) xi=-10.0;  
	      if (xi>1.0) xi=10.0;  
 
	      //	      cout << "JJH: pt = " << pt << ", pz = " << GenProtCand_pz[nGenProtCand] << ", t = " << t << ", xi = " << xi << endl; 
 
	      GenProtCand_t[nGenProtCand] = t; 
	      GenProtCand_xi[nGenProtCand] = xi; 

	      if((GenProtCand_pz[nGenProtCand] > 2500.0) || (GenProtCand_pz[nGenProtCand] < -2500.0))
		{ 
		} 

	      nGenProtCand++;
	    }
	  if(((*p)->pdg_id() == 2000013) || ((*p)->pdg_id() == -2000013))
	    {
	      GenSmuonCand_px[nGenSmuonCand] = MCPar_px[nMCPar];
	      GenSmuonCand_py[nGenSmuonCand] = MCPar_py[nMCPar];
	      GenSmuonCand_pz[nGenSmuonCand] = MCPar_pz[nMCPar];
	      GenSmuonCand_eta[nGenSmuonCand] = MCPar_eta[nMCPar];
	      GenSmuonCand_phi[nGenSmuonCand] = MCPar_phi[nMCPar];
	      GenSmuonCand_charge[nGenSmuonCand] = 
		(int)(MCPar_pdgid[nMCPar]/2000013);
	      GenSmuonCand_eta[nGenSmuonCand] = MCPar_eta[nMCPar];
	      GenSmuonCand_mass[nGenSmuonCand] = MCPar_mass[nMCPar];
	      nGenSmuonCand++;
	    }
	  if(((*p)->pdg_id() == 1000022))
	    {
	      GenNeutCand_px[nGenNeutCand] = MCPar_px[nMCPar];
	      GenNeutCand_py[nGenNeutCand] = MCPar_py[nMCPar];
	      GenNeutCand_pz[nGenNeutCand] = MCPar_pz[nMCPar];
	      GenNeutCand_eta[nGenNeutCand] = MCPar_eta[nMCPar];
	      GenNeutCand_phi[nGenNeutCand] = MCPar_phi[nMCPar];
	      GenNeutCand_eta[nGenNeutCand] = MCPar_eta[nMCPar];
	      GenNeutCand_mass[nGenNeutCand] = MCPar_mass[nMCPar];
	      nGenNeutCand++;
	    }
	  //	  if(((*p)->pdg_id() == 22) && (nGenPhotCand < PHOTMAX)) 
	  //	    { 
	  //	      GenPhotCand_px[nGenPhotCand] = MCPar_px[nMCPar]; 
	  //	      GenPhotCand_py[nGenPhotCand] = MCPar_py[nMCPar]; 
	  //	      GenPhotCand_pz[nGenPhotCand] = MCPar_pz[nMCPar]; 
	  //	      GenPhotCand_eta[nGenPhotCand] = MCPar_eta[nMCPar]; 
	  //	      GenPhotCand_phi[nGenPhotCand] = MCPar_phi[nMCPar]; 
	  //	      GenPhotCand_eta[nGenPhotCand] = MCPar_eta[nMCPar]; 
	  //	      GenPhotCand_e[nGenPhotCand] = MCPar_e[nMCPar]; 
	  //	      GenPhotCand_p[nGenPhotCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+  
	  //						 (MCPar_py[nMCPar]*MCPar_py[nMCPar])+ 
	  //						 (MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 
	  //	      nGenPhotCand++; 
	  //             } 

	}
      nMCPar++;
    }
  }

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByLabel("photons",photons);

  edm::Handle<reco::MuonCollection> muons;
  //  iEvent.getByLabel("paramMuons", "ParamGlobalMuons", muons);
  iEvent.getByLabel("muons", muons);
  reco::MuonCollection::const_iterator muon;
  nMuonCand=0;
  for( muon = muons->begin(); muon != muons->end() && nMuonCand<MUONMAX; ++ muon ) {
     MuonCand_p[nMuonCand]=muon->p();
     MuonCand_px[nMuonCand]=muon->px();
     MuonCand_py[nMuonCand]=muon->py();
     MuonCand_pz[nMuonCand]=muon->pz();
     MuonCand_pt[nMuonCand]=muon->pt();
     MuonCand_eta[nMuonCand]=muon->eta();
     MuonCand_phi[nMuonCand]=muon->phi();
     MuonCand_charge[nMuonCand]=muon->charge();

     MuonCandTrack_p[nMuonCand] = muon->track()->p();

     MuonCand_vtxx[nMuonCand]=muon->vertex().x(); 
     MuonCand_vtxy[nMuonCand]=muon->vertex().y();  
     MuonCand_vtxz[nMuonCand]=muon->vertex().z();  

     // Isolation
     //     MuonCand_hcalisor3[nMuonCand]=muon->getIsolationR03().hadEt;
     //     MuonCand_ecalisor3[nMuonCand]=muon->getIsolationR03().emEt; 
     //     MuonCand_trkisor3[nMuonCand]=muon->getIsolationR03().nTracks; 
     //     MuonCand_hcalisor5[nMuonCand]=muon->getIsolationR05().hadEt; 
     //     MuonCand_ecalisor5[nMuonCand]=muon->getIsolationR05().emEt;  
     //     MuonCand_trkisor5[nMuonCand]=muon->getIsolationR05().nTracks;  


     if(MuonCand_pt[nMuonCand] > highestptmuon)
       highestptmuon = MuonCand_pt[nMuonCand];

     //     cout << MuonCand_p[nMuonCand] << ", " << MuonCand_charge[nMuonCand] << ", "  << MuonCand_eta[nMuonCand] << endl;
     nMuonCand++;
   }

  //2c: jets
  Handle<reco::CaloJetCollection> jets;
  iEvent.getByLabel("sisCone5CaloJets",jets);
  reco::CaloJetCollection::const_iterator jet;
  nJetCand=0;
  nJetCandE20=0;
  for( jet = jets->begin(); jet != jets->end() && nJetCand<JETMAX; ++ jet ) {
    JetCand_e[nJetCand]=jet->energy();
    if(JetCand_e[nJetCand] > 20.0)
      nJetCandE20++;
    JetCand_px[nJetCand]=jet->px();
    JetCand_py[nJetCand]=jet->py();
    JetCand_pz[nJetCand]=jet->pz();
    JetCand_phi[nJetCand]=jet->phi();
    JetCand_eta[nJetCand]=jet->eta();
    nJetCand++;
  }
  
  //2c: tracks and vertexing
  Handle<reco::TrackCollection> tracks; 
  iEvent.getByLabel("generalTracks",tracks); 
  reco::TrackCollection::const_iterator track; 
  nTrackCand=0; 

  edm::ESHandle<TransientTrackBuilder> theVtx;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theVtx);
  //  vector < reco::TransientTrack > mutrks;
  vector<TransientTrack> transmutrks; 
  reco::TrackCollection * mutrks = new reco::TrackCollection;


  // First get muon tracks
  bool isMuon = false;
  for( track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track ) 
    { 
      isMuon = false;
      for(int j = 0;j < nMuonCand; j++)
	{
	  if(MuonCandTrack_p[j] == track->p())
	    {
	      isMuon = true;
	      mutrks->push_back( *track );
	      TransientTrack tmptrk = (*theVtx).build( *track );
	      transmutrks.push_back( tmptrk );
	    }
	}
    }

  // If 2 muons, make a vertex
  if(transmutrks.size() > 1) 
    { 
      KalmanVertexFitter fitter(true); 
      TransientVertex mumuVertex = fitter.vertex(transmutrks); 
      MuMu_vtxx = mumuVertex.position().x(); 
      MuMu_vtxy = mumuVertex.position().y(); 
      MuMu_vtxz = mumuVertex.position().z(); 
    } 
  else 
    { 
      MuMu_vtxx = 0; 
      MuMu_vtxy = 0; 
      MuMu_vtxz = 0; 
    } 

  // Now go back and look at other tracks
  isMuon = false;
  for( track = tracks->begin(); track != tracks->end() && nTrackCand<TRACKMAX; ++ track )  
    {  
      isMuon = false; 
      for(int j = 0;j < nMuonCand; j++) 
        { 
          if(MuonCandTrack_p[j] == track->p()) 
            { 
              isMuon = true; 
	    }
	}

      if(isMuon == false)
	{
	  TrackCand_p[nTrackCand]=track->p();
	  TrackCand_pt[nTrackCand]=track->pt();  
	  TrackCand_px[nTrackCand]=track->px(); 
	  TrackCand_py[nTrackCand]=track->py(); 
	  TrackCand_pz[nTrackCand]=track->pz(); 
	  TrackCand_phi[nTrackCand]=track->phi(); 
	  TrackCand_eta[nTrackCand]=track->eta(); 
	  TrackCand_vtxx[nTrackCand]=track->vertex().x();
	  TrackCand_vtxy[nTrackCand]=track->vertex().y(); 
	  TrackCand_vtxz[nTrackCand]=track->vertex().z(); 

	  if(fabs(TrackCand_vtxz[nTrackCand] - MuMu_vtxz) < 1.0)
	    MuMu_extratracksz10mm++;
          if(fabs(TrackCand_vtxz[nTrackCand] - MuMu_vtxz) < 0.5) 
            MuMu_extratracksz5mm++; 
          if(fabs(TrackCand_vtxz[nTrackCand] - MuMu_vtxz) < 0.2) 
            MuMu_extratracksz2mm++; 
          if(fabs(TrackCand_vtxz[nTrackCand] - MuMu_vtxz) < 0.1) 
            MuMu_extratracksz1mm++; 
	  
	  nTrackCand++; 
	}
    }

   // calculate something and fill a histogram:
   for(int imu=0; imu<nMuonCand; imu++){
     for(int jmu=imu+1; jmu<nMuonCand; jmu++){
       if(MuonCand_charge[imu]*MuonCand_charge[jmu]<0){//opposite charge muons
	 double mass = pow(MuonCand_p[imu]+MuonCand_p[jmu],2);
	 mass-=pow(MuonCand_px[imu]+MuonCand_px[jmu],2);
	 mass-=pow(MuonCand_py[imu]+MuonCand_py[jmu],2);
	 mass-=pow(MuonCand_pz[imu]+MuonCand_pz[jmu],2);
	 MuMu_mass = sqrt(mass);
	 //	 cout << "MuMu_mass = " << MuMu_mass << endl;
       }
     }
   }

   // calculate something and fill a histogram:
   for(int imu=0; imu<nGenMuonCand; imu++){
     for(int jmu=imu+1; jmu<nGenMuonCand; jmu++){
       //       if(GenMuonCand_charge[imu]*GenMuonCand_charge[jmu]<0){//opposite charge muons
       double mass = pow(GenMuonCand_p[imu]+GenMuonCand_p[jmu],2);
       mass-=pow(GenMuonCand_px[imu]+GenMuonCand_px[jmu],2);
       mass-=pow(GenMuonCand_py[imu]+GenMuonCand_py[jmu],2);
       mass-=pow(GenMuonCand_pz[imu]+GenMuonCand_pz[jmu],2);
       GenMuMu_mass = sqrt(mass);
       //	 std::cout << "Gen mumu mass = " << GenMuMu_mass << std::endl;
     }
   }

   // calculate something and fill a histogram:
   //   for(int imu=0; imu<nGenProtCand; imu++){
   //     for(int jmu=imu+1; jmu<nGenProtCand; jmu++){
   //       double mass = pow(GenProtCand_p[imu]+GenProtCand_p[jmu],2);
   //       mass-=pow(GenProtCand_px[imu]+GenProtCand_px[jmu],2);
   //       mass-=pow(GenProtCand_py[imu]+GenProtCand_py[jmu],2);
   //       mass-=pow(GenProtCand_pz[imu]+GenProtCand_pz[jmu],2);
   //       GenProPro_mass = sqrt(mass);
   //     }
   //} 
 
   if(nMuonCand == 2)
     {
       // psuedo trigger for FastSim
       //       if(highestptmuon > 19.0)
       //	 {
       //	   cout << "Muon passed" << endl;
	   smalltree->Fill();
       //	 }
     }
}

