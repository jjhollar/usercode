
// -*- C++ -*-
//
// Package:    GammaGammaSleptonSlepton
// Class:      GammaGammaSleptonSlepton
// 

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
#include "DataFormats/TrackReco/interface/Track.h"   
#include "DataFormats/TrackReco/interface/TrackFwd.h"  
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"    
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/FP420Cluster/interface/RecoFP420.h" 
#include "DataFormats/FP420Cluster/interface/RecoCollectionFP420.h"
#include "DataFormats/FP420Cluster/interface/TrackFP420.h"  
#include "DataFormats/FP420Cluster/interface/TrackCollectionFP420.h" 
#include "DataFormats/FP420Cluster/interface/ClusterCollectionFP420.h" 
#include "DataFormats/FP420Cluster/interface/ClusterFP420.h"  

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>


using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;

class GammaGammaSleptonSlepton : public edm::EDAnalyzer {
public:
  explicit GammaGammaSleptonSlepton(const edm::ParameterSet&);
  ~GammaGammaSleptonSlepton();
  
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
private:
  int nEvt;// used to count the number of events

  
  // to be used for root output tree
  TFile *thefile;
  TTree *smalltree;
  int nMCPar;
  int MCPARMAX;// used to set maximum of arrays
  int MCPar_status[50];
  double MCPar_px[50];
  double MCPar_py[50];
  double MCPar_pz[50];
  double MCPar_e[50];
  double MCPar_phi[50];
  double MCPar_eta[50];
  double MCPar_mass[50];
  int MCPar_pdgid[50];
  int nEleCand;
  int ELEMAX;// used to set maximum of arrays
  double EleCand_px[10];
  double EleCand_py[10];
  double EleCand_pz[10];
  double EleCand_e[10];
  int nMuonCand;
  int MUONMAX;// used to set maximum of arrays
  double MuonCand_px[10];
  double MuonCand_py[10];
  double MuonCand_pz[10];
  double MuonCand_p[10];
  double MuonCand_eta[10];
  double MuonCand_pt[10];
  double MuonCand_phi[10];
  double MuonCand_e[10];
  int MuonCand_charge[10];
  int nJetCand;
  int JETMAX;// used to set maximum of arrays
  double JetCand_px[30];
  double JetCand_py[30];
  double JetCand_pz[30];
  double JetCand_e[30];
  double JetCand_eta[30];
  double JetCand_phi[30];
  double eventWeight;
  // a histogram

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
  int nGenProtCand;
  double GenProtCand_px[10];
  double GenProtCand_py[10];
  double GenProtCand_pz[10];
  double GenProtCand_p[10];
  double GenProtCand_eta[10];
  double GenProtCand_pt[10];
  double GenProtCand_phi[10];
  double GenProtCand_e[10];
  int GenProtCand_charge[10];
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

  int nProtCand;
  double ProtCand_e[10];
  double ProtCand_x[10];
  double ProtCand_y[10];
  int ProtCand_direction[10];

  int nGenPhotCand; 
  double GenPhotCand_px[10]; 
  double GenPhotCand_py[10]; 
  double GenPhotCand_pz[10]; 
  double GenPhotCand_p[10]; 
  double GenPhotCand_eta[10]; 
  double GenPhotCand_pt[10]; 
  double GenPhotCand_phi[10]; 
  double GenPhotCand_e[10]; 

  
  double GenMuMu_mass;
  double MuMu_mass;
  double GenProPro_mass;
};


GammaGammaSleptonSlepton::GammaGammaSleptonSlepton(const edm::ParameterSet& iConfig)
{
  nEvt=0;
  MCPARMAX=50;
  ELEMAX=10;
  MUONMAX=10;
  JETMAX=30;
  
  thefile = new TFile("gammagamma.anal.root","recreate");
  thefile->cd();

  smalltree= new TTree("ntp1","ntp1");
  smalltree->Branch("nMCPar",&nMCPar,"nMCPar/I");
  smalltree->Branch("MCPar_status",&MCPar_status,"MCPar_status[nMCPar]/I");
  smalltree->Branch("MCPar_px",MCPar_px,"MCPar_px[nMCPar]/D");
  smalltree->Branch("MCPar_py",MCPar_py,"MCPar_py[nMCPar]/D");
  smalltree->Branch("MCPar_pz",MCPar_pz,"MCPar_pz[nMCPar]/D");
  smalltree->Branch("MCPar_e",MCPar_e,"MCPar_e[nMCPar]/D");
  smalltree->Branch("MCPar_phi",MCPar_phi,"MCPar_phi[nMCPar]/D");
  smalltree->Branch("MCPar_eta",MCPar_eta,"MCPar_eta[nMCPar]/D");
  smalltree->Branch("MCPar_pdgid",MCPar_pdgid,"MCPar_pdgid[nMCPar]/I");
  smalltree->Branch("MCPar_mass",MCPar_mass,"MCPar_mass[nMCPar]/D");

//   smalltree->Branch("nEleCand",&nEleCand,"nEleCand/I");
//   smalltree->Branch("EleCand_px",EleCand_px,"EleCand_px[nEleCand]/D");
//   smalltree->Branch("EleCand_py",EleCand_py,"EleCand_py[nEleCand]/D");
//   smalltree->Branch("EleCand_pz",EleCand_pz,"EleCand_pz[nEleCand]/D");
//   smalltree->Branch("EleCand_e",EleCand_e,"EleCand_e[nEleCand]/D");

  smalltree->Branch("nJetCand",&nJetCand,"nJetCand/I");
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

  smalltree->Branch("nGenPhotCand",&nGenPhotCand,"nGenPhotCand/I"); 
  smalltree->Branch("GenPhotCand_px",GenPhotCand_px,"GenPhotCand_px[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_py",GenPhotCand_py,"GenPhotCand_py[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_pz",GenPhotCand_pz,"GenPhotCand_pz[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_p",GenPhotCand_p,"GenPhotCand_p[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_pt",GenPhotCand_pt,"GenPhotCand_pt[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_eta",GenPhotCand_eta,"GenPhotCand_eta[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_phi",GenPhotCand_phi,"GenPhotCand_phi[nGenPhotCand]/D"); 
  smalltree->Branch("GenPhotCand_e",GenPhotCand_e,"GenPhotCand_e[nGenPhotCand]/D"); 
  
  smalltree->Branch("nProtCand",&nProtCand,"nProtCand/I");
  smalltree->Branch("ProtCand_e",ProtCand_e,"ProtCand_e[nProtCand]/D"); 
  smalltree->Branch("ProtCand_x",ProtCand_x,"ProtCand_x[nProtCand]/D"); 
  smalltree->Branch("ProtCand_y",ProtCand_y,"ProtCand_y[nProtCand]/D"); 
  smalltree->Branch("ProtCand_direction",ProtCand_direction,"ProtCand_direction[nProtCand]/I");

  smalltree->Branch("GenMuMu_mass",&GenMuMu_mass,"GenMuMu_mass/D");
  smalltree->Branch("MuMu_mass",&MuMu_mass,"MuMu_mass/D");
  smalltree->Branch("GenProPro_mass",&GenProPro_mass,"GenProPro_mass/D");
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
  
  
  // step 1: fill some basic MC information into the root tree
  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel("source", evt);
  
  //   std::cout << "Got Event" << std::endl;
  
  nMCPar=0;
  nGenMuonCand = 0;
  nGenProtCand = 0;
  nGenSmuonCand = 0;
  nGenNeutCand = 0;
  nGenPhotCand = 0;
  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(evt->GetEvent()));
  
  //   std::cout << "Got GenEvent" << std::endl;
  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end() && nMCPar<MCPARMAX ; ++p ) {
    if ( abs((*p)->pdg_id()) !=0 && (abs((*p)->pdg_id())<30 || abs((*p)->pdg_id() == 2212) || abs((*p)->pdg_id()) > 1000000) && nMCPar<MCPARMAX){
      MCPar_status[nMCPar]=(*p)->status();
      MCPar_px[nMCPar]=(*p)->momentum().x();
      MCPar_py[nMCPar]=(*p)->momentum().y();
      MCPar_pz[nMCPar]=(*p)->momentum().z();
      //       MCPar_e[nMCPar]=(*p)->momentum().e();//(*p)->energy() does not work!!;
      MCPar_phi[nMCPar]=atan2(MCPar_py[nMCPar],MCPar_px[nMCPar]);
      MCPar_eta[nMCPar] = (*p)->momentum().eta();
      MCPar_pdgid[nMCPar]=(*p)->pdg_id();
      MCPar_mass[nMCPar]=(*p)->momentum().m();
      MCPar_e[nMCPar] = sqrt(MCPar_mass[nMCPar]*MCPar_mass[nMCPar] + (MCPar_px[nMCPar]*MCPar_px[nMCPar] + MCPar_py[nMCPar]*MCPar_py[nMCPar] + MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 
       
      if((*p)->status() < 3)
	{
	  if(((*p)->pdg_id() == 13) || ((*p)->pdg_id() == -13))
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
	  if((*p)->pdg_id() == 2212)
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
	  if((*p)->pdg_id() == 22) 
	    { 
	      GenPhotCand_px[nGenPhotCand] = MCPar_px[nMCPar]; 
	      GenPhotCand_py[nGenPhotCand] = MCPar_py[nMCPar]; 
	      GenPhotCand_pz[nGenPhotCand] = MCPar_pz[nMCPar]; 
	      GenPhotCand_eta[nGenPhotCand] = MCPar_eta[nMCPar]; 
	      GenPhotCand_phi[nGenPhotCand] = MCPar_phi[nMCPar]; 
	      GenPhotCand_eta[nGenPhotCand] = MCPar_eta[nMCPar]; 
	      GenPhotCand_e[nGenPhotCand] = MCPar_e[nMCPar]; 
	      GenPhotCand_p[nGenPhotCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+  
						 (MCPar_py[nMCPar]*MCPar_py[nMCPar])+ 
						 (MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 
	      nGenPhotCand++; 
             } 
	  
	}
      nMCPar++;
    }
  }

  nProtCand = 0;
  Handle<RecoCollectionFP420> recocollSrc; 
  iEvent.getByLabel("FP420Reco" , recocollSrc); 
  int dn0=3;  
  for (int number_detunits=1; number_detunits<dn0; number_detunits++) 
    { 
      int StID = number_detunits; // =1 for +420m set up , =2 for -420m set up 
      
      std::vector<RecoFP420> zcollector; 
      zcollector.clear(); 
      
      RecoCollectionFP420::Range toutputRange; 
      toutputRange = recocollSrc->get(StID); 
      
      // fill output in zcollector vector (for may be sorting? or other checks) 
      RecoCollectionFP420::ContainerIterator sort_begin = toutputRange.first; 
      RecoCollectionFP420::ContainerIterator sort_end = toutputRange.second; 

      for ( ;sort_begin != sort_end; ++sort_begin )
	{ 
	  zcollector.push_back(*sort_begin); 
	} 

      vector<RecoFP420>::const_iterator simHitIter = zcollector.begin(); 
      vector<RecoFP420>::const_iterator simHitIterEnd = zcollector.end(); 
      // loop in #tracks 
      for (;simHitIter != simHitIterEnd; ++simHitIter) 
	{ 
	  const RecoFP420 itrack = *simHitIter; 
	  
	  std::cout << "   itrack.direction():    " << itrack.direction() << std::endl; 
	  std::cout << "   itrack.e0():    " << itrack.e0() << std::endl; 
	  std::cout << "   itrack.x0():    " << itrack.x0() << std::endl; 
	  std::cout << "   itrack.y0():    " << itrack.y0() << std::endl; 
	  std::cout << "   itrack.tx0():    " << itrack.tx0() << std::endl; 
	  std::cout << "   itrack.ty0():    " << itrack.ty0() << std::endl; 
	  std::cout << "   itrack.q20():    " << itrack.q20() << std::endl; 

	  ProtCand_e[nProtCand] = itrack.e0();
	  ProtCand_x[nProtCand] = itrack.x0(); 
          ProtCand_y[nProtCand] = itrack.y0(); 
	  ProtCand_direction[nProtCand] = itrack.direction();

	  nProtCand++;

	  long double tx,ty,th,pfi,pfigrad,eta1= 0,recovtxX,recovtxY; 

	  recovtxX=itrack.x0()/1000.; recovtxY=itrack.y0()/1000.; // goes to mm 
	  tx = itrack.tx0() / 1000000.; ty = itrack.ty0() / 1000000.; // goes to radians 
	  th  = std::sqrt((tx*tx) + (ty*ty)); 
	  pfi     = std::atan2(tx,ty); if (pfi < 0.) pfi += (2.0 * 3.14159); 
	  eta1 = -std::log(std::tan(th/2)); 

	  long double Recop1 = 7000.-itrack.e0(); 
	  double Recoxi1 = 1.-Recop1/7000.; 
	  // cross-checks: 
	  long double ccc = (1.-std::cos(th)); 
	  long double Recot = 2.*7000.*Recop1*ccc; 
	  long double Recop2 = 7000.; 
	  if(ccc !=0.) Recop2 = itrack.q20()/(2.*7000.*ccc); 
	  double ddddddp = Recop1-Recop2; 
	  double ddddddt = itrack.q20()-Recot; 
	  double Recoxi2 = 1.-Recop2/7000.; 
	}// track loop 

    }

   Handle<reco::MuonCollection> muons;
   iEvent.getByLabel("muons",muons);
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
     nMuonCand++;
   }
   //2c: jets
   Handle<reco::CaloJetCollection> jets;
   iEvent.getByLabel("iterativeCone5CaloJets",jets);
   reco::CaloJetCollection::const_iterator jet;
   nJetCand=0;
   for( jet = jets->begin(); jet != jets->end() && nJetCand<JETMAX; ++ jet ) {
     JetCand_e[nJetCand]=jet->energy();
     JetCand_px[nJetCand]=jet->px();
     JetCand_py[nJetCand]=jet->py();
     JetCand_pz[nJetCand]=jet->pz();
     JetCand_phi[nJetCand]=jet->phi();
     JetCand_eta[nJetCand]=jet->eta();
     nJetCand++;
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
//       }
     }
   }  
   // calculate something and fill a histogram:
   for(int imu=0; imu<nGenProtCand; imu++){
     for(int jmu=imu+1; jmu<nGenProtCand; jmu++){
	 double mass = pow(GenProtCand_p[imu]+GenProtCand_p[jmu],2);
	 mass-=pow(GenProtCand_px[imu]+GenProtCand_px[jmu],2);
	 mass-=pow(GenProtCand_py[imu]+GenProtCand_py[jmu],2);
	 mass-=pow(GenProtCand_pz[imu]+GenProtCand_pz[jmu],2);
	 GenProPro_mass = sqrt(mass);
     }
   }  
 
  smalltree->Fill();
}


//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(GammaGammaSleptonSlepton);
