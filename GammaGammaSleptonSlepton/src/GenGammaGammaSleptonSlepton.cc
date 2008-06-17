
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

#include "UserCode/GammaGammaSleptonSlepton/interface/GenGammaGammaSleptonSlepton.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>


using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;


GenGammaGammaSleptonSlepton::GenGammaGammaSleptonSlepton(const edm::ParameterSet& iConfig)
{
  nEvt=0;
  MCPARMAX=1000;
  PROTMAX=100;
  MUONMAX = 10;
  ELEMAX = 10;

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

  smalltree->Branch("nGenMuonCand",&nGenMuonCand,"nGenMuonCand/I");
  smalltree->Branch("GenMuonCand_px",GenMuonCand_px,"GenMuonCand_px[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_py",GenMuonCand_py,"GenMuonCand_py[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_pz",GenMuonCand_pz,"GenMuonCand_pz[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_p",GenMuonCand_p,"GenMuonCand_p[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_pt",GenMuonCand_pt,"GenMuonCand_pt[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_eta",GenMuonCand_eta,"GenMuonCand_eta[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_phi",GenMuonCand_phi,"GenMuonCand_phi[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_e",GenMuonCand_e,"GenMuonCand_e[nGenMuonCand]/D");
  smalltree->Branch("GenMuonCand_moth",GenMuonCand_moth,"GenMuonCand_moth[nGenMuonCand]/I"); 

  smalltree->Branch("nGenEleCand",&nGenEleCand,"nGenEleCand/I"); 
  smalltree->Branch("GenEleCand_px",GenEleCand_px,"GenEleCand_px[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_py",GenEleCand_py,"GenEleCand_py[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_pz",GenEleCand_pz,"GenEleCand_pz[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_p",GenEleCand_p,"GenEleCand_p[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_pt",GenEleCand_pt,"GenEleCand_pt[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_eta",GenEleCand_eta,"GenEleCand_eta[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_phi",GenEleCand_phi,"GenEleCand_phi[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_e",GenEleCand_e,"GenEleCand_e[nGenEleCand]/D"); 
  smalltree->Branch("GenEleCand_moth",GenEleCand_moth,"GenEleCand_moth[nGenEleCand]/I");

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
  smalltree->Branch("GenMuMu_dphi",&GenMuMu_dphi,"GenMuMu_dphi/D");
  smalltree->Branch("GenMuMu_dpt",&GenMuMu_dpt,"GenMuMu_dpt/D");
  smalltree->Branch("GenElEl_mass",&GenElEl_mass,"GenElEl_mass/D");
  smalltree->Branch("GenElEl_dphi",&GenElEl_dphi,"GenElEl_dphi/D");
  smalltree->Branch("GenElEl_dpt",&GenElEl_dpt,"GenElEl_dpt/D");
  smalltree->Branch("GenProPro_mass",&GenProPro_mass,"GenProPro_mass/D");

}


GenGammaGammaSleptonSlepton::~GenGammaGammaSleptonSlepton() 
{  
  thefile->Write();
  thefile->Close();
}

void GenGammaGammaSleptonSlepton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // this analyzer produces a small root file with basic candidates and some MC information
  // some additional print statements
  nEvt++;
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;

  // step 1: fill some basic MC information into the root tree
  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel("source", evt);
 
  nMCPar=0;
  nGenEleCand = 0;
  nGenMuonCand = 0;
  nGenProtCand = 0;
  nGenSmuonCand = 0;
  nGenNeutCand = 0;
  nGenPhotCand = 0;

  GenMuMu_mass = -1.0; GenMuMu_dphi = -1.0; GenMuMu_dpt = -1.0;
  GenElEl_mass = -1.0; GenElEl_dphi = -1.0; GenElEl_dpt = -1.0;

  double highesteproton = 0.0;

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

       MCPar_motherid[nMCPar] = -1;

      if((*p)->status() < 3)
	{
	  //	  HepMC::GenParticle* mother = (*((*p)->production_vertex()->particles_in_const_begin()));  
	  //	  MCPar_motherid[nMCPar] = mother->pdg_id();   

	  if((((*p)->pdg_id() == 13) || ((*p)->pdg_id() == -13)) && (nGenMuonCand < MUONMAX))
	    {
	      GenMuonCand_px[nGenMuonCand] = MCPar_px[nMCPar];
	      GenMuonCand_py[nGenMuonCand] = MCPar_py[nMCPar];
	      GenMuonCand_pz[nGenMuonCand] = MCPar_pz[nMCPar];
	      GenMuonCand_phi[nGenMuonCand] = MCPar_phi[nMCPar];
	      GenMuonCand_pt[nGenMuonCand] = sqrt((GenMuonCand_px[nGenMuonCand] * GenMuonCand_px[nGenMuonCand]) +
						  (GenMuonCand_py[nGenMuonCand] * GenMuonCand_py[nGenMuonCand]));
	      GenMuonCand_charge[nGenMuonCand] = 
		(int)(MCPar_pdgid[nMCPar]/13);
	      GenMuonCand_eta[nGenMuonCand] = MCPar_eta[nMCPar];
	      GenMuonCand_e[nGenMuonCand] = MCPar_e[nMCPar];
	      GenMuonCand_p[nGenMuonCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+ 
						 (MCPar_py[nMCPar]*MCPar_py[nMCPar])+
						 (MCPar_pz[nMCPar]*MCPar_pz[nMCPar]));
	      GenMuonCand_moth[nGenMuonCand] = MCPar_motherid[nMCPar];

	      nGenMuonCand++;
	    }

	  if((((*p)->pdg_id() == 11) || ((*p)->pdg_id() == -11)) && (nGenEleCand < ELEMAX)) 
	    { 
	      GenEleCand_px[nGenEleCand] = MCPar_px[nMCPar]; 
	      GenEleCand_py[nGenEleCand] = MCPar_py[nMCPar]; 
	      GenEleCand_pz[nGenEleCand] = MCPar_pz[nMCPar]; 
	      GenEleCand_phi[nGenEleCand] = MCPar_phi[nMCPar]; 
	      GenEleCand_pt[nGenEleCand] = sqrt((GenEleCand_px[nGenEleCand] * GenEleCand_px[nGenEleCand]) + 
						(GenEleCand_py[nGenEleCand] * GenEleCand_py[nGenEleCand])); 
	      GenEleCand_charge[nGenEleCand] =  
		(int)(MCPar_pdgid[nMCPar]/13); 
	      GenEleCand_eta[nGenEleCand] = MCPar_eta[nMCPar]; 
	      GenEleCand_e[nGenEleCand] = MCPar_e[nMCPar]; 
	      GenEleCand_p[nGenEleCand] = sqrt((MCPar_px[nMCPar]*MCPar_px[nMCPar])+  
					       (MCPar_py[nMCPar]*MCPar_py[nMCPar])+ 
					       (MCPar_pz[nMCPar]*MCPar_pz[nMCPar])); 
              GenEleCand_moth[nGenEleCand] = MCPar_motherid[nMCPar]; 
	      
	      nGenEleCand++; 
	    } 

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

       double dphi = abs(GenMuonCand_phi[0]-GenMuonCand_phi[1]); 
       if(dphi > 3.14159) 
	 dphi = (2.0*3.14159) - dphi;        
       GenMuMu_dphi = dphi; 

       GenMuMu_dpt = GenMuonCand_pt[0]-GenMuonCand_pt[1];

     }
   }
   
   for(int imu=0; imu<nGenEleCand; imu++){ 
     for(int jmu=imu+1; jmu<nGenEleCand; jmu++){ 
       //       if(GenEleCand_charge[imu]*GenEleCand_charge[jmu]<0){//opposite charge muons 
       double mass = pow(GenEleCand_p[imu]+GenEleCand_p[jmu],2); 
       mass-=pow(GenEleCand_px[imu]+GenEleCand_px[jmu],2); 
       mass-=pow(GenEleCand_py[imu]+GenEleCand_py[jmu],2); 
       mass-=pow(GenEleCand_pz[imu]+GenEleCand_pz[jmu],2); 
       GenElEl_mass = sqrt(mass); 
       //        std::cout << "Gen mumu mass = " << GenMuMu_mass << std::endl; 

       double dphi = abs(GenEleCand_phi[0]-GenEleCand_phi[1]); 
       if(dphi > 3.14159) 
	 dphi = (2.0*3.14159) - dphi;        
       GenElEl_dphi = dphi; 

       GenElEl_dpt = GenEleCand_pt[0] - GenEleCand_pt[1];
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
 
   if(nGenMuonCand == 2 || nGenEleCand == 2)
     {
	   smalltree->Fill();
     }
}

