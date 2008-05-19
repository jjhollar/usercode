
// -*- C++ -*-
//
// Package:    DimuonPlotter
// Class:      DimuonPlotter
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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>


using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;

class DimuonPlotter : public edm::EDAnalyzer {
public:
  explicit DimuonPlotter(const edm::ParameterSet&);
  ~DimuonPlotter();
  
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
private:
  int nEvt;// used to count the number of events

  
  // to be used for root output tree
  TFile *thefile;
  TTree *smalltree;

  edm::TriggerNames trigNames ;

  int nDimuonCand, nTrkDimuonCand, nSADimuonCand;
  int nTrig;
  int DIMUONMAX;// used to set maximum of arrays
  double DimuonCand_goodmumu_mass[10];
  double DimuonCand_mumuonetrack_mass[10]; 
  double DimuonCand_mumuonesamuon_mass[10]; 
  double MuonCand_pt1[10], MuonCand_pt2[10];
  double MuonCand_eta1[10], MuonCand_eta2[10];
  int TriggerBits[10];
};


DimuonPlotter::DimuonPlotter(const edm::ParameterSet& iConfig)
{
  nEvt=0;
  DIMUONMAX=10;
  nTrig = 4;

  thefile = new TFile("dimuon.plot.root","recreate");
  thefile->cd();

  smalltree= new TTree("ntp1","ntp1");

  smalltree->Branch("nDimuonCand",&nDimuonCand,"nDimuonCand/I");
  smalltree->Branch("nTrkDimuonCand",&nTrkDimuonCand,"nTrkDimuonCand/I"); 
  //  smalltree->Branch("nSADimuonCand",&nSADimuonCand,"nSADimuonCand/I"); 

  smalltree->Branch("nTrig",&nTrig,"nTrig/I");
  smalltree->Branch("TriggerBits",TriggerBits,"TriggerBits[nTrig]/I");
  smalltree->Branch("DimuonCand_goodmumu_mass",DimuonCand_goodmumu_mass,"DimuonCand_goodmumu_mass[nDimuonCand]/D");
  smalltree->Branch("DimuonCand_mumuonetrack_mass",DimuonCand_mumuonetrack_mass,"DimuonCand_mumuonetrack_mass[nTrkDimuonCand]/D"); 
  smalltree->Branch("MuonCand_pt1",MuonCand_pt1,"MuonCand_pt[nDimuonCand]/D");
  smalltree->Branch("MuonCand_eta1",MuonCand_eta1,"MuonCand_eta[nDimuonCand]/D");
  smalltree->Branch("MuonCand_pt2",MuonCand_pt2,"MuonCand_pt2[nDimuonCand]/D"); 
  smalltree->Branch("MuonCand_eta2",&MuonCand_eta2,"MuonCand_eta2[nDimuonCand]/D"); 
  //  smalltree->Branch("DimuonCand_mumuonesamuon_mass",DimuonCand_mumuonesamuon_mass,"DimuonCand_mumuonesamuon_mass[nSADimuonCand]/D");  

}


DimuonPlotter::~DimuonPlotter() 
{  
  thefile->Write();
  thefile->Close();
}

void DimuonPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // this analyzer produces a small root file with basic candidates and some MC information
  // some additional print statements
  nEvt++;
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;
  
  //   std::cout << "Got Event" << std::endl;

  edm::Handle<edm::TriggerResults> hltResults ;
  //  edm::Handle<l1extra::L1ParticleMapCollection> l1Results ;

  iEvent.getByLabel(InputTag("TriggerResults::HLT"),hltResults) ;

  //  iEvent.getByLabel("l1extraParticleMap",l1Results) ;
  trigNames.init(*hltResults) ;

  //  cout << "# of triggers = " << trigNames.size() << endl;

  for (unsigned int i=0; i<trigNames.size(); i++) 
    {
      //      cout << "\tTrigger = " << trigNames.triggerNames().at(i) << endl;
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt5" )     
	{
	  if ( hltResults->accept(i) )
	    TriggerBits[0] = 1;
	  else
	    TriggerBits[0] = 0;
	}
      else if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt7x7" )
	{
          if ( hltResults->accept(i) ) 
            TriggerBits[1] = 1; 
          else 
            TriggerBits[1] = 0; 
	}
      else if ( trigNames.triggerNames().at(i) == "HLT1MuonIso" )
        { 
          if ( hltResults->accept(i) ) 
            TriggerBits[2] = 1; 
          else 
            TriggerBits[2] = 0; 
        } 
      else if ( trigNames.triggerNames().at(i) == "HLT2MuonNonIso" )
        {  
          if ( hltResults->accept(i) ) 
            TriggerBits[3] = 1; 
          else 
            TriggerBits[3] = 0; 
        }  


    }

   Handle<reco::CompositeCandidateCollection> dimuons;
   iEvent.getByLabel("dimuons",dimuons);
   reco::CompositeCandidateCollection::const_iterator dimuon;
   nDimuonCand=0;
   for( dimuon = dimuons->begin(); dimuon != dimuons->end() && nDimuonCand<DIMUONMAX; ++ dimuon ) {
     DimuonCand_goodmumu_mass[nDimuonCand]=dimuon->mass();
     MuonCand_pt1[nDimuonCand] = dimuon->daughter(0)->pt();
     MuonCand_pt2[nDimuonCand] = dimuon->daughter(1)->pt();
     MuonCand_eta1[nDimuonCand] = dimuon->daughter(0)->eta();
     MuonCand_eta2[nDimuonCand] = dimuon->daughter(1)->eta();
     nDimuonCand++;
   }

   Handle<reco::CompositeCandidateCollection> trkdimuons; 
   iEvent.getByLabel("dimuonsOneTrack",trkdimuons); 
   reco::CompositeCandidateCollection::const_iterator trkdimuon; 
   nTrkDimuonCand=0; 
   for( trkdimuon = trkdimuons->begin(); trkdimuon != trkdimuons->end() && nTrkDimuonCand<DIMUONMAX; ++ trkdimuon ) { 
     DimuonCand_mumuonetrack_mass[nTrkDimuonCand]=trkdimuon->mass(); 
     nTrkDimuonCand++; 
   } 

   //   Handle<reco::CompositeCandidateCollection> sadimuons;  
   //   iEvent.getByLabel("dimuonsOneStandAloneMuon",sadimuons);  
   //   reco::CompositeCandidateCollection::const_iterator sadimuon;  
   //   nSADimuonCand=0;  
   //   for( sadimuon = sadimuons->begin(); sadimuon != sadimuons->end() && nSADimuonCand<DIMUONMAX; ++ sadimuon ) {  
   //     DimuonCand_mumuonesamuon_mass[nSADimuonCand]=sadimuon->mass();  
   //     nSADimuonCand++;  
   //   }  

 
  smalltree->Fill();
}


//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(DimuonPlotter);
