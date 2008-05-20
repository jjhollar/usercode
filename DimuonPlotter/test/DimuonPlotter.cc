
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

  TH1F *m_mumu;
  TH1F *m_mutrk;
  TH1F *m_ups_mumu;
  TH1F *m_ups_mutrk;
  TH1F *m_jpsi_mumu;
  TH1F *m_jpsi_mutrk;
  TH1F *m_z_mumu;
  TH1F *m_z_mutrk;
  TH1F *triggerbits;

  edm::TriggerNames trigNames ;

  int nDimuonCand, nTrkDimuonCand, nSADimuonCand;
  int nTrig;
  int DIMUONMAX;// used to set maximum of arrays
  double DimuonCand_goodmumu_mass[10];
  double DimuonCand_mumuonetrack_mass[10]; 
  double DimuonCand_mumuonesamuon_mass[10]; 
  double MuonCand_pt1[10], MuonCand_pt2[10];
  double MuonCand_eta1[10], MuonCand_eta2[10];

};


DimuonPlotter::DimuonPlotter(const edm::ParameterSet& iConfig)
{
  nEvt=0;
  DIMUONMAX=10;
  nTrig = 4;

  thefile = new TFile("dimuon.plot.root","recreate");
  thefile->cd();

  m_mumu = new TH1F("m_mumu","m_mumu",200,0,200);
  m_mutrk = new TH1F("m_mutrk","m_mutrk",200,0,200); 
  m_jpsi_mumu = new TH1F("m_jpsi_mumu","m_jpsi_mumu",20,1,5); 
  m_jpsi_mutrk = new TH1F("m_jpsi_mutrk","m_jpsi_mutrk",20,1,5);  
  m_ups_mumu = new TH1F("m_ups_mumu","m_ups_mumu",40,7,11); 
  m_ups_mutrk = new TH1F("m_ups_mutrk","m_ups_mutrk",40,7,11);  
  m_z_mumu = new TH1F("m_z_mumu","m_z_mumu",100,60,120);  
  m_z_mutrk = new TH1F("m_z_mutrk","m_z_mutrk",100,60,120);   
  triggerbits = new TH1F("triggerbits","triggerbits",4,0,4);
  triggerbits->GetXaxis()->SetBinLabel(1,"HLT1MuonPrescalePt5");
  triggerbits->GetXaxis()->SetBinLabel(2,"HLT1MuonPrescalePt7x7"); 
  triggerbits->GetXaxis()->SetBinLabel(3,"HLT1MuonIso"); 
  triggerbits->GetXaxis()->SetBinLabel(4,"HLT2MuonNonIso"); 

}


DimuonPlotter::~DimuonPlotter() 
{  
  m_mumu->Write();
  m_mutrk->Write(); 
  m_jpsi_mumu->Write(); 
  m_jpsi_mutrk->Write();  
  m_ups_mumu->Write();  
  m_ups_mutrk->Write();   
  m_z_mumu->Write();  
  m_z_mutrk->Write();   
  triggerbits->Write();

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
	    triggerbits->Fill(0);
	}
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt7x7" )
	{
          if ( hltResults->accept(i) ) 
	    triggerbits->Fill(1);
	}
      if ( trigNames.triggerNames().at(i) == "HLT1MuonIso" )
        { 
          if ( hltResults->accept(i) ) 
	    triggerbits->Fill(2);
        } 
      if ( trigNames.triggerNames().at(i) == "HLT2MuonNonIso" )
        {  
          if ( hltResults->accept(i) ) 
	    triggerbits->Fill(3);
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

     m_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
     m_jpsi_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]); 
     m_ups_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);  
     m_z_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);  

     nDimuonCand++;
   }

   Handle<reco::CompositeCandidateCollection> trkdimuons; 
   iEvent.getByLabel("dimuonsOneTrack",trkdimuons); 
   reco::CompositeCandidateCollection::const_iterator trkdimuon; 
   nTrkDimuonCand=0; 
   for( trkdimuon = trkdimuons->begin(); trkdimuon != trkdimuons->end() && nTrkDimuonCand<DIMUONMAX; ++ trkdimuon ) { 
     DimuonCand_mumuonetrack_mass[nTrkDimuonCand]=trkdimuon->mass(); 

     m_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
     m_jpsi_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
     m_ups_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]); 
     m_z_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]); 

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

 
}


//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(DimuonPlotter);
