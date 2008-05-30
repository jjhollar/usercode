
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
#include <TObject.h>
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
  double readefficiency(double pt, double eta);

private:
  int nEvt;// used to count the number of events

  TFile *theefffile;
  
  // to be used for root output tree
  TFile *thefile;

  TH1F *m_mumu;
  TH1F *m_mutrk;
  TH1F *m_musamu; 

  TH1F *m_ups_mumu;
  TH1F *m_ups_mutrk;
  TH1F *m_ups_musamu;

  TH1F *m_jpsi_mumu;
  TH1F *m_jpsi_mutrk;
  TH1F *m_jpsi_musamu;

  TH1F *m_z_mumu;
  TH1F *m_z_mutrk;
  TH1F *m_z_musamu;

  TH1F *m_mumu_effcor;

  TH1F *pt_mu;
  TH1F *pt_mu_effcor;
  TH1F *pt_low_mu;
  TH1F *pt_low_mu_effcor;
  TH1F *eta_mu;
  TH1F *eta_mu_effcor;
  TH1F *pt_mu_mutrk;
  TH1F *eta_mu_mutrk;

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

  int nbinspt;
  double ptlow;
  double pthigh;
  double binwidth;
};


DimuonPlotter::DimuonPlotter(const edm::ParameterSet& iConfig)
{
  nEvt=0;
  DIMUONMAX=10;
  nTrig = 4;

  nbinspt = 20; ptlow = 0.0; pthigh = 100.0; binwidth = (pthigh-ptlow)/nbinspt;

  theefffile = TFile::Open("muon_eff_test.root");
  thefile = new TFile("dimuon.plot.root","recreate");
  thefile->cd();

  m_mumu = new TH1F("m_mumu","m_mumu",500,0,200);
  m_mutrk = new TH1F("m_mutrk","m_mutrk",500,0,200); 
  m_musamu = new TH1F("m_musamu","m_musamu",500,0,200);

  m_jpsi_mumu = new TH1F("m_jpsi_mumu","m_jpsi_mumu",100,1,5); 
  m_jpsi_mutrk = new TH1F("m_jpsi_mutrk","m_jpsi_mutrk",100,1,5);  
  m_jpsi_musamu = new TH1F("m_jpsi_musamu","m_jpsi_musamu",100,1,5);   

  m_ups_mumu = new TH1F("m_ups_mumu","m_ups_mumu",100,7,12); 
  m_ups_mutrk = new TH1F("m_ups_mutrk","m_ups_mutrk",100,7,12);  
  m_ups_musamu = new TH1F("m_ups_musamu","m_ups_musamu",100,7,12);   

  m_z_mumu = new TH1F("m_z_mumu","m_z_mumu",100,60,120);  
  m_z_mutrk = new TH1F("m_z_mutrk","m_z_mutrk",100,60,120);   
  m_z_musamu = new TH1F("m_z_musamu","m_z_musamu",100,60,120);    

  eta_mu = new TH1F("eta_mu","eta_mu",50,-2.5,2.5);
  pt_mu = new TH1F("pt_mu","pt_mu",100,0,100);
  eta_mu_mutrk = new TH1F("eta_mu_mutrk","eta_mu_mutrk",50,-2.5,2.5); 
  pt_mu_mutrk = new TH1F("pt_mu_mutrk","pt_mu_mutrk",100,0,100); 
  m_mumu_effcor = new TH1F("m_mumu_effcor","m_mumu_effcor",500,0,200);
  pt_mu_effcor = new TH1F("pt_mu_effcor","pt_mu_effcor",100,0,100);
  eta_mu_effcor = new TH1F("eta_mu_effcor","eta_mu_effcor",50,-2.5,2.5); 
  pt_low_mu = new TH1F("pt_low_mu","pt_low_mu",20,0,20);
  pt_low_mu_effcor = new TH1F("pt_low_mu_effcor","pt_low_mu_effcor",20,0,20);

  triggerbits = new TH1F("triggerbits","triggerbits",6,0,6);
  triggerbits->GetXaxis()->SetBinLabel(1,"HLT1MuonPrescalePt3");
  triggerbits->GetXaxis()->SetBinLabel(2,"HLT1MuonPrescalePt5");
  triggerbits->GetXaxis()->SetBinLabel(3,"HLT1MuonPrescalePt7x7"); 
  triggerbits->GetXaxis()->SetBinLabel(4,"HLT1MuonIso"); 
  triggerbits->GetXaxis()->SetBinLabel(5,"HLT1MuonNonIso15");   
  triggerbits->GetXaxis()->SetBinLabel(6,"HLT2MuonNonIso"); 
}


DimuonPlotter::~DimuonPlotter() 
{  
  theefffile->Close();
  thefile->cd();

  m_mumu->Write();
  m_mutrk->Write(); 
  m_musamu->Write();
  
  m_jpsi_mumu->Write(); 
  m_jpsi_mutrk->Write();  
  m_jpsi_musamu->Write();

  m_ups_mumu->Write();  
  m_ups_mutrk->Write();   
  m_ups_musamu->Write();

  m_z_mumu->Write();  
  m_z_mutrk->Write();   
  m_z_musamu->Write();

  m_mumu_effcor->Write();
  triggerbits->Write();
  pt_mu->Write();
  eta_mu->Write();
  pt_mu_effcor->Write();
  eta_mu_effcor->Write();
  pt_low_mu->Write();
  pt_low_mu_effcor->Write();

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
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt3" )      
        { 
          if ( hltResults->accept(i) ) 
            triggerbits->Fill(0); 
        } 
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt5" )     
	{
	  if ( hltResults->accept(i) )
	    triggerbits->Fill(1);
	}
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt7x7" )
	{
          if ( hltResults->accept(i) ) 
	    triggerbits->Fill(2);
	}
      if ( trigNames.triggerNames().at(i) == "HLT1MuonIso" )
        { 
          if ( hltResults->accept(i) ) 
	    triggerbits->Fill(3);
        } 
      if ( trigNames.triggerNames().at(i) == "HLT1MuonNonIso15" ) 
        {  
          if ( hltResults->accept(i) )  
            triggerbits->Fill(4); 
        }  
      if ( trigNames.triggerNames().at(i) == "HLT2MuonNonIso" )
        {  
          if ( hltResults->accept(i) ) 
	    triggerbits->Fill(5);
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

     cout << "Candidate with mass = " << dimuon->mass() << endl;
     double mu1weight = readefficiency(dimuon->daughter(0)->pt(),dimuon->daughter(0)->eta());
     double mu2weight = readefficiency(dimuon->daughter(1)->pt(),dimuon->daughter(1)->eta()); 
     double dimuweight = (1.0/(mu1weight*mu2weight));
     cout << "Weight = " << (1.0/mu1weight) << ", " << (1.0/mu2weight) << ", " << dimuweight << endl;

     if(mu1weight == 0 || mu2weight == 0)
       {
	 mu1weight = 1.0;
	 mu2weight = 1.0;
	 dimuweight = 1.0;
       }

     pt_mu_effcor->Fill(MuonCand_pt1[nDimuonCand],1.0/mu1weight);
     pt_mu_effcor->Fill(MuonCand_pt2[nDimuonCand],1.0/mu2weight); 
     eta_mu_effcor->Fill(MuonCand_eta1[nDimuonCand],1.0/mu1weight); 
     eta_mu_effcor->Fill(MuonCand_eta2[nDimuonCand],1.0/mu2weight); 

     m_mumu_effcor->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight);
     m_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
     m_jpsi_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]); 
     m_ups_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);  
     m_z_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);  
     pt_mu->Fill(MuonCand_pt1[nDimuonCand]);
     pt_mu->Fill(MuonCand_pt2[nDimuonCand]);
     eta_mu->Fill(MuonCand_eta1[nDimuonCand]); 
     eta_mu->Fill(MuonCand_eta2[nDimuonCand]); 

     pt_low_mu->Fill(MuonCand_pt1[nDimuonCand]);
     pt_low_mu->Fill(MuonCand_pt2[nDimuonCand]); 
     pt_low_mu_effcor->Fill(MuonCand_pt1[nDimuonCand],1.0/mu1weight);
     pt_low_mu_effcor->Fill(MuonCand_pt2[nDimuonCand],1.0/mu2weight); 

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

     //     pt_mu_mutrk->Fill(dimuon->daughter(0)->pt()); 
     //     pt_mu_mutrk->Fill(dimuon->daughter(1)->pt()); 
     //     eta_mu_mutrk->Fill(dimuon->daughter(0)->eta());  
     //     eta_mu_mutrk->Fill(dimuon->daughter(1)->eta());  

     nTrkDimuonCand++; 
   } 

   //   Handle<reco::CompositeCandidateCollection> sadimuons;  
   //   iEvent.getByLabel("dimuonsOneStandAloneMuon",sadimuons);  
   //   reco::CompositeCandidateCollection::const_iterator sadimuon;  
   //   nSADimuonCand=0;  
   //   for( sadimuon = sadimuons->begin(); sadimuon != sadimuons->end() && nSADimuonCand<DIMUONMAX; ++ sadimuon ) {  
   //     DimuonCand_mumuonesamuon_mass[nSADimuonCand]=sadimuon->mass();  
   //
   //     m_musamu->Fill(DimuonCand_mumuonesamuon_mass[nSADimuonCand]); 
   //     m_jpsi_musamu->Fill(DimuonCand_mumuonesamuon_mass[nSADimuonCand]); 
   //     m_ups_musamu->Fill(DimuonCand_mumuonesamuon_mass[nSADimuonCand]);  
   //     m_z_musamu->Fill(DimuonCand_mumuonesamuon_mass[nSADimuonCand]);  
   //
   //
   //     nSADimuonCand++;  
   //   }   
}

double DimuonPlotter::readefficiency(double pt, double eta)
{
  double theeff = 0.0;

  //This is testing the readback of Tag&Probe efficiencies
  //  theefffile->cd();
  //  TH1F *hpteff = new TH1F();
  //  hpteff = (TH1F *)theefffile->Get("heff_Pt");
  //  double binhigh = 0.0;
  //  double binlow = 0.0;

  //  for(int i = 0;i < hpteff->GetNbinsX(); i++)
  //    {
  //      if(pt > (i * binwidth) && pt < ((i+1)*binwidth))
  //	theeff = hpteff->GetBinContent(i+1);
  //    }

  return theeff;
}

//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(DimuonPlotter);
