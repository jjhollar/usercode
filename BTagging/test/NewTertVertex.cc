
// -*- C++ -*-
//
// Package:    NewTertVertex
// Class:      NewTertVertex
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

#include "RecoVertex/TertiaryTracksVertexFinder/interface/TertiaryTracksVertexFinder.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
//#include "RecoBTag/SoftLepton/plugins/findProductIDByLabel.h"
#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CondFormats/JetMETObjects/interface/SimpleMCJetCorrector.h"

#include "DataFormats/BTauReco/interface/BaseTagInfo.h"
#include "DataFormats/Common/interface/View.h"
#include "RecoBTag/MCTools/interface/AddTertiaryTracks.h"

#include "RecoVertex/TertiaryTracksVertexFinder/interface/TransientTrackInVertices.h"
#include "RecoVertex/TertiaryTracksVertexFinder/interface/VertexMass.h"

#include "RecoBTag/SecondaryVertex/interface/TrackSelector.h"
#include "RecoBTag/SecondaryVertex/interface/TrackSorting.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoBTag/SecondaryVertex/interface/VertexFilter.h"
#include "RecoBTag/SecondaryVertex/interface/VertexSorting.h"

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>


using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;

class NewTertVertex : public edm::EDAnalyzer {
public:
  explicit NewTertVertex(const edm::ParameterSet&);
  ~NewTertVertex();
  
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
private:
  TFile *thefile;
  TTree *thetree;

  InputTag primVtxLabel;
  InputTag jetTrackLabel;
  string outfile;
  JetFlavourIdentifier jfi;
  SimpleMCJetCorrector *jetcorrector;

  edm::ParameterSet               vtxRecoPSet;

  AddTertiaryTracks trackadder;

  double BBmass;

  double jetpx,jetpy,jetpz,jete,jeteta,jetphi,jetpt,jetcorr,jetpartonpx,jetpartonpy,jetpartonpz,jetpartonpt,jetpartoneta,jetpartonphi;
  int jetflavour, jetntracks, jetprimvtxntrk;

  double primvtxposx, primvtxposy, primvtxposz;

  int jetnsecvtx, jetnfreetrks;
  double jetsvbtag, jetsmbtag, jetjpbtag, jettrkbtag;

  double jetbestsvtxchi2, jetbestsvtxx, jetbestsvtxy, jetbestsvtxz, jetbestsvtxmass, jetbestsvtxdxy;
  int jetbestsvtxntrk;
  
  double jetfreetrksvtxmass, jetfreetrksvtxposx, jetfreetrksvtxposy, jetfreetrksvtxposz; 
  double jetfreetrksvtxchi2, jetfreetrksvtxrxy, jetfreetrksvtxrxyz;
  int jetfreetrksvtxntrk, jetfreetrkshasvtx;
  int jetbestsvtxtrks1add, jetbestsvtxtrks3add, jetbestsvtxtrks5add;

  int jetnum, evnum, evjetcount, evnfreetrks, evpvtxntrk;
};


NewTertVertex::NewTertVertex(const edm::ParameterSet& iConfig)
{
  jetTrackLabel      = iConfig.getParameter<InputTag>("jetTracks");
  primVtxLabel       = iConfig.getParameter<InputTag>("primaryVertex");
  outfile            = iConfig.getParameter<string>("outfilename");
  jfi = JetFlavourIdentifier(iConfig.getParameter<edm::ParameterSet>("jetIdParameters"));
  jetcorrector = new SimpleMCJetCorrector(); 
  jetcorrector->init("CMSSW_152_iterativeCone5.txt");   

  vtxRecoPSet        = iConfig.getParameter<edm::ParameterSet>("vertexReco");
  trackadder = AddTertiaryTracks();
  trackadder.configure(vtxRecoPSet);
 
  thefile = new TFile(outfile.c_str(),"recreate");
  thefile->cd();
  thetree= new TTree("ntp1","ntp1");
  
  //  thetree->Branch("jetnum",&jetnum,"jetnum/I");
  thetree->Branch("evjetcount",&evjetcount,"evjetcount/I");
  thetree->Branch("evnum",&evnum,"evnum/I");
  thetree->Branch("evpvtxntrk",&evpvtxntrk,"evpvtxntrk/I");
  thetree->Branch("evnfreetrks",&evnfreetrks,"evnfreetrks/I");
  thetree->Branch("jetntracks",&jetntracks,"jetntracks/I");
  thetree->Branch("jetnsecvtx",&jetnsecvtx,"jetnsecvtx/I");

  thetree->Branch("jetpx",&jetpx,"jetpx/D");
  thetree->Branch("jetpy",&jetpy,"jetpy/D"); 
  thetree->Branch("jetpz",&jetpz,"jetpz/D"); 
  thetree->Branch("jetpt",&jetpt,"jetpt/D");
  thetree->Branch("jete",&jete,"jete/D");  
  thetree->Branch("jeteta",&jeteta,"jeteta/D"); 
  thetree->Branch("jetphi",&jetphi,"jetphi/D");  

  thetree->Branch("jetpartonpx",&jetpartonpx,"jetpartonpx/D"); 
  thetree->Branch("jetpartonpy",&jetpartonpy,"jetpartonpy/D");  
  thetree->Branch("jetpartonpz",&jetpartonpz,"jetpartonpz/D");  
  thetree->Branch("jetpartonpt",&jetpartonpt,"jetpartonpt/D"); 
  thetree->Branch("jetpartoneta",&jetpartoneta,"jetpartoneta/D");  
  thetree->Branch("jetpartonphi",&jetpartonphi,"jetpartonphi/D");   

  thetree->Branch("jetflavour",&jetflavour,"jetflavour/I");
  thetree->Branch("jettrkbtag",&jettrkbtag,"jettrkbtag/D"); 
  thetree->Branch("jetsmbtag",&jetsmbtag,"jetsmbtag/D"); 
  thetree->Branch("jetjpbtag",&jetjpbtag,"jetjpbtag/D");
  thetree->Branch("jetcorr",&jetcorr,"jetcorr/D"); 
  thetree->Branch("jetprimvtxntrk",&jetprimvtxntrk,"jetprimvtxntrk/I");

  thetree->Branch("jetbestsvtxdxy",&jetbestsvtxdxy,"jetbestsvtxdxy/D");
  thetree->Branch("jetbestsvtxx",&jetbestsvtxx,"jetbestsvtxx/D"); 
  thetree->Branch("jetbestsvtxy",&jetbestsvtxy,"jetbestsvtxy/D"); 
  thetree->Branch("jetbestsvtxz",&jetbestsvtxz,"jetbestsvtxz/D");  
  thetree->Branch("jetbestsvtxchi2",&jetbestsvtxchi2,"jetbestsvtxchi2/D");  
  thetree->Branch("jetbestsvtxntrk",&jetbestsvtxntrk,"jetbestsvtxntrk/I");   
  thetree->Branch("jetbestsvtxmass",&jetbestsvtxmass,"jetbestsvtxmass/D");   
  thetree->Branch("jetbestsvtxtrks1add",&jetbestsvtxtrks1add,"jetbestsvtxtrks1add/I");
  thetree->Branch("jetbestsvtxtrks3add",&jetbestsvtxtrks3add,"jetbestsvtxtrks3add/I"); 
  thetree->Branch("jetbestsvtxtrks5add",&jetbestsvtxtrks5add,"jetbestsvtxtrks5add/I"); 

  thetree->Branch("jetnfreetrks",&jetnfreetrks,"jetnfreetrks/I");

  thetree->Branch("primvtxposx",&primvtxposx,"primvtxposx/D");
  thetree->Branch("primvtxposy",&primvtxposy,"primvtxposy/D"); 
  thetree->Branch("primvtxposz",&primvtxposz,"primvtxposz/D"); 

  thetree->Branch("jetfreetrksvtxmass",&jetfreetrksvtxmass,"jetfreetrksvtxmass/D");
  thetree->Branch("jetfreetrksvtxposx",&jetfreetrksvtxposx,"jetfreetrksvtxposx/D");
  thetree->Branch("jetfreetrksvtxposy",&jetfreetrksvtxposy,"jetfreetrksvtxposy/D"); 
  thetree->Branch("jetfreetrksvtxposz",&jetfreetrksvtxposz,"jetfreetrksvtxposz/D"); 
  thetree->Branch("jetfreetrksvtxntrk",&jetfreetrksvtxntrk,"jetfreetrksvtxntrk/I"); 
  thetree->Branch("jetfreetrksvtxchi2",&jetfreetrksvtxchi2,"jetfreetrksvtxchi2/D");
  thetree->Branch("jetfreetrksvtxrxy",&jetfreetrksvtxrxy,"jetfreetrksvtxrxy/D");
  thetree->Branch("jetfreetrksvtxrxyz",&jetfreetrksvtxrxyz,"jetfreetrksvtxrxyz/D"); 
  thetree->Branch("jetfreetrkshasvtx",&jetfreetrkshasvtx,"jetfreetrkshasvtx/I"); 


  evnum = 0;
}


NewTertVertex::~NewTertVertex() 
{  
  thefile->Write();
  thefile->Close();
}

void NewTertVertex::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  cout <<"\n============================================================================\n";
  //  cout <<"\nNew Event\n";
  //  cout <<"\n============================================================================\n";

  jfi.readEvent(iEvent);
  vector<MCParton> partonList = jfi.getListOfPartons();
		  
  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel("source", evt);
  //  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(evt->GetEvent())); 

  // Get the Jet collection from the event
  edm::Handle<reco::CaloJetCollection> pJets;
  iEvent.getByLabel("iterativeCone5CaloJets",pJets);
  //  const reco::CaloJetCollection* jets = pJets.product();
  //  reco::CaloJetCollection::const_iterator jet;

  // Get jet probability tags
  edm::Handle<reco::JetTagCollection> jpTagHandle;
  iEvent.getByLabel("jetBProbabilityBJetTags",jpTagHandle);
  const reco::JetTagCollection & jpTags = *(jpTagHandle.product());

  // Get track counting tags
  edm::Handle<reco::JetTagCollection> trkTagHandle; 
  iEvent.getByLabel("trackCountingHighEffBJetTags", trkTagHandle); 
  const reco::JetTagCollection & trkTags = *(trkTagHandle.product()); 

  // Get Soft muon tags
  edm::Handle<reco::JetTagCollection> smbTagHandle;  
  iEvent.getByLabel("softMuonBJetTags", smbTagHandle);  
  const reco::JetTagCollection & smbTags = *(smbTagHandle.product());  

  // Get secondary vertex infos
  edm::Handle<reco::SecondaryVertexTagInfoCollection> secvtxTagHandle;
  iEvent.getByLabel("secondaryVertexTagInfos", secvtxTagHandle);
  const reco::SecondaryVertexTagInfoCollection & secvtxTags = *(secvtxTagHandle.product());

  // Make Transient tracks from RECO tracks
  edm::Handle<reco::TrackCollection> tks;
  //  iEvent.getByLabel("pixelTracks", tks);
  iEvent.getByLabel("generalTracks", tks); 
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  //do the conversion:
  vector<TransientTrack> mytracks = (*theB).build(tks);
  std::vector<reco::TransientTrack> evfreetrks;
  std::vector<reco::TransientTrack> jetfreetrks;
  std::vector<reco::TransientTrack>::const_iterator evfreetrk;
  std::vector<reco::TransientTrack>::const_iterator jetfreetrk;
  std::vector<reco::SecondaryVertexTagInfo> jjhsecvtxs;

  // Jet<->tracks
  ProductID jets_id;
  Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
  iEvent.getByLabel("ic5JetTracksAssociatorAtVertex", jetTracksAssociation);
  reco::JetTracksAssociationCollection::const_iterator assocjet;

  // Track quality
  edm::Handle<TrackIPTagInfoCollection> trackIPTagInfos; 
  iEvent.getByLabel("impactParameterTagInfos", trackIPTagInfos); 

  // Primary vertex
  Handle<reco::VertexCollection> pvs;
  iEvent.getByLabel(primVtxLabel,pvs);
  reco::Vertex pv = pvs->front();
  primvtxposx = pv.x();
  primvtxposy = pv.y();
  primvtxposz = pv.z();

  // Vertexer
  TertiaryTracksVertexFinder fitter; 
  ConfigurableVertexReconstructor finder(vtxRecoPSet);

  edm::ESHandle<TransientTrackBuilder> builder; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder ); 
  vector < reco::TransientTrack > trks;
  vector < reco::TransientTrack >::const_iterator trk;  

  jetnum = jetTracksAssociation->size();
  evjetcount = 0;
  evpvtxntrk = 0;
  jetflavour = -1;
  evnfreetrks = 0;

  VertexMass theVertexMass;
  Vertex::trackRef_iterator pviter;  
  TrackRefVector::iterator sviter;

  math::XYZTLorentzVector jetparton;

  evfreetrks.clear();  

  // Get all tracks in the event not associated with a primary or secondary vertex
  for(evfreetrk = mytracks.begin(); evfreetrk != mytracks.end(); ++evfreetrk)  
    {  
      bool freetrack = true;
      //      cout << "Track with pT = " << evfreetrk->track().pt() << endl;
      for(pviter = pv.tracks_begin(); pviter != pv.tracks_end();++pviter)
	{
	  evpvtxntrk++;
	  if((*pviter)->px() == evfreetrk->track().px())
	    {
	      freetrack = false;
	    }
	}
      for(int svtxiter = 0; svtxiter < secvtxTags.size(); ++svtxiter)
      	{
	  TrackRefVector thetmptrks = (secvtxTags[svtxiter]).vertexTracks();

	  for(sviter = thetmptrks.begin(); sviter != thetmptrks.end(); ++ sviter)
	    {
	      if((*sviter)->px() == evfreetrk->track().px()) 
		{
		  freetrack = false;
		}
	    }
	}

      if(freetrack == true)
	{
	  evfreetrks.push_back( *evfreetrk );   
	  evnfreetrks++;
	}
    }  

  TrackIPTagInfo::SortCriteria    sortCriterium;
  sortCriterium = TrackSorting::getCriterium("sip3dSig");

  // Now look at track quality


  // Now look at jets
  //  for(assocjet = jetTracksAssociation->begin();assocjet != jetTracksAssociation->end(); ++assocjet)
  //    {
  for(TrackIPTagInfoCollection::const_iterator assocjets = 
        trackIPTagInfos->begin(); assocjets != trackIPTagInfos->end(); 
      ++assocjets)  
    {

      std::vector<std::size_t> indices = 
        assocjets->sortedIndexes(sortCriterium); 
       
      const TrackRefVector &trackRefs = 
        assocjets->sortedTracks(indices); 
 
      edm::RefToBase<Jet> assocjet = assocjets->jet(); 


      jetnfreetrks = 0;
      jetprimvtxntrk = 0;
      jetfreetrkshasvtx = 0;
      jetfreetrksvtxmass = 0.0;
      jetfreetrksvtxposx = 0.0;
      jetfreetrksvtxposy = 0.0;
      jetfreetrksvtxposz = 0.0;
      jetfreetrksvtxntrk = 0;
      jetfreetrksvtxchi2 = 0.0;
      jetfreetrksvtxrxy = 0.0;
      jetfreetrksvtxrxyz = 0.0;
      jetbestsvtxtrks1add = 0; 
      jetbestsvtxtrks3add = 0;
      jetbestsvtxtrks5add = 0;

      trks.clear();
      jetfreetrks.clear();
      bool jetfreetrack = true;

      // Get the tracks associated with this jet
      for(unsigned int i = 0; i < indices.size(); i++)  
        { 
	  trks.push_back ( builder->build ( trackRefs[i] ) ); 
	  //      for ( edm::RefVector < reco::TrackCollection >::const_iterator i=assocjet->second.begin();
	  //	    i!=assocjet->second.end() ; ++i )
	  //	  trks.push_back ( builder->build ( *i ) );
	}

      // See if any of the tracks are not attached to a primary or secondary vertex
      for(trk = trks.begin(); trk != trks.end(); ++trk) 
	{
	  jetfreetrack = true;

	  for(int svtxiter = 0; svtxiter < secvtxTags.size(); ++svtxiter) 
	    { 
	      TrackRefVector thetmptrks = (secvtxTags[svtxiter]).vertexTracks(); 
 
	      for(sviter = thetmptrks.begin(); sviter != thetmptrks.end(); ++ sviter) 
		{ 
		  if((*sviter)->pt() == trk->track().pt())
		    {
		      jetfreetrack = false;
		    }
		} 
	    } 
	  for(pviter = pv.tracks_begin(); pviter != pv.tracks_end();++pviter) 
	    { 
	      if((*pviter)->pt() == trk->track().pt()) 
		{
		  jetprimvtxntrk++;
		  jetfreetrack = false;
		}
	    } 

	  // If they're unattached, add to the list of free tracks
	  if(jetfreetrack == true)
	    {
	      jetfreetrks.push_back( * trk ); 
	      jetnfreetrks++;
	    }
	} 

      // Get jet and b-tagging information
      jetpx = assocjet->px();
      jetpy = assocjet->py();
      jetpz = assocjet->pz();
      jetpt = assocjet->pt();
      jete = assocjet->energy();
      jeteta = assocjet->eta();
      jetphi = assocjet->phi();

      jetcorr = jetcorrector->correctionXYZT(jetpx,jetpy,jetpz,jete); 

      for (int i = 0; i != smbTags.size(); ++i) 
	{
	  if((smbTags[i].first->eta() == jeteta) && (smbTags[i].first->phi() == jetphi))
	    jetsmbtag = smbTags[i].second;
	}
      
      for (int i = 0; i != trkTags.size(); ++i)  
	{ 
	  if((trkTags[i].first->eta() == jeteta) && (trkTags[i].first->phi() == jetphi)) 
	    jettrkbtag = trkTags[i].second; 
	} 
      
      for (int i = 0; i != jpTags.size(); ++i)
	{
	  if((jpTags[i].first->eta() == jeteta) && (jpTags[i].first->phi() == jetphi))
	    jetjpbtag = jpTags[i].second;
	}
      
      // Get jet truth
      BTagMCTools::JetFlavour jetFlavour = jfi.identifyBasedOnPartons(* (assocjet));
      jetflavour = jetFlavour.flavour();
      jetparton = jetFlavour.underlyingParton4Vec();
      jetpartonpx = jetparton.x();
      jetpartonpy = jetparton.y();
      jetpartonpz = jetparton.z();
      jetpartonpt = jetparton.t();
      jetpartoneta = jetparton.eta();
      jetpartonphi = jetparton.phi();
      
      // Now get secondary vertex information
      jetntracks = trks.size(); 
      jetnsecvtx = 0;

      for (int svtxiter = 0; svtxiter < secvtxTags.size(); ++svtxiter) 
	{ 
	  double thejpt = secvtxTags[svtxiter].jet()->pt();  
	  if(thejpt == jetpt)
	    { 
	      jjhsecvtxs.push_back( (secvtxTags[svtxiter]) );

	      double esum=0., pxsum=0., pysum=0., pzsum=0.;

	      if(secvtxTags[svtxiter].nVertices() != 0) 
		{ 
		  jetbestsvtxdxy = secvtxTags[svtxiter].flightDistance(0, true).value();
		  jetbestsvtxntrk = secvtxTags[svtxiter].nVertexTracks(0); 
		  jetbestsvtxchi2 = (secvtxTags[svtxiter]).secondaryVertex(0).normalizedChi2(); 
		  jetbestsvtxx = (secvtxTags[svtxiter]).secondaryVertex(0).x(); 
		  jetbestsvtxy = (secvtxTags[svtxiter]).secondaryVertex(0).y(); 
		  jetbestsvtxz = (secvtxTags[svtxiter]).secondaryVertex(0).z(); 

		  //		  jetbestsvtxmass = theVertexMass((secvtxTags[svtxiter]).secondaryVertex(0));
		  // Get the vertex mass
		  TrackRefVector thetmpstrks = (secvtxTags[svtxiter]).vertexTracks(0); 
		  for(sviter = thetmpstrks.begin(); sviter != thetmpstrks.end(); ++ sviter)
		    {
		      double px = (*sviter)->px();
		      double py = (*sviter)->py(); 
		      double pz = (*sviter)->pz(); 
		      pxsum+=px;
                      pysum+=py;
                      pzsum+=pz;
		      esum += sqrt(px*px + py*py + pz*pz + 0.13957*0.13957);
		    }

		  jetbestsvtxmass = sqrt(esum*esum - (pxsum*pxsum + pysum*pysum + pzsum*pzsum));
		  jetnsecvtx++;
		} 
	      else
		{
		  jetbestsvtxmass = -1;
                  jetbestsvtxdxy = -1;
                  jetbestsvtxntrk = -1;
                  jetbestsvtxchi2 = -1;
                  jetbestsvtxx = 0;
                  jetbestsvtxy = 0;
                  jetbestsvtxz = 0;
		  jetnsecvtx = 0;
		}
	    } 
	} 

      // Now feed the list of free tracks to the trackadder, 
      // and look for tertiary vertexes
      trackadder.attachtracks(jetfreetrks,pv,jjhsecvtxs);

      jetbestsvtxtrks1add = trackadder.getsvtxtrks1();
      jetbestsvtxtrks3add = trackadder.getsvtxtrks3();
      jetbestsvtxtrks5add = trackadder.getsvtxtrks5();

      jetfreetrkshasvtx = trackadder.hastertvtx();

      if(jetfreetrkshasvtx == 1)
	{
	  jetfreetrksvtxmass = trackadder.getmass();
	  jetfreetrksvtxposx = trackadder.vertex().x();
	  jetfreetrksvtxposy = trackadder.vertex().y();
	  jetfreetrksvtxposz = trackadder.vertex().z();
	  jetfreetrksvtxntrk = trackadder.getnterttracks();
	  jetfreetrksvtxchi2 = trackadder.getchi2();
	  jetfreetrksvtxrxy = trackadder.getrxy();
	  jetfreetrksvtxrxyz = trackadder.getrxyz();
	}
      else
	{
          jetfreetrksvtxmass = 0;
          jetfreetrksvtxposx = 0;
          jetfreetrksvtxposy = 0;
          jetfreetrksvtxposz = 0;
          jetfreetrksvtxntrk = 0;
          jetfreetrksvtxchi2 = 0;
          jetfreetrksvtxrxy = 0;
          jetfreetrksvtxrxyz = 0;
	}
      
      thetree->Fill();

      evjetcount++;

    }

  evnum++;
}
  

//define this as a plug-in
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(NewTertVertex);
