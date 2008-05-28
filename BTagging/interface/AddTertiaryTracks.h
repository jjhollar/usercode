#ifndef AddTertiaryTracks_H
#define AddTertiaryTracks_H

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoBTag/MCTools/interface/MCParton.h"
#include "RecoBTag/MCTools/interface/JetFlavour.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/VertexReco/interface/VertexFwd.h" 
#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h" 
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/Vector3D.h" 
#include "DataFormats/Math/interface/LorentzVector.h" 
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h" 
#include "DataFormats/GeometryVector/interface/GlobalVector.h" 
#include "DataFormats/GeometryVector/interface/VectorUtil.h" 

class AddTertiaryTracks {
public:

  typedef math::XYZPoint Point;

  AddTertiaryTracks();

  void attachtracks(const std::vector<reco::TransientTrack>  freetracks, reco::Vertex & pv, std::vector<reco::SecondaryVertexTagInfo>  jjhsecvtxs);

  double getmass() const {return vtxmass_;}

  int hastertvtx() const {return vtxvalid_;}

  int getnterttracks() const {return vtxntrk_;}

  double getchi2() const {return vtxchi2_;}

  void configure(const edm::ParameterSet& iConfigPset);

  GlobalPoint vertex() const {return GlobalPoint(posx_,posy_,posz_);}

  double getrxyz() const {return vtxrxyz_;}

  double getrxy() const {return vtxrxy_;}

  int getsvtxtrks1() const {return nsvtxtrks1_;}
  int getsvtxtrks3() const {return nsvtxtrks3_;} 
  int getsvtxtrks5() const {return nsvtxtrks5_;} 

private:
  double vtxmass_;
  double posx_;
  double posy_;
  double posz_;
  double vtxrxy_;
  double vtxrxyz_;
  int vtxntrk_;
  int vtxvalid_;
  double vtxchi2_;
  int nsvtxtrks1_;
  int nsvtxtrks3_;
  int nsvtxtrks5_;

  edm::ParameterSet  vtxConfigPSet_;
};

#endif
