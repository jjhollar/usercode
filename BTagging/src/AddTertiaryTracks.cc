#include "RecoBTag/MCTools/interface/AddTertiaryTracks.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "RecoBTag/MCTools/interface/MCParticleInfo.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <Math/GenVector/VectorUtil.h>
#include "DataFormats/Math/interface/Vector.h"
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
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h" 

using namespace edm;
using namespace std;
using namespace HepMC;
using namespace reco;
using namespace math;

AddTertiaryTracks::AddTertiaryTracks()
{
  vtxmass_ = 0.0;
  posx_ = 0.0;
  posy_ = 0.0;
  posz_ = 0.0;
  vtxvalid_ = 0;
  vtxntrk_ = 0;
  vtxchi2_ = 0.0;
  vtxrxy_ = 0.0;
  vtxrxyz_ = 0.0;
  nsvtxtrks1_ = 0; nsvtxtrks3_ = 0; nsvtxtrks5_ = 0;
}

void 
AddTertiaryTracks::configure(const edm::ParameterSet& iConfigPset)
{
  vtxConfigPSet_ = iConfigPset;
}

void 
AddTertiaryTracks::attachtracks(const std::vector<reco::TransientTrack>  freetracks, reco::Vertex & pv, std::vector<reco::SecondaryVertexTagInfo> jjhsecvtxs)
{
  double esum=0., pxsum=0., pysum=0., pzsum=0., mass=0., chitwo=0.,xpos=0.,ypos=0.,zpos=0.; 
  double rxy=0.,rxyz=0.; 
  double doca = 0., docaxy = 0., docaz = 0.; 
  int hasvtx = 0; int valid = 0;
  int ntracks = 0;
  int naddsvtxtrkspt1 = 0; 
  int naddsvtxtrkspt3 = 0;
  int naddsvtxtrkspt5 = 0; 

  double thePionMass = 0.13957; 

  // Look for tracks close to secondary vertex
  for (std::vector<reco::TransientTrack>::const_iterator iter = freetracks.begin();
       iter != freetracks.end(); ++iter)
    {
      for(int svtxiter = 0; svtxiter < jjhsecvtxs.size(); ++svtxiter)
	{
	  if(jjhsecvtxs[svtxiter].nVertices() != 0)  
	    {  
	      Point svtxpos = Point((jjhsecvtxs[svtxiter]).secondaryVertex(0).x(),(jjhsecvtxs[svtxiter]).secondaryVertex(0).y(),(jjhsecvtxs[svtxiter]).secondaryVertex(0).z());
	      double docaxy = (iter->track()).dxy(svtxpos);
	      double docaz = (iter->track()).dz(svtxpos);
	      double doca = sqrt(docaxy*docaxy + docaz*docaz);
	      if(doca < 0.1)
		naddsvtxtrkspt1++; 
	      if(doca < 0.3)
		naddsvtxtrkspt3++;
	      if(doca < 0.5)
		naddsvtxtrkspt5++; 
	      //	      cout << "Doca to secondary vertex at (" << svtxpos.x() << ", " << svtxpos.y() << ", " << svtxpos.z() << ") is " << doca << endl;
	      //	      cout << "\tDocaxy = " << docaxy << ", Docaz = " << docaz << endl;
	    }
	}
    }


  ConfigurableVertexReconstructor vertexer(vtxConfigPSet_); 
  vector<TransientVertex> tertvertexes = vertexer.vertices(freetracks); 
  vector<TransientTrack> otracks;

  if(tertvertexes.size() > 0)
    {
      otracks = (tertvertexes[0]).originalTracks();  

      hasvtx = 1;

      for (std::vector<reco::TransientTrack>::const_iterator it=otracks.begin();
      	   it!=otracks.end();++it)
	{
	  reco::TransientTrack track = *it; 
	  
	  double px = track.impactPointState().globalMomentum().x(); 
	  double py = track.impactPointState().globalMomentum().y(); 
	  double pz = track.impactPointState().globalMomentum().z(); 
	  
	  pxsum += px; 
	  pysum += py; 
	  pzsum += pz; 
	  esum += sqrt(px*px + py*py + pz*pz + thePionMass*thePionMass); 

	  ntracks++;
	} 
      
      mass = sqrt(esum*esum - (pxsum*pxsum + pysum*pysum + pzsum*pzsum));
      chitwo = (tertvertexes[0]).normalisedChiSquared();  
      valid = (tertvertexes[0]).isValid();
      xpos = (tertvertexes[0]).position().x(); 
      ypos = (tertvertexes[0]).position().y(); 
      zpos = (tertvertexes[0]).position().z(); 
      rxy = sqrt(xpos*xpos + ypos*ypos);
      rxyz = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

    }

  //  vtxvalid_ = hasvtx;
  vtxmass_ = mass;
  vtxchi2_ = chitwo;
  vtxntrk_ = ntracks; 
  vtxrxy_ = rxy;
  vtxrxyz_ = rxyz;
  vtxvalid_ = valid;
  posx_ = xpos;
  posy_ = ypos;
  posz_ = zpos;
  nsvtxtrks1_ = naddsvtxtrkspt1;
  nsvtxtrks3_ = naddsvtxtrkspt3;
  nsvtxtrks5_ = naddsvtxtrkspt5;
}



