#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayloadFromTable.h"

#include "CondFormats/PhysicsToolsObjects/interface/PerformanceWorkingPoint.h"

#include <TFile.h>  
#include <TH1D.h>  
#include <TH2F.h>

class PhysicsPerformanceDBWriterFromTPHist : public edm::EDAnalyzer
{
public:
  PhysicsPerformanceDBWriterFromTPHist(const edm::ParameterSet&);
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) {}
  virtual void endJob() {}
  ~PhysicsPerformanceDBWriterFromTPHist() {}

private:
  std::string inputHistFile;
  std::string rec1,rec2;
  boost::uint64_t iovBegin,iovEnd;  

  TFile *f1;
  TH2F *heff;
};

PhysicsPerformanceDBWriterFromTPHist::PhysicsPerformanceDBWriterFromTPHist
  (const edm::ParameterSet& p)
{
  inputHistFile = p.getUntrackedParameter<std::string>("inputHistFile");
  rec1 = p.getUntrackedParameter<std::string>("RecordPayload");
  rec2 = p.getUntrackedParameter<std::string>("RecordWP");
  iovBegin = p.getParameter<boost::uint64_t>("IOVBegin"); 
  iovEnd = p.getParameter<boost::uint64_t>("IOVEnd"); 
}

void PhysicsPerformanceDBWriterFromTPHist::beginJob(const edm::EventSetup&)
{
  //
  // read object from file
  //
  f1 = TFile::Open("tphist.root");

  //
  // File Format is
  // - tagger name
  // - cut
  // - concrete class name
  // - how many results and how many binning
  // - their values
  // - vector<float>
  //

  std::ifstream in;
  in.open(inputHistFile.c_str());
  std::string tagger;
  float cut;
 
  std::string concreteType;
  std::string comment;
  std::vector<float> pl;
  int stride;

  in >> tagger;
  std::cout << "WP Tagger is "<<tagger<<std::endl;
  
  in >> cut;
  std::cout << "WP Cut is "<<cut<<std::endl;

  in >> concreteType;
  std::cout << "concrete Type is "<<concreteType<<std::endl;

  //  return ;

  // read # of results

  int nres, nbin;
  in >> nres;
  in >> nbin;
  std::cout <<" Results: " << nres<<" Binning variables: "<<nbin<<std::endl;

  stride = nres+nbin*2;
  
  int number=0;
  
  std::vector<PerformanceResult::ResultType> res;
  std::vector<BinningVariables::BinningVariablesType> bin;

  while (number<nres && !in.eof()) {
    int tmp;
    in>> tmp;
    res.push_back((PerformanceResult::ResultType)(tmp));
    number++;
  }
  if (number != nres){
    std::cout <<" Table not well formed"<<std::endl;
  }
  number=0;
  while (number<nbin && !in.eof()) {
    int tmp;
    in>> tmp;
    bin.push_back((BinningVariables::BinningVariablesType)(tmp));
    number++;
  }
  if (number != nbin){
    std::cout <<" Table not well formed"<<std::endl;
  }

  number=0;
  while (!in.eof()){
    float temp;
    in >> temp;
    std::cout <<" Intersing "<<temp<< " in position "<<number<<std::endl;
    number++;
    pl.push_back(temp);
  }

  //
  // CHECKS
  //
  if (stride != nbin*2+nres){
    std::cout <<" Table not well formed"<<std::endl;
  }
  if ((number % stride) != 0){
    std::cout <<" Table not well formed"<<std::endl;
  }

  in.close();


  /*  for (int k=0;k<(number/stride); k++){
    for (int j=0; j<stride; j++){
      std::cout << "Pos["<<k<<","<<j<<"] = "<<pl[k*stride+j]<<std::endl;
    }
  }
  */

  //
  // now create pl etc etc
  //

  PerformanceWorkingPoint * wp = new PerformanceWorkingPoint(cut, tagger);

  PerformancePayloadFromTable * btagpl = 0;

  if (concreteType == "PerformancePayloadFromTable"){
    btagpl = new PerformancePayloadFromTable(res, bin, stride, pl);
  }else{
    std::cout <<" Non existing request: " <<concreteType<<std::endl;
  }
  
  std::cout <<" Created the "<<concreteType <<" object"<<std::endl;
  
  edm::Service<cond::service::PoolDBOutputService> s;
  if (s.isAvailable())
    {
      if (s->isNewTagRequest(rec1))
	{
	  s->createNewIOV<PerformancePayload>(btagpl,
					      //						  s->beginOfTime(),
					      //						  s->endOfTime(),
					      iovBegin,
					      iovEnd,
						  rec1);
	}
      else
	{
	  
	  s->appendSinceTime<PerformancePayload>(btagpl,
						 // JUST A STUPID PATCH
						 //						     111,
						 iovBegin,
						     rec1);
	}
    }

  // write also the WP
  
  if (s.isAvailable())
    {
      if (s->isNewTagRequest(rec2))
	{
	  s->createNewIOV<PerformanceWorkingPoint>(wp,
						   //					    s->beginOfTime(),
						   //					    s->endOfTime(),
						   iovBegin,
						   iovEnd,
					    rec2);
	}
      else
	{
	  
	  s->appendSinceTime<PerformanceWorkingPoint>(wp,
						      /// JUST A STUPID PATCH
						      //					       111,
						      iovBegin,
					       rec2);
	}
    }
  
  
  
}

DEFINE_FWK_MODULE(PhysicsPerformanceDBWriterFromTPHist);
