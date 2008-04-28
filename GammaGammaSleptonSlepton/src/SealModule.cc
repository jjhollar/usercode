#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "UserCode/GammaGammaSleptonSlepton/interface/GammaGammaSleptonSlepton.h"
#include "UserCode/GammaGammaSleptonSlepton/interface/GenGammaGammaSleptonSlepton.h"

DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(GammaGammaSleptonSlepton);
DEFINE_FWK_MODULE(GenGammaGammaSleptonSlepton);
