////////////////////////////////////////////////////////////////////////
// Class: SystToolsEventWeight
// Module Type: producer
// File: SystToolsEventWeight_module.cc
//
// Generated at Mon Feb 21 09:36:11 2015 by Jaesung Kim (jae.sung.kim.3426@gmail.com)
//
////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Assns.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h"

#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interpreters/ParamHeaderHelper.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "nusystematics/utility/response_helper.hh"

#include "larnusystematics/utility/fhiclsimple_tools.h"
#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "Framework/EventGen/EventRecord.h"

#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Run.h"


namespace sbn::evwgh {

class SystToolsEventWeight : public art::EDProducer {
public:
  explicit SystToolsEventWeight(fhicl::ParameterSet const& p);

  SystToolsEventWeight(SystToolsEventWeight const &) = delete;
  SystToolsEventWeight(SystToolsEventWeight &&) = delete;
  SystToolsEventWeight& operator = (SystToolsEventWeight const&) = delete;
  SystToolsEventWeight& operator = (SystToolsEventWeight&&) = delete;

private:
  void produce(art::Event& e) override;
  void beginRun(art::Run& run) override;

private:
  std::string fGenieModuleLabel;
  bool fAllowMissingTruth;
  bool fDebugMode;

  nusyst::response_helper fParamHeaderHelper;
  //systtools::provider_list_t fSystProviders;
  //systtools::param_header_map_t fParamHeaderMap;

  /*
  systtools::provider_list_t fSystProviders; -> res_helper.syst_providers
  systtools::param_header_map_t fParamHeaderMap; -> resp_helper.fHeaders
  systtools::ParamHeaderHelper fParamHeaderHelper; -> resp_helper itself
  */

};


SystToolsEventWeight::SystToolsEventWeight(fhicl::ParameterSet const& p)
  : EDProducer{p}
{

  fGenieModuleLabel = p.get<std::string>("GeneratorModuleLabel", "generator");
  fAllowMissingTruth = p.get<bool>("AllowMissingTruth");
  fDebugMode = p.get<bool>("DebugMode", false);

  fhicl::ParameterSet f_syst_provider_config = p.get<fhicl::ParameterSet>("generated_systematic_provider_configuration");
  fhiclsimple::ParameterSet fs_syst_provider_config = nusyst::FhiclToFhiclSimple(f_syst_provider_config);

  MF_LOG_DEBUG("SystToolsEventWeight") << "Configuring ISystProvider";

/*
  std::function< std::unique_ptr<systtools::ISystProviderTool>(fhiclsimple::ParameterSet const &fs) > InstanceBuilder =
  [](fhiclsimple::ParameterSet const &fs) -> std::unique_ptr<systtools::ISystProviderTool> {
    std::cout << "[JSKIMDEBUG][InstanceBuilder] called" << std::endl;
    return art::make_tool<systtools::ISystProviderTool>( nusyst::FhiclSimpleToFhicl(fs) );
  };

  // note: typedef std::vector<std::unique_ptr<ISystProviderTool>> provider_list_t;
  fSystProviders = systtools::ConfigureISystProvidersFromParameterHeaders<systtools::ISystProviderTool>(
    fs_syst_provider_config,
    InstanceBuilder,
    "syst_providers"
  );

  std::cout << "[JSKIMDEBUG][SystToolsEventWeight::SystToolsEventWeight] HERE 2" << std::endl;

  fParamHeaderMap = systtools::BuildParameterHeaders(fSystProviders);
  fParamHeaderHelper = systtools::ParamHeaderHelper(fParamHeaderMap);
*/

  
  fParamHeaderHelper.LoadProvidersAndHeaders(fs_syst_provider_config);
  //fSystProviders = fParamHeaderHelper.syst_providers;
  //fParamHeaderMap = fParamHeaderHelper.fHeaders;

/*
  systtools::provider_list_t fSystProviders; -> fParamHeaderHelper.syst_providers
  systtools::param_header_map_t fParamHeaderMap; -> fParamHeaderHelper.fHeaders
  systtools::ParamHeaderHelper fParamHeaderHelper; ->fParamHeaderHelper itself
*/

  MF_LOG_DEBUG("SystToolsEventWeight") << "ISystProvider is configuered";

  produces<std::vector<sbn::evwgh::EventWeightMap> >();
  produces<art::Assns<simb::MCTruth, sbn::evwgh::EventWeightMap> >();
  produces<std::vector<sbn::evwgh::EventWeightParameterSet>, art::InRun>();

}


void SystToolsEventWeight::produce(art::Event& e) {

  auto mcwghvec = std::make_unique<std::vector<sbn::evwgh::EventWeightMap> >();
  auto wghassns = std::make_unique<art::Assns<simb::MCTruth, sbn::evwgh::EventWeightMap> >();

  art::PtrMaker<sbn::evwgh::EventWeightMap> makeWeightPtr(e);

  std::vector<art::Ptr<simb::MCTruth> > mclist;
  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle= e.getHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);
  art::Handle<std::vector<simb::GTruth>> gTruthHandle = e.getHandle<std::vector<simb::GTruth>>(fGenieModuleLabel);

  size_t NEventUnits = mcTruthHandle->size();
  if (mcTruthHandle->size() != gTruthHandle->size()) {
    NEventUnits = std::min(mcTruthHandle->size(), gTruthHandle->size());
  }
  std::vector<std::unique_ptr<genie::EventRecord>> gheps; // TODO If we alread have genie::EventRecord, use that instead of making it through RetrieveGHEP
  for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
    gheps.emplace_back(evgb::RetrieveGHEP(mcTruthHandle->at(eu_it),
                                          gTruthHandle->at(eu_it)));
  }

  // Prooceed even with missing handle if we want to require the MCTruth to be
  // found, so that an exception will be thrown explaining the problem.
  if(mcTruthHandle.isValid() || !fAllowMissingTruth){

    art::fill_ptr_vector(mclist, mcTruthHandle);

    if(fDebugMode){
      std::cout << "[SystToolsEventWeight::produce] mclist.size() = " << mclist.size() << "\n"
                << "[SystToolsEventWeight::produce] fParamHeaderHelper.GetSystProvider().size() = " << fParamHeaderHelper.GetSystProvider().size() << std::endl;
    }

    for( auto &sp : fParamHeaderHelper.GetSystProvider() ){

      // TODO this part must be improved;
      // sp is "std::unique_ptr<ISystProviderTool", and this parent class does not have functions like GetEventResponse()
      std::unique_ptr<systtools::EventResponse> syst_resp = sp->GetEventResponse(gheps);
      //std::unique_ptr<systtools::EventResponse> syst_resp = ( dynamic_cast<nusyst::IGENIESystProvider_tool*>(sp.get()) )->GetEventResponse(gheps);

      if( !syst_resp ){
        throw cet::exception{ "SystToolsEventWeight" }
          << "Got nullptr systtools::EventResponse from provider "
          << sp->GetFullyQualifiedName() << "\n";
      }

      // syst_resp->size() is (Number of MCTruth)
      // Each index corresponds to each of MCTruth
      int nMCTruthIndex = 0;
      if(fDebugMode) std::cout << "[SystToolsEventWeight::produce]   syst_resp.size() (= Number of MCTruth) of this SystProvider = " << syst_resp->size() << std::endl;

      // Looping over syst_resp is identical to looping over MCTruth
      for(systtools::event_unit_response_t const& resp: *syst_resp) {
        // resp.size() corresponds to number of knobs we altered;
        // e.g., MaCCQE, MaCCRES, MvCCRE -> resp.size() = 3
        if(fDebugMode){
          std::cout << "[SystToolsEventWeight::produce]     sp->GetSystMetaData().size() (expected) = " << sp->GetSystMetaData().size() << "\n"
                    << "[SystToolsEventWeight::produce]     resp.size() of this syst_resp (produced) = " << resp.size() << std::endl;
        }

        // Below check is not valid if there is a dependent dial
        /*
        if(sp->GetSystMetaData().size()!=resp.size()){
          throw cet::exception{ "SystToolsEventWeight" } 
            << "sp->GetFullyQualifiedName() = " << sp->GetFullyQualifiedName() << "\n"
            << "We expect to have " << sp->GetSystMetaData().size() << " knobs for this SystProvider, but "
                    << resp.size() << " are produced. "
                    << "Probably this particular event is not relevant to this systematic variation.\n";
        }
        */

        sbn::evwgh::EventWeightMap mcwgh;

        for( auto const& r: resp){

          // responses.size() is the number of universes
          systtools::SystParamHeader const &sph = fParamHeaderHelper.GetHeader( r.pid );
          std::string prettyName = sph.prettyName;

          if(sph.isResponselessParam) continue;

          if(fDebugMode){
            std::cout << "[SystToolsEventWeight::produce]       pid of this resp = " << r.pid << "\n"
                      << "[SystToolsEventWeight::produce]       prettyName of this resp = " << prettyName << "\n"
                      << "[SystToolsEventWeight::produce]       paramVariations.size() of this resp (expected) = " << sph.paramVariations.size() << "\n"
                      << "[SystToolsEventWeight::produce]       responses.size() of this resp (produced) = " << r.responses.size() << std::endl;
          }
          if(sph.paramVariations.size()!=r.responses.size()){
            throw cet::exception{ "SystToolsEventWeight" }
              << "prettyName of this resp = " << prettyName << "\n"
              << "We expect to have " << sph.paramVariations.size() << " universes, but "
                      << r.responses.size() << " are produced. "
                      << "Probably this particular event is not relevant to this systematic variation.\n";
          }

          // r.responses : std::vector<double>
          // std::map<std::string, std::vector<float> > EventWeightMap
          mcwgh[sp->GetFullyQualifiedName()+"_"+prettyName].assign(r.responses.cbegin(), r.responses.cend());

        } // END loop over knobs

        mcwghvec->push_back(std::move(mcwgh));
        art::Ptr<sbn::evwgh::EventWeightMap> wghPtr = makeWeightPtr(mcwghvec->size() - 1);
        wghassns->addSingle(mclist.at(nMCTruthIndex), wghPtr);

        nMCTruthIndex++;

      } // END loop over MCTruth particles

    } // END systprovider loop

  }

  e.put(std::move(mcwghvec));
  e.put(std::move(wghassns));

}


void SystToolsEventWeight::beginRun(art::Run& run) {

  auto p = std::make_unique<std::vector<sbn::evwgh::EventWeightParameterSet> >();
  for( auto &sp : fParamHeaderHelper.GetSystProvider() ) {
    MF_LOG_INFO("SystToolsEventWeight") << "sp->GetToolType() = " << sp->GetToolType() << "\n"
                                  << "sp->GetFullyQualifiedName() = " << sp->GetFullyQualifiedName() << "\n"
                                  << "sp->GetInstanceName() = " << sp->GetInstanceName() << "\n"
                                  << "Printing each SystParamHeader of this ISystProviderTool..";
    // Note: typedef std::vector<SystParamHeader> SystMetaData;
    auto const& smd = sp->GetSystMetaData();

    // make a map of responsless-response params
    std::map<systtools::paramId_t, std::vector<systtools::paramId_t>> map_resp_to_respless;
    for( auto &sph : smd ){
      if(sph.isResponselessParam){
        auto it = map_resp_to_respless.find( sph.responseParamId );
        if( it != map_resp_to_respless.end() ){
          it->second.push_back( sph.systParamId );
        }
        else{
          map_resp_to_respless[sph.responseParamId] = {};
          map_resp_to_respless[sph.responseParamId].push_back( sph.systParamId );
        }
      }
    }
    if (fDebugMode) {
      for(const auto& it: map_resp_to_respless){
        const auto& sph = systtools::GetParam(smd, it.first);
        std::cout << "Found a dependent dial: " << sph.prettyName << std::endl;
        for(const auto& depdialid: it.second){
          const auto& sph_dep = systtools::GetParam(smd, depdialid);
          std::cout << "- dep dial: " << sph_dep.prettyName << std::endl;
        }
      }
    }

    for( auto &sph : smd ){

      // responsless
      if(sph.isResponselessParam){
        if (fDebugMode) {
          std::cout << "Responsless dial found: " << sph.prettyName << ", thus skipping" << std::endl;
        }
        continue;
      }

      sbn::evwgh::EventWeightParameterSet fParameterSet;
      std::string rwmode = "";

      auto it = map_resp_to_respless.find( sph.systParamId );
      if(it!=map_resp_to_respless.end()){

        for(const auto depdialid: it->second){
          const auto& sph_dep = systtools::GetParam(smd, depdialid);

          if(rwmode=="") rwmode = sph_dep.isRandomlyThrown ? "multisim" : "multisigma";
          else{
            if(rwmode!= (sph_dep.isRandomlyThrown ? "multisim" : "multisigma")){
              throw cet::exception{ "SystToolsEventWeight" }
                << sph.prettyName << " depends on other dials, but the remode are different between the deps dials\n";
            }
          }

          std::vector<float> widths_dep { sph_dep.paramVariations.begin(), sph_dep.paramVariations.end() };
          fParameterSet.AddParameter(sph_dep.prettyName, std::move(widths_dep));

        }



      }
      else{

        rwmode = sph.isRandomlyThrown ? "multisim" : "multisigma";

        std::vector<float> widths { sph.paramVariations.begin(), sph.paramVariations.end() };

        fParameterSet.AddParameter(sph.prettyName, std::move(widths));

      }

      if(rwmode==""){
        throw cet::exception{ "SystToolsEventWeight" }
          << "rwmode not set for " << sph.prettyName << "\n";
      }

      if (fDebugMode) {
        std::cout << "  sph.prettyName = " << sph.prettyName << ", rwmode = " << rwmode << " is added to the header" << std::endl;
      }

      fParameterSet.Configure(sp->GetFullyQualifiedName()+"_"+sph.prettyName, rwmode, sph.paramVariations.size());
      fParameterSet.FillKnobValues();

      p->push_back( std::move(fParameterSet) );

    }
  }

  run.put(std::move(p));

}

}  // namespace sbn::evwgh

DEFINE_ART_MODULE(sbn::evwgh::SystToolsEventWeight)

